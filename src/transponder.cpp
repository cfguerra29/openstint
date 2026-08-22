#include "transponder.hpp"

#include <bit>
#include <string>

extern "C" {
#include <fec.h>
}
#include <liquid/liquid.h>

static void *viterbi_decoder;
static crc_scheme crc8_scheme = LIQUID_CRC_8;

void init_transponders() {
    viterbi_decoder = create_viterbi29(32);
}

int decode_openstint(const uint8_t *softbits, uint32_t *transponder_id) {
    uint8_t decoded[4];

    init_viterbi29(viterbi_decoder, 0);
    update_viterbi29_blk(viterbi_decoder, const_cast<uint8_t*>(softbits), 32+8); // khm...
    chainback_viterbi29(viterbi_decoder, decoded, 32, 0);
    
    *transponder_id = (static_cast<uint32_t>(decoded[0]) << 16) | (static_cast<uint32_t>(decoded[1]) << 8) | static_cast<uint32_t>(decoded[2]);
    return crc_validate_message(crc8_scheme, decoded, 3, decoded[3]);
}

int decode_rc3(const uint8_t *softbits, uint32_t *transponder_id, uint8_t *status_code) {
    // RC3 use a K=24, r=1/2 convolutional encoder with polynoms 0xEEC20F and 0xEEC20D
    // Decoding this properly with error correction must have some unknown trick. However,
    // we can do non-trivial decoding as well.
    // Note, the generating polynoms differ only in one bit (0xf=1111 vs 0xd=1101). 
    // This reduces the complexity to:
    // bit0 = parity(SHREG & 0xEEC20C) ^ SHREG[1] ^ SHREG[0]
    // bit1 = parity(SHREG & 0xEEC20C) ^            SHREG[0]
    // As such, we know: 
    // bit0 ^ bit1 = SHREG[1]
    // SHREG[0] = bit0 ^ parity(SHREG & 0xEEC20C) ^ SHREG[1]
    // SHREG[0] = bit1 ^ parity(SHREG & 0xEEC20C)
    // we can use these for a makeshift error correction
    
    uint64_t shreg = 0; // shift register
    bool last_ok = true; // last 2-bits were successfully decoded
    
    // Before encoding, the 24 bit transponder id is scrambled with extra 8 bits,
    // resulting in 32 bits. This is further appended with 0x00 so viterbi-decoder (???)
    // can process it. The rate=1/2 encoder generates 2x40=80 bits in total.
    int sym = 0, prev_sym = 0;  // differential decoder
    for (int i=0; i<80; i+=2) {
        int p = std::popcount(shreg & 0xEEC20C) % 2; // parity bit from SHREG

        // differential-BPSK is decoded here as well (sym ^ prev_sym magic):
        sym = (softbits[i+0] > 127) ? 1 : 0;
        int b0 = sym ^ prev_sym;
        prev_sym = (softbits[i+1] > 127) ? 1 : 0;
        int b1 = prev_sym ^ sym; // decode softbit bit1

        int shreg1 = (shreg & 2) ? 1 : 0; // shift register last-1 bit
        // two estimates for SHREG[0] (should be equal):
        int shreg0p0 = p ^ shreg1 ^ b0;
        int shreg0p1 = p ^ b1;
        if (last_ok) { // no error correction for SHREG[1] is needed
            last_ok = (shreg0p0 == shreg0p1);
            if (last_ok) {
                shreg |= shreg0p0; // high certainty, write bit to SHREG
            }
        } else { // must correct SHREG[1] based on bit0 and bit1
            int shreg1p = b0 ^ b1; // SHREG[1] guesstimate; see top comment
            // SHREG[0]'  = b0 ^ PAR(...) ^ SHREG[1] = b0 ^ PAR(...) ^ (b0 ^ b1) = PAR(...) ^ b1 = shreg0p1
            shreg |= (shreg1p << 1) | shreg0p1;
            last_ok = true;
        }
        shreg <<= 1; // no matter if we have the last bit correctly, shift it
    }

    // error detection
    shreg >>= 1;
    uint32_t trail = static_cast<uint32_t>(shreg & 0xff);
    uint32_t message = static_cast<uint32_t>((shreg>>8) & 0xffffffff);

    // Example: your transponder id is "1234567", or "0b00010010_11010110_10000111".
    // To get the pre-encoded message, reverse the binary word and break into 
    // chunks of 3. You'll end up with 24/3=8 chunks:
    // 111 000 010 110 101 101 001 000
    // Suffix each chunk with a bit from a "status code",
    // in my RC4 hybrid transponder, it was 00000101
    // 1110 0000 0100 1100 1010 1011 0010 0001
    uint32_t tid = 0; // transponder_id
    uint8_t status = 0; // status code
    for (int i=0; i<32; i++) {
        uint32_t bitmask = (1 << i);
        uint32_t bit = (message & bitmask) ? 1 : 0;

        if (i % 4 != 0) { // every 4th is status bit
            tid = (tid << 1) | bit;
        } else {
            status = (status << 1) | bit;
        }
    }
    *transponder_id = tid;
    *status_code = status;

    // the last byte must be zero (tail==0 error check)
    return (trail == 0);
}

// Vostok transponders reuse the RC3 preamble but carry a completely different
// payload: no convolutional code at all, just the plain 24 bit id sent twice
// with running XOR checksums. Frame layout of the 80 payload symbols, after
// differential (DBPSK) demodulation, as 10 bytes:
//
//   byte 0-1  header (0x631A on the unit this was reverse engineered from)
//   byte 2-4  transponder id, 24 bit big-endian
//   byte 5    running XOR of bytes 0..4
//   byte 6-8  transponder id again
//   byte 9    running XOR of bytes 0..8
//
// Verified against a live capture: 15 consecutive frames of a Vostok labelled
// 5112266 decoded bit-identical, id and both checksums matching.
//
// The id must appear twice and both checksums must hold: 40 bits of constraint,
// so a false positive is ~2^-40. That is far stronger than the 8 bit tail check
// RC3 relies on, which is why this may safely run as a fallback after decode_rc3.
int decode_vostok(const uint8_t *softbits, uint32_t *transponder_id) {
    uint8_t bytes[10] = {0};
    int prev_sym = 0;
    for (int i = 0; i < 80; i++) {
        const int sym = (softbits[i] > 127) ? 1 : 0;
        const int bit = sym ^ prev_sym;   // differential decode
        prev_sym = sym;
        bytes[i / 8] = static_cast<uint8_t>((bytes[i / 8] << 1) | bit);
    }

    const uint32_t id1 = (static_cast<uint32_t>(bytes[2]) << 16) |
                         (static_cast<uint32_t>(bytes[3]) << 8) | bytes[4];
    const uint32_t id2 = (static_cast<uint32_t>(bytes[6]) << 16) |
                         (static_cast<uint32_t>(bytes[7]) << 8) | bytes[8];
    if (id1 != id2) { return 0; }

    uint8_t xsum = 0;
    for (int i = 0; i < 5; i++) { xsum ^= bytes[i]; }
    if (xsum != bytes[5]) { return 0; }
    for (int i = 5; i < 9; i++) { xsum ^= bytes[i]; }
    if (xsum != bytes[9]) { return 0; }

    if (id1 == 0 || id1 >= 10000000) { return 0; }

    *transponder_id = id1;
    return 1;
}

void AmbRcBlacklist::process(uint64_t timestamp, uint8_t status_code, uint32_t transponder_id) {
    // not a candidate status/validation message:
    if ((status_code & 0xf8) != 0xf8) return; // not an AmbRc message
    if ((status_code & 0x07) == 0) return; // no counter set

    uint8_t msb8 = static_cast<uint8_t>((transponder_id >> 16) & 0xff);
    if (status_code != 0xff) { // status/validation message for sure
        msb8_timestamps[msb8] = timestamp;
    } else { // might be a status/validation message
        auto it = msb8_timestamps.find(msb8);
        if (it != msb8_timestamps.end() && timestamp <= (it->second + 250000ul)) {
            banned_transponders.insert(transponder_id);
        }
    }
}

bool AmbRcBlacklist::check_banned(uint32_t transponder_id) const {
    return banned_transponders.contains(transponder_id);
}
