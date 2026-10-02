#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <vector>

namespace baysor {

class Xoshiro256pp {
public:
    using result_type = std::uint64_t;

    explicit Xoshiro256pp(std::uint64_t seed_value = 1) {
        seed(seed_value);
    }

    void seed(std::uint64_t seed_value) {
        std::vector<std::uint32_t> seed_words;
        do {
            seed_words.push_back(static_cast<std::uint32_t>(seed_value & 0xffffffffULL));
            seed_value >>= 32;
        } while (seed_value != 0);

        const auto digest = sha256(
            reinterpret_cast<const unsigned char*>(seed_words.data()),
            seed_words.size() * sizeof(std::uint32_t)
        );

        for (int i = 0; i < 4; ++i) {
            state_[i] = load_le64(digest.data() + 8 * i);
        }

        if ((state_[0] | state_[1] | state_[2] | state_[3]) == 0) {
            state_[0] = 0x9e3779b97f4a7c15ULL; // GCOVR_EXCL_LINE: a SHA-256 digest is never all-zero
        }
    }

    static constexpr result_type min() {
        return std::numeric_limits<result_type>::min();
    }

    static constexpr result_type max() {
        return std::numeric_limits<result_type>::max();
    }

    result_type operator()() {
        return next_uint64();
    }

    double rand_float64() {
        return static_cast<double>(next_uint64() >> 11) * 0x1.0p-53;
    }

private:
    inline static constexpr std::uint32_t kSha256Init[8] = {
        0x6a09e667U, 0xbb67ae85U, 0x3c6ef372U, 0xa54ff53aU,
        0x510e527fU, 0x9b05688cU, 0x1f83d9abU, 0x5be0cd19U
    };

    inline static constexpr std::uint32_t kSha256Round[64] = {
        0x428a2f98U, 0x71374491U, 0xb5c0fbcfU, 0xe9b5dba5U,
        0x3956c25bU, 0x59f111f1U, 0x923f82a4U, 0xab1c5ed5U,
        0xd807aa98U, 0x12835b01U, 0x243185beU, 0x550c7dc3U,
        0x72be5d74U, 0x80deb1feU, 0x9bdc06a7U, 0xc19bf174U,
        0xe49b69c1U, 0xefbe4786U, 0x0fc19dc6U, 0x240ca1ccU,
        0x2de92c6fU, 0x4a7484aaU, 0x5cb0a9dcU, 0x76f988daU,
        0x983e5152U, 0xa831c66dU, 0xb00327c8U, 0xbf597fc7U,
        0xc6e00bf3U, 0xd5a79147U, 0x06ca6351U, 0x14292967U,
        0x27b70a85U, 0x2e1b2138U, 0x4d2c6dfcU, 0x53380d13U,
        0x650a7354U, 0x766a0abbU, 0x81c2c92eU, 0x92722c85U,
        0xa2bfe8a1U, 0xa81a664bU, 0xc24b8b70U, 0xc76c51a3U,
        0xd192e819U, 0xd6990624U, 0xf40e3585U, 0x106aa070U,
        0x19a4c116U, 0x1e376c08U, 0x2748774cU, 0x34b0bcb5U,
        0x391c0cb3U, 0x4ed8aa4aU, 0x5b9cca4fU, 0x682e6ff3U,
        0x748f82eeU, 0x78a5636fU, 0x84c87814U, 0x8cc70208U,
        0x90befffaU, 0xa4506cebU, 0xbef9a3f7U, 0xc67178f2U
    };

    std::uint64_t state_[4] = {0, 0, 0, 0};

    static inline std::uint64_t rotl(const std::uint64_t x, int k) {
        return (x << k) | (x >> (64 - k));
    }

    static inline std::uint32_t rotr32(std::uint32_t x, int k) {
        return (x >> k) | (x << (32 - k));
    }

    static std::uint64_t load_le64(const unsigned char* p) {
        return (static_cast<std::uint64_t>(p[0])      ) |
               (static_cast<std::uint64_t>(p[1]) <<  8) |
               (static_cast<std::uint64_t>(p[2]) << 16) |
               (static_cast<std::uint64_t>(p[3]) << 24) |
               (static_cast<std::uint64_t>(p[4]) << 32) |
               (static_cast<std::uint64_t>(p[5]) << 40) |
               (static_cast<std::uint64_t>(p[6]) << 48) |
               (static_cast<std::uint64_t>(p[7]) << 56);
    }

    static std::uint32_t load_be32(const unsigned char* p) {
        return (static_cast<std::uint32_t>(p[0]) << 24) |
               (static_cast<std::uint32_t>(p[1]) << 16) |
               (static_cast<std::uint32_t>(p[2]) << 8) |
               static_cast<std::uint32_t>(p[3]);
    }

    static void store_be32(std::uint32_t value, unsigned char* out) {
        out[0] = static_cast<unsigned char>(value >> 24);
        out[1] = static_cast<unsigned char>(value >> 16);
        out[2] = static_cast<unsigned char>(value >> 8);
        out[3] = static_cast<unsigned char>(value);
    }

    static std::array<unsigned char, 32> sha256(
        const unsigned char* data, std::size_t len
    ) {
        std::vector<unsigned char> msg(data, data + len);
        msg.push_back(0x80U);
        while ((msg.size() % 64) != 56) {
            msg.push_back(0x00U);
        }

        const std::uint64_t bit_len = static_cast<std::uint64_t>(len) * 8U;
        for (int shift = 56; shift >= 0; shift -= 8) {
            msg.push_back(static_cast<unsigned char>(bit_len >> shift));
        }

        std::uint32_t h[8];
        for (int i = 0; i < 8; ++i) {
            h[i] = kSha256Init[i];
        }

        std::uint32_t w[64];
        for (std::size_t chunk = 0; chunk < msg.size(); chunk += 64) {
            for (int i = 0; i < 16; ++i) {
                w[i] = load_be32(msg.data() + chunk + 4 * i);
            }
            for (int i = 16; i < 64; ++i) {
                const std::uint32_t s0 =
                    rotr32(w[i - 15], 7) ^ rotr32(w[i - 15], 18) ^ (w[i - 15] >> 3);
                const std::uint32_t s1 =
                    rotr32(w[i - 2], 17) ^ rotr32(w[i - 2], 19) ^ (w[i - 2] >> 10);
                w[i] = w[i - 16] + s0 + w[i - 7] + s1;
            }

            std::uint32_t a = h[0];
            std::uint32_t b = h[1];
            std::uint32_t c = h[2];
            std::uint32_t d = h[3];
            std::uint32_t e = h[4];
            std::uint32_t f = h[5];
            std::uint32_t g = h[6];
            std::uint32_t hh = h[7];

            for (int i = 0; i < 64; ++i) {
                const std::uint32_t s1 =
                    rotr32(e, 6) ^ rotr32(e, 11) ^ rotr32(e, 25);
                const std::uint32_t ch = (e & f) ^ ((~e) & g);
                const std::uint32_t temp1 = hh + s1 + ch + kSha256Round[i] + w[i];
                const std::uint32_t s0 =
                    rotr32(a, 2) ^ rotr32(a, 13) ^ rotr32(a, 22);
                const std::uint32_t maj = (a & b) ^ (a & c) ^ (b & c);
                const std::uint32_t temp2 = s0 + maj;

                hh = g;
                g = f;
                f = e;
                e = d + temp1;
                d = c;
                c = b;
                b = a;
                a = temp1 + temp2;
            }

            h[0] += a;
            h[1] += b;
            h[2] += c;
            h[3] += d;
            h[4] += e;
            h[5] += f;
            h[6] += g;
            h[7] += hh;
        }

        std::array<unsigned char, 32> digest{};
        for (int i = 0; i < 8; ++i) {
            store_be32(h[i], digest.data() + 4 * i);
        }
        return digest;
    }

    std::uint64_t next_uint64() {
        const std::uint64_t result = rotl(state_[0] + state_[3], 23) + state_[0];
        const std::uint64_t t = state_[1] << 17;

        state_[2] ^= state_[0];
        state_[3] ^= state_[1];
        state_[1] ^= state_[2];
        state_[0] ^= state_[3];

        state_[2] ^= t;
        state_[3] = rotl(state_[3], 45);

        return result;
    }
};

} // namespace baysor
