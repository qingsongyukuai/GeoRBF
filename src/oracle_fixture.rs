use std::fmt::Write;

pub(crate) fn verify_fixture_identity(
    fixture: &str,
    expected_sha256: &str,
    case_id: &str,
    source: &str,
    output_sha256: &str,
) {
    verify_artifact_identity(fixture, expected_sha256);
    assert!(
        fixture.contains(case_id),
        "fixture must retain its stable CaseId"
    );
    assert!(
        fixture.contains(source),
        "fixture must retain its provenance source"
    );
    assert!(
        fixture.contains(output_sha256),
        "fixture must retain its independently generated output hash"
    );
    assert!(fixture.contains("\"working_decimal_digits\": 120"));
    assert!(fixture.contains("\"rounding\": \"ROUND_HALF_EVEN\""));
}

pub(crate) fn verify_artifact_identity(artifact: &str, expected_sha256: &str) {
    assert_eq!(sha256_hex(artifact.as_bytes()), expected_sha256);
}

pub(crate) fn hex_scalar(fixture: &str, key: &str) -> f64 {
    hex_values(fixture, key, 1)[0]
}

pub(crate) fn hex_values(fixture: &str, key: &str, count: usize) -> Vec<f64> {
    let key_marker = format!("\"{key}\":");
    let start = fixture
        .find(&key_marker)
        .unwrap_or_else(|| panic!("missing fixture key {key}"));
    let mut remaining = &fixture[start + key_marker.len()..];
    let marker = "\"f64_hex\": \"";
    let mut values = Vec::with_capacity(count);
    for _ in 0..count {
        let value_start = remaining
            .find(marker)
            .unwrap_or_else(|| panic!("missing hexadecimal f64 below key {key}"))
            + marker.len();
        remaining = &remaining[value_start..];
        let value_end = remaining
            .find('"')
            .unwrap_or_else(|| panic!("unterminated hexadecimal f64 below key {key}"));
        values.push(parse_hex_f64(&remaining[..value_end]));
        remaining = &remaining[value_end + 1..];
    }
    values
}

fn parse_hex_f64(value: &str) -> f64 {
    let (negative, magnitude) = value
        .strip_prefix('-')
        .map_or((false, value), |magnitude| (true, magnitude));
    if magnitude == "0x0.0p+0" {
        return 0.0;
    }
    let normalized = magnitude
        .strip_prefix("0x1.")
        .unwrap_or_else(|| panic!("unsupported hexadecimal f64 {value}"));
    let (fraction, exponent) = normalized
        .split_once('p')
        .unwrap_or_else(|| panic!("missing exponent in hexadecimal f64 {value}"));
    assert!(fraction.len() <= 13, "f64 fraction is wider than 52 bits");
    let fraction_bits = u64::from_str_radix(fraction, 16)
        .unwrap_or_else(|_| panic!("invalid hexadecimal fraction in {value}"))
        << (4 * (13 - fraction.len()));
    let unbiased_exponent = exponent
        .parse::<i32>()
        .unwrap_or_else(|_| panic!("invalid binary exponent in {value}"));
    assert!((-1022..=1023).contains(&unbiased_exponent));
    let sign = if negative { 1_u64 << 63 } else { 0 };
    let exponent_bits = ((unbiased_exponent + 1023) as u64) << 52;
    f64::from_bits(sign | exponent_bits | fraction_bits)
}

fn sha256_hex(input: &[u8]) -> String {
    const INITIAL: [u32; 8] = [
        0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a, 0x510e527f, 0x9b05688c, 0x1f83d9ab,
        0x5be0cd19,
    ];
    const ROUND: [u32; 64] = [
        0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4,
        0xab1c5ed5, 0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe,
        0x9bdc06a7, 0xc19bf174, 0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f,
        0x4a7484aa, 0x5cb0a9dc, 0x76f988da, 0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7,
        0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967, 0x27b70a85, 0x2e1b2138, 0x4d2c6dfc,
        0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85, 0xa2bfe8a1, 0xa81a664b,
        0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070, 0x19a4c116,
        0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
        0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7,
        0xc67178f2,
    ];

    let bit_length = (input.len() as u64) * 8;
    let mut padded = input.to_vec();
    padded.push(0x80);
    while padded.len() % 64 != 56 {
        padded.push(0);
    }
    padded.extend_from_slice(&bit_length.to_be_bytes());

    let mut state = INITIAL;
    for chunk in padded.chunks_exact(64) {
        let mut words = [0_u32; 64];
        for (index, bytes) in chunk.chunks_exact(4).enumerate() {
            words[index] = u32::from_be_bytes([bytes[0], bytes[1], bytes[2], bytes[3]]);
        }
        for index in 16..64 {
            let s0 = words[index - 15].rotate_right(7)
                ^ words[index - 15].rotate_right(18)
                ^ (words[index - 15] >> 3);
            let s1 = words[index - 2].rotate_right(17)
                ^ words[index - 2].rotate_right(19)
                ^ (words[index - 2] >> 10);
            words[index] = words[index - 16]
                .wrapping_add(s0)
                .wrapping_add(words[index - 7])
                .wrapping_add(s1);
        }

        let [mut a, mut b, mut c, mut d, mut e, mut f, mut g, mut h] = state;
        for index in 0..64 {
            let sum1 = h
                .wrapping_add(e.rotate_right(6) ^ e.rotate_right(11) ^ e.rotate_right(25))
                .wrapping_add((e & f) ^ ((!e) & g))
                .wrapping_add(ROUND[index])
                .wrapping_add(words[index]);
            let sum0 = (a.rotate_right(2) ^ a.rotate_right(13) ^ a.rotate_right(22))
                .wrapping_add((a & b) ^ (a & c) ^ (b & c));
            h = g;
            g = f;
            f = e;
            e = d.wrapping_add(sum1);
            d = c;
            c = b;
            b = a;
            a = sum1.wrapping_add(sum0);
        }
        for (slot, value) in state.iter_mut().zip([a, b, c, d, e, f, g, h]) {
            *slot = slot.wrapping_add(value);
        }
    }

    let mut output = String::with_capacity(64);
    for word in state {
        write!(&mut output, "{word:08x}").expect("writing to a String cannot fail");
    }
    output
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fixture_hash_consumer_matches_the_sha256_standard_vector() {
        assert_eq!(
            sha256_hex(b"abc"),
            "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"
        );
    }
}
