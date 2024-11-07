use num_bigint::{BigUint, RandBigInt};

// output => n^exponet mod p
pub fn exponetiate(n: &BigUint, exponent: &BigUint, p: &BigUint) -> BigUint {
    n.modpow(exponent, p)
}

// output => s = k - c * x mod q
pub fn solve(k: &BigUint, c:&BigUint, x:&BigUint, q:&BigUint) -> BigUint {
    if *k >= c * x {
        return (k - c * x).modpow(&BigUint::from(1u32), q);
    }

    return q - (c * x - k).modpow(&BigUint::from(1u32), q);
}

// cond1: r1 = alpha^s * y1^c mod p
// cond2: r2 = beta^s * y2^c mod p
pub fn verify(r1: &BigUint, r2:&BigUint, alpha: &BigUint, beta:&BigUint, y1:&BigUint, y2:&BigUint, s:&BigUint, c:&BigUint, p:&BigUint) -> bool {
    let cond1 = *r1 == (alpha.modpow(s, p) * y1.modpow(c, p)).modpow(&BigUint::from(1u32), p);
    let cond2 = *r2 == (beta.modpow(s, p) * y2.modpow(c, p)).modpow(&BigUint::from(1u32), p);

    cond1 && cond2
}

pub fn generate_random_below(bound: &BigUint) -> BigUint {
    let mut rng = rand::thread_rng();
    rng.gen_biguint_below(bound)
}


#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_example() {
        let alpha = BigUint::from(4u32);
        let beta = BigUint::from(9u32);
        let p = BigUint::from(23u32);
        let q = BigUint::from(11u32);

        let x = BigUint::from(6u32);
        let k = BigUint::from(7u32);

        let c: BigUint = BigUint::from(4u32);

        let y1 = exponetiate(&alpha, &x, &p);
        let y2 = exponetiate(&beta, &x, &p);
        assert_eq!(y1, BigUint::from(2u32));
        assert_eq!(y2, BigUint::from(3u32));

        let r1 = exponetiate(&alpha, &k, &p);
        let r2 = exponetiate(&beta, &k, &p);
        assert_eq!(r1, BigUint::from(8u32));
        assert_eq!(r2, BigUint::from(4u32));

        let s = solve(&k, &c, &x, &q);
        assert_eq!(s, BigUint::from(5u32));

        let result = verify(&r1, &r2, &alpha, &beta, &y1, &y2, &s, &c, &p);
        assert!(result);

        // Fake Secret
        let x_fake = BigUint::from(7u32);
        let s_fake = solve(&k, &c, &x_fake, &q);

        let result = verify(&r1, &r2, &alpha, &beta, &y1, &y2, &s_fake, &c, &p);
        assert!(!result);

    }

    #[test]
    fn test_example_with_random_numbers() {
        let alpha = BigUint::from(4u32);
        let beta = BigUint::from(9u32);
        let p = BigUint::from(23u32);
        let q = BigUint::from(11u32);

        let x = BigUint::from(6u32);
        let k = generate_random_below(&q);

        let c: BigUint = generate_random_below(&q);

        let y1 = exponetiate(&alpha, &x, &p);
        let y2 = exponetiate(&beta, &x, &p);
        assert_eq!(y1, BigUint::from(2u32));
        assert_eq!(y2, BigUint::from(3u32));

        let r1 = exponetiate(&alpha, &k, &p);
        let r2 = exponetiate(&beta, &k, &p);

        let s = solve(&k, &c, &x, &q);

        let result = verify(&r1, &r2, &alpha, &beta, &y1, &y2, &s, &c, &p);
        assert!(result);
    }
}