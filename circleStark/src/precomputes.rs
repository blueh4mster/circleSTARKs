use crate::{circle::CirclePoint, inverse, CircleImpl, G, Z};
use lambdaworks_math::field::element::FieldElement;
use lambdaworks_math::field::fields::mersenne31::field::Mersenne31Field as M31;

const TOP_DOMAIN_SIZE: u32 = 1 << 24; // 2**24

const modulus: u32 = 1; // will be changed

pub fn log2(x: u32) -> usize {
    (x as f64).log2() as usize
}
//  Generator point
pub fn generator_point(G: CirclePoint) -> CirclePoint {
    let mut ans: CirclePoint = G;
    for _ in log2(TOP_DOMAIN_SIZE)..log2(modulus + 1) - 1 {
        ans = ans.double();
    }
    ans
}

// Get the subdomains
pub fn get_subdomains() -> Vec<CirclePoint> {
    let mut top_domain = CirclePoint::zeroes(2);
    top_domain[1] = G();
    for i in 1..log2(TOP_DOMAIN_SIZE * 2) {
        let doubled: Vec<CirclePoint> = top_domain.iter().map(|x| x.double()).collect();

        let mut new_domain = CirclePoint::zeroes(1 << (i + 1));

        // Fill even indices with doubled points
        for (idx, point) in doubled.iter().enumerate() {
            new_domain[idx * 2] = point.clone();
        }

        // Fill odd indices with doubled points + G
        for (idx, point) in doubled.iter().enumerate() {
            new_domain[idx * 2 + 1] = point.clone() + G().clone();
        }

        top_domain = new_domain;
    }

    //  only odd-indexed points
    top_domain.into_iter().step_by(2).collect();

    let mut sub_domains = CirclePoint::zeroes((TOP_DOMAIN_SIZE as usize) * 2);

    // Copy the top domain into the second half of sub_domains
    sub_domains[(TOP_DOMAIN_SIZE as usize)..].clone_from_slice(&top_domain);

    // Iterate backwards through the domain sizes
    for i in (0..=log2(TOP_DOMAIN_SIZE) - 1).rev() {
        let start_idx = 1 << i;
        let end_idx = 1 << (i + 1);
        let source_start = end_idx.clone();
        let source_end = 1 << (i + 3) * 3;

        // Get the source slice and double it
        let doubled = &sub_domains[source_start..source_end]
            .iter()
            .map(|x| x.double())
            .collect();

        // Copy the doubled points to their destination
        sub_domains[start_idx..end_idx].clone_from_slice(&doubled);
    }

    sub_domains
}

pub fn inverse_x(sub_domains: Vec<CirclePoint>) -> Vec<FieldElement<M31>> {
    sub_domains.iter().map(|c| c.inverse_x()).collect()
}

pub fn inverse_y(sub_domains: Vec<CirclePoint>) -> Vec<FieldElement<M31>> {
    sub_domains.iter().map(|c| c.inverse_y()).collect()
}
// invx = 1 / sub_domains.x
// invy = 1 / sub_domains.y

// rbos = cp.zeros(TOP_DOMAIN_SIZE * 2, dtype=cp.uint32)
// for i in range(log2(TOP_DOMAIN_SIZE)):
//     rbos[2**i:2**(i+1)] = reverse_bit_order(cp.arange(2**i))

// folded_rbos = cp.zeros(TOP_DOMAIN_SIZE * 2, dtype=cp.uint32)
// for i in range(log2(TOP_DOMAIN_SIZE)):
//     folded_rbos[2**i:2**(i+1)] = folded_reverse_bit_order(cp.arange(2**i))
