use crate::circle::{CircleImpl, CirclePoint, G, MODULUS};
use crate::utils::log2;
use crate::utils::{folded_reverse_bit_order, reverse_bit_order};
use lambdaworks_math::field::element::FieldElement;
use lambdaworks_math::field::fields::mersenne31::field::Mersenne31Field as M31;

const TOP_DOMAIN_SIZE: u32 = 1 << 24; // 2**24
const LOG2_TOP_DOMAIN_SIZE: usize = 24;

//  Generator point
pub fn generator_point(G: CirclePoint) -> CirclePoint {
    let mut ans: CirclePoint = G;
    for _ in LOG2_TOP_DOMAIN_SIZE..log2(MODULUS + 1) - 1 {
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
        let doubled: Vec<CirclePoint> = sub_domains[source_start..source_end]
            .iter()
            .map(|x| x.double())
            .collect();
        // .collect::Vec<_>();

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

pub fn compute_bit_orders() -> (Vec<u32>, Vec<u32>) {
    // Initialize arrays
    let mut rbos = vec![0u32; TOP_DOMAIN_SIZE as usize * 2];
    let mut folded_rbos = vec![0u32; TOP_DOMAIN_SIZE as usize * 2];

    // Compute RBO (Reversed Bit Order)
    for i in 0..LOG2_TOP_DOMAIN_SIZE {
        let start = 1 << i;
        let end = 1 << (i + 1);

        let rbo_slice = reverse_bit_order(start as u32, i + 1);
        rbos[start..end].copy_from_slice(&rbo_slice);
    }

    // Compute Folded RBO
    for i in 0..LOG2_TOP_DOMAIN_SIZE {
        let start = 1 << i;
        let end = 1 << (i + 1);

        let folded_rbo_slice = folded_reverse_bit_order(start as u32, i + 1);
        folded_rbos[start..end].copy_from_slice(&folded_rbo_slice);
    }

    (rbos, folded_rbos)
}
