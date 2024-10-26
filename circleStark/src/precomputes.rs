use crate::{circle::CirclePoint, CircleImpl};
use lambdaworks_math::field::fields::mersenne31;

const TOP_DOMAIN_SIZE: u32 = 1 << 24; // 2**24

const modulus: u32 = 1; // will be changed

//  Generator point
pub fn generator_point(G: CirclePoint) -> CirclePoint {
    let mut ans: CirclePoint = G;
    for _ in (TOP_DOMAIN_SIZE as u32).ilog2()..(modulus + 1).ilog2() - 1 {
        ans = ans.double();
    }
    ans
}

// Get the subdomains

pub fn get_subdomains(){
    let top_domain=
}
// def get_subdomains():
//     # Compute the points in the largest-size domain
//     top_domain = Point.zeros(2)
//     top_domain[1] = G
//     for i in range(1, log2(TOP_DOMAIN_SIZE * 2)):
//         doubled = top_domain.double()
//         top_domain = Point.zeros(2**(i+1))
//         top_domain[::2] = doubled
//         top_domain[1::2] = doubled + G
//     top_domain = top_domain[1::2]

//     # Compute an array that contains the top domain and all smaller-size domains
//     sub_domains = Point.zeros(TOP_DOMAIN_SIZE * 2)
//     sub_domains[TOP_DOMAIN_SIZE:] = top_domain
//     for i in range(log2(TOP_DOMAIN_SIZE)-1, -1, -1):
//         sub_domains[2**i:2**(i+1)] = sub_domains[2**(i+1):(2**i)*3].double()
//     return sub_domains

// sub_domains = get_subdomains()

// invx = 1 / sub_domains.x
// invy = 1 / sub_domains.y

// rbos = cp.zeros(TOP_DOMAIN_SIZE * 2, dtype=cp.uint32)
// for i in range(log2(TOP_DOMAIN_SIZE)):
//     rbos[2**i:2**(i+1)] = reverse_bit_order(cp.arange(2**i))

// folded_rbos = cp.zeros(TOP_DOMAIN_SIZE * 2, dtype=cp.uint32)
// for i in range(log2(TOP_DOMAIN_SIZE)):
//     folded_rbos[2**i:2**(i+1)] = folded_reverse_bit_order(cp.arange(2**i))
