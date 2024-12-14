use lambdaworks_math::field::element::FieldElement;
use lambdaworks_math::field::fields::mersenne31::field::Mersenne31Field as M31;

// from zorch.m31 import (
//     M31, ExtendedM31, Point, modulus, zeros_like, Z, G
// )
// use crate::utils::cp
use crate::circle::{CircleImpl, CirclePoint, MODULUS};
// use crate::fast_fft::bary_eval;
use crate::precomputes::get_subdomains;
// from fast_fft import bary_eval

// Unsure of domain's Type
pub fn line_function(
    p1: CirclePoint,
    p2: CirclePoint,
    domain: &[CirclePoint],
) -> Vec<FieldElement<M31>> {
    let a = p2.get_y() - p1.get_y();
    let b = p1.get_x() - p2.get_x();
    let c = p2.get_x() * p1.get_y() - p1.get_x() * p2.get_y();
    let res: Vec<FieldElement<M31>> = domain
        .iter()
        .map(|d| a * d.get_x() + b * d.get_y() + c)
        .collect();
    return res;
}

pub fn interpolant(
    p1: CirclePoint,
    p2: CirclePoint,
    v1: &[Vec<FieldElement<M31>>],
    v2: &[Vec<FieldElement<M31>>],
    domain: &[CirclePoint],
) -> Vec<Vec<FieldElement<M31>>> {
    let mut result = Vec::with_capacity(v1.len());

    for j in 0..v1[0].len() {
        let mut column_result = Vec::with_capacity(domain.len()); // format is different take care

        for i in 0..domain.len() {
            let dx = p2.get_x() - p1.get_x();
            let dy = p2.get_y() - p1.get_y();

            // Compute inverse distance (simplified for demonstration)
            let invdist = (dx * dx + dy * dy).inv().unwrap();

            let dot = (domain[i].get_x() - p1.get_x()) * dx + (domain[i].get_y() - p1.get_y()) * dy;

            let interpolated_val = v1[0][j] + (v2[0][j] - v1[0][j]) * dot * invdist;

            column_result.push(interpolated_val);
        }
        result.push(column_result);
    }

    result
}


pub fn public_args_to_vanish_and_interp(domain_size:u32,indices:&[CirclePoint],vals,)
// def public_args_to_vanish_and_interp(domain_size,
//                                      indices,
//                                      vals,
//                                      out_domain=None):
//     assert len(indices) % 2 == 0
//     next_power_of_2 = 2**(len(indices)-1).bit_length() * 2
//     assert next_power_of_2 < domain_size
//     lines = []
//     eval_domain = sub_domains[next_power_of_2: next_power_of_2*2]
//     if out_domain is not None:
//         eval_domain = Point(
//             M31.append(eval_domain.x, out_domain.x),
//             M31.append(eval_domain.y, out_domain.y)
//         )
//     vpoly = M31.zeros(eval_domain.shape) + 1
//     points = sub_domains[domain_size + cp.array(indices)]
//     for i in range(0, len(indices), 2):
//         lines.append(
//             line_function(points[i], points[i+1], eval_domain)
//         )
//         vpoly = vpoly * lines[-1]
//     interp = vpoly.__class__.zeros((eval_domain.shape[0],) + vals.shape[1:])
//     for i in range(0, len(indices), 2):
//         vpoly_adjusted = (
//             (vpoly / lines[i//2])
//             .reshape((eval_domain.shape[0],) + (1,) * (len(vals.shape) - 1))
//         )
//         y1 = bary_eval(vpoly_adjusted[:next_power_of_2], points[i])
//         y2 = bary_eval(vpoly_adjusted[:next_power_of_2], points[i+1])
//         I = interpolant(
//             points[i], vals[i] / y1,
//             points[i+1], vals[i+1] / y2,
//             eval_domain,
//         )
//         interp += vpoly_adjusted * I
//     if out_domain is not None:
//         return vpoly[next_power_of_2:], interp[next_power_of_2:]
//     else:
//         return vpoly, interp
