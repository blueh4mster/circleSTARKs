use lambdaworks_math::field::element::FieldElement;
use lambdaworks_math::field::fields::mersenne31::field::Mersenne31Field as M31;
//use ndarray at last maybe
use crate::circle::{CircleImpl, CirclePoint};
// use crate::fast_fft::bary_eval;
use crate::precomputes::get_subdomains;

pub fn line_function(
    p1: CirclePoint,
    p2: CirclePoint,
    domain: &Vec<CirclePoint>,
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

pub fn public_args_to_vanish_and_interp(
    domain_size: u32, // might be usize
    indices: &[u32],
    vals: &[Vec<FieldElement<M31>>],
    out_domain: Option<CirclePoint>,
) -> (Vec<FieldElement<M31>>, Vec<Vec<FieldElement<M31>>>) {
    assert!(indices.len() % 2 == 0, "Indices array has odd length");
    let next_power_of_2 = 2usize.pow((indices.len() - 1).ilog2() + 1) as u32;
    assert!(next_power_of_2 < domain_size, "domain size is large");

    let mut lines: Vec<FieldElement<M31>> = Vec::new();

    let sub_domains = get_subdomains();
    let mut eval_domain: Vec<CirclePoint> = sub_domains
        .iter()
        .skip(next_power_of_2 as usize)
        .take(next_power_of_2 as usize)
        .cloned()
        .collect();
    if let Some(out_domain) = out_domain {
        eval_domain.push(out_domain);
    }

    // Initialize vpoly
    let mut vpoly: Vec<FieldElement<M31>> = vec![FieldElement::new(1); eval_domain.len()];

    let points: Vec<CirclePoint> = indices
        .iter()
        .map(|i| sub_domains[(domain_size + i) as usize])
        .collect();

    for i in (0..indices.len()).step_by(2) {
        let mut line = line_function(points[i], points[i + 1], &eval_domain); // this might be that we only need a single line , oonstead of all the lines
        lines.append(&mut line);

        vpoly = vpoly
            .iter()
            .zip(lines.clone())
            .map(|(v, l)| FieldElement::<M31>::new(v.to_raw() * l.to_raw()))
            .collect();
    }

    // Initialize interpolation
    let mut interp = vec![vec![FieldElement::new(0); vals[0].len()]; eval_domain.len()];
    for i in (0..indices.len()).step_by(2) {
        let mut vpoly_adjusted = Vec::with_capacity(eval_domain.len());
        for j in 0..eval_domain.len() {
            // Divide vpoly by corresponding line (simplified)
            let divided_val = vpoly[j] / (lines[i / 2]);
            vpoly_adjusted.push(divided_val);
        }

        let y1 = FieldElement::new(1); //bary_eval
        let y2 = FieldElement::new(1); // bary_eval needs to be replaced

        // Compute interpolant
        let interpolant_vals = interpolant(
            points[i],
            points[i + 1],
            &[vec![vals[i][0] / (&y1)]],
            &[vec![vals[i + 1][0] / (&y2)]],
            &eval_domain,
        );

        // Update interpolation
        for j in 0..eval_domain.len() {
            for k in 0..vals[0].len() {
                // Multiply vpoly_adjusted by interpolant
                let mul_val = vpoly_adjusted[j] * interpolant_vals[0][j];
                interp[j][k] = interp[j][k] + mul_val;
            }
        }
    }
    if out_domain.is_some() {
        return (
            vpoly
                .iter()
                .skip(next_power_of_2 as usize)
                .cloned()
                .collect(),
            interp
                .iter()
                .skip(next_power_of_2 as usize)
                .cloned()
                .collect(),
        );
    }
    return (vpoly, interp);
}
