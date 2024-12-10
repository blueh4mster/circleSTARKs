use std::ops::Add;

use crate::circle::{scalar_division, scalar_multiply, subtract, CircleImpl, CirclePoint, MODULUS};
use lambdaworks_math::field::element::FieldElement;
use lambdaworks_math::field::fields::mersenne31::field::Mersenne31Field as M31;

// all the functions are specific to M31
fn get_generator(field_modulus: u32) -> CirclePoint {
    for y in 2..field_modulus {
        let y_pt: FieldElement<M31> = FieldElement::new(y);
        let x_pt: FieldElement<M31> = FieldElement::new(1 - y * y).sqrt().unwrap().0; // check docs for once , if its the first value or not
        let mut point = <CirclePoint as CircleImpl>::new_with_field_elements(x_pt, y_pt);
        let point_cloned = point.clone();
        for _ in 0..(f64::log2((field_modulus + 1) as f64) as usize - 1) {
            point = point.double();
        }

        if point != <CirclePoint as CircleImpl>::new(1, 0) {
            return point_cloned;
        }
    }
    panic!("Could not find generator");
}

fn get_initial_domain_of_size(field_modulus: u32, size: usize) -> Vec<CirclePoint> {
    assert!(size < field_modulus as usize);
    let mut g = get_generator(field_modulus);

    for _ in 0..(f64::log2((field_modulus + 1) as f64 / size as f64) as usize - 1) {
        g = g.double();
    }

    let gx2 = g.double();
    let mut domain = vec![g];
    for _ in 1..size {
        let temp_domain_last = domain.last().unwrap();
        let res = gx2.add(*temp_domain_last);
        domain.push(res);
    }
    domain
}

fn get_single_domain_value(field_modulus: u32, size: usize, index: usize) -> CirclePoint {
    assert!(size < field_modulus as usize);
    let mut g = get_generator(field_modulus);

    for _ in 0..(f64::log2((field_modulus + 1) as f64 / size as f64) as usize - 1) {
        g = g.double(); // cause we need to double the field elements , not according to CirclePoint numbers like in (FieldElement<M31>,FieldElement<M31>)
    }

    scalar_multiply(g, (2 * index + 1) as u32)
}

fn halve_domain(domain: &[CirclePoint], preserve_length: bool) -> Vec<FieldElement<M31>> {
    let new_length = if preserve_length {
        domain.len()
    } else {
        domain.len() / 2
    };
    domain
        .iter()
        .take(new_length)
        .map(|c| c.get_x().clone())
        .collect()
}

fn halve_single_domain_value(value: &CirclePoint) -> FieldElement<M31> {
    value.get_x().clone()
}

fn fft(vals: &[CirclePoint], domain: Option<&[CirclePoint]>) -> Vec<CirclePoint> {
    if vals.len() == 1 {
        return vals.to_vec();
    }

    let domain = match domain {
        Some(d) => d.to_vec(),
        None => get_initial_domain_of_size(MODULUS, vals.len()),
    };
    let half_domain = halve_domain(&domain, false); //get the half length

    let (f0, f1) = if let (x, y) = (
        domain.get(0).unwrap().get_x(),
        domain.get(0).unwrap().get_y(),
    ) {
        let (left, right) = vals.split_at(vals.len() / 2);
        let right_reversed: Vec<CirclePoint> = right.iter().rev().cloned().collect();
        let f0: Vec<_> = left
            .iter()
            .zip(&right_reversed)
            .map(|(&l, &r)| (scalar_division(l.add(r), 2)))
            .collect();
        let f1: Vec<_> = left
            .iter()
            .zip(&right_reversed)
            .zip(domain.iter())
            .map(|((&l, &r), d)| scalar_division(subtract(l, r), scalar_multiply(*d, 2)))
            .collect();
        (f0, f1)
    } else {
        let (left, right) = vals.split_at(vals.len() / 2);
        let right_reversed: Vec<_> = right.iter().rev().cloned().collect();
        let f0: Vec<_> = left
            .iter()
            .zip(&right_reversed)
            .map(|(&l, &r)| (l + r) / CirclePoint::from(2.0))
            .collect();
        let f1: Vec<_> = left
            .iter()
            .zip(&right_reversed)
            .zip(domain.iter())
            .map(|((&l, &r), &x)| (l - r) / CirclePoint::from(2.0 * x))
            .collect();
        (f0, f1)
    };

    let mut result = vec![CirclePoint::zero(); domain.len()];
    let half_result = fft(&f0, Some(&half_domain));
    for (i, &val) in half_result.iter().enumerate() {
        result[i * 2] = val;
    }
    let half_result = fft(&f1, Some(&half_domain));
    for (i, &val) in half_result.iter().enumerate() {
        result[i * 2 + 1] = val;
    }

    result
}

fn inv_fft(vals: &[CirclePoint], domain: Option<&[CirclePoint]>) -> Vec<CirclePoint> {
    if vals.len() == 1 {
        return vals.to_vec();
    }

    let domain = match domain {
        Some(d) => d.to_vec(),
        None => get_initial_domain_of_size(MODULUS, vals.len()),
    };
    let half_domain = halve_domain(&domain, false);

    let f0 = inv_fft(
        &vals.iter().step_by(2).cloned().collect::<Vec<_>>(),
        &half_domain,
    );
    let f1 = inv_fft(
        &vals.iter().skip(1).step_by(2).cloned().collect::<Vec<_>>(),
        &half_domain,
    );

    let (left, right): (Vec<_>, Vec<_>) = if let Some((x, y)) = domain.get(0) {
        f0.iter()
            .zip(&f1)
            .zip(domain.iter())
            .map(|((&l, &r), &(_, y))| ((l + y * r), (l - y * r)))
            .unzip()
    } else {
        f0.iter()
            .zip(&f1)
            .zip(domain.iter())
            .map(|((&l, &r), &x)| ((l + x * r), (l - x * r)))
            .unzip()
    };

    let mut result = left;
    result.extend(right.iter().rev().cloned());
    result
}
