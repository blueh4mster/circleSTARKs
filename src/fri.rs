use crate::merkle::{merkelize, hash, verify_branch, get_branch};
use crate::fft::{fft, inv_fft, get_initial_domain_of_size, halve_domain, get_single_domain_value, halve_single_domain_value};
use crate::circle::{div, scalar_division, scalar_multiply, subtract, CircleImpl, CirclePoint, MODULUS};
use crate::utils::folded_reverse_bit_order;
use std::ops::Add;
use lambdaworks_math::field::element::FieldElement;
use lambdaworks_math::field::fields::mersenne31::field::Mersenne31Field as M31;

const BASE_CASE_SIZE : u32= 128;
const FOLDS_PER_ROUND : u32 = 3;
const BASE : u32 = 2;
const FOLD_SIZE_RATIO : u32= BASE.pow(FOLDS_PER_ROUND);
const NUM_CHALLENGES : u32 = 80;

pub fn extend_trace(field: u32, trace: &Vec<CirclePoint>) -> Vec<CirclePoint>{
    let small_domain = get_initial_domain_of_size(field, trace.len());
    let coeffs = fft(trace, Some(&small_domain));
    let big_domain = get_initial_domain_of_size(field, trace.len()*2);
    let res = inv_fft(trace, Some(&big_domain));
    res
}

pub fn line_function(P1: CirclePoint, P2: CirclePoint, domain: &[CirclePoint]) ->Vec<FieldElement<M31>>{
    let x1 = P1.get_x();
    let x2 = P2.get_x();
    let y1 = P1.get_y();
    let y2 = P2.get_y();
    let denominator = x1*y2 - x2*y1;
    let a = (y2-y1) / denominator;
    let b = (x1-x2) / denominator;
    let negative : FieldElement<M31> = FieldElement::new(MODULUS-1);
    let c = negative * (a*x1 + b*y1);
    domain.iter().map(|cpoint| a*(cpoint.get_x()) + b*cpoint.get_y() + c).collect()
}

pub fn rbo(values: &Vec<CirclePoint>) -> Vec<CirclePoint> {
    if values.len() == 1{
        return values.to_vec();
    }
    let x : Vec<CirclePoint> = values.iter().cloned().step_by(2).collect();
    let y : Vec<CirclePoint> = values.iter().cloned().skip(1).step_by(2).collect();
    let L = rbo(&x);
    let R = rbo(&y);
    let mut res = L;
    res.extend(R);
    return res;
}

pub fn fold_reverse_bit_order(values: &Vec<CirclePoint>) -> Vec<CirclePoint>{
    let l : Vec<CirclePoint>= values.iter().cloned().step_by(2).collect();
    let r : Vec<CirclePoint>= values.iter().cloned().skip(1).step_by(2).rev().collect();
    let L = rbo(&l);
    let R = rbo(&r);
    let mut o = vec![None; values.len()];
    let mut i = 0;
    for x in L{
        o[i] = Some(x);
        i+=2;
    }
    let mut j = 1;
    for x in R{
        o[j] = Some(x);
        j+=2;
    }
    let res: Vec<CirclePoint> = o.into_iter().filter_map(|u| u).collect();
    res
}

pub fn rbo_index_to_original(length: usize, index: usize) -> usize {
    let sub = format!("{:b}", length + index)[3..(format!("{:b}", length + index).len() - 1)]
        .chars()
        .rev()
        .collect::<String>();
    let sub2 = usize::from_str_radix(sub.as_str(), 2);

    if index % 2 == 0{
        return sub2.unwrap();
    } else {
        return length - 1 - (sub2.unwrap())*2;
    }
}

pub fn fold(mut values: Vec<CirclePoint>, coeff: u32, mut domain: Vec<CirclePoint>) -> (Vec<CirclePoint>, Vec<CirclePoint>) {
    for i in 0..FOLDS_PER_ROUND{
        let left : Vec<CirclePoint>= values.iter().cloned().step_by(2).collect();
        let right : Vec<CirclePoint>= values.iter().cloned().skip(1).step_by(2).collect();
        let mut f0: Vec<Option<CirclePoint>> = vec![None; values.len()/2];
        for j in 0..(values.len()/2) {
            let L = left[j];
            let R = right[j];
            // instead of L+R, use circle point addition
            f0[j] = Some(scalar_division(L+R, 2));
        }
        
        let mut domain_tmp :Vec<CirclePoint> = domain.iter().cloned().step_by(2).collect();
        let mut f1 : Vec<Option<CirclePoint>> = vec![None; f0.len()];
        for k in 0..(values.len()/2){
            let L = left[k];
            let R = right[k];
            let x = domain_tmp[k];
                // not sure about the type of x, might have to use scalar multiply 
            f1[k] = Some(div(subtract(L, R),scalar_multiply(x, 2)));
        }
        let mut vals = vec![None;f0.len()];
        for i in 0..f0.len(){
            let f0val = f0[i].unwrap();
            let f1val = f1[i].unwrap();
            vals[i] = Some(f0val.add(scalar_multiply(f1val, coeff)));
        }
        let vals2 : Vec<CirclePoint> = vals.into_iter().filter_map(|e| e).collect();
        values = vals2;
        let domain2= halve_domain(&domain_tmp, true);
        domain = domain2;
    }
    return (values, domain);
}

pub fn get_challenges(root: &[u8], domain_size: u32, num_challenges: usize) -> Vec<u32>{
    let mut challenge_data = Vec::new();
    for i in 0..((num_challenges + 7) / 8) {
        let mut hash_input = Vec::new();
        hash_input.extend_from_slice(root);
        hash_input.push((i / 256) as u8);
        hash_input.push((i % 256) as u8);
        let hash_output = hash(&hash_input);
        challenge_data.extend_from_slice(&hash_output);
    }

    (0..num_challenges)
        .map(|i| {
            let start = i * 4;
            let end = start + 4;
            let value = u32::from_le_bytes(
                challenge_data[start..end]
                    .try_into()
                    .expect("Invalid slice length"),
            );
            value % domain_size
        }).collect()
}

pub fn is_rbo_low_degree(evaluations: &Vec<CirclePoint>, domain: &Vec<CirclePoint>) -> bool{
    let halflen = evaluations.len()/2;
    let o = fft(&fold_reverse_bit_order(evaluations), Some(&fold_reverse_bit_order(domain)));
    let zero = FieldElement::<M31>::zero();
    return o[halflen..].iter().all(|&c| c.get_x() == zero && c.get_y() == zero );
}

//need to_bytes and from_bytes functions for CirclePoint

pub fn chunkify(values: &Vec<CirclePoint>) -> Vec<&[u8]>{
    let mut v = Vec::new();
    for i in 0..values.len(){
        let element = values[i..i + (FOLD_SIZE_RATIO as usize)]
        .iter()
        .flat_map(|&x| x.to_le_bytes())
        .collect();
        v.push(element);
    }
    v
}

pub fn unchunkify(field: &FieldElement<M31>, data: &[u8]) -> Vec<FieldElement<M31>>
{
    data.chunks(16).map(|chunk| field.from_bytes(chunk)).collect()
}