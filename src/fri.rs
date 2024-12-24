use crate::merkle::{merkelize, hash, verify_branch, get_branch};
use crate::fft::{fft, inv_fft, get_initial_domain_of_size, halve_domain, get_single_domain_value, halve_single_domain_value};
use crate::circle::{div, scalar_division, scalar_multiply, subtract, CircleImpl, CirclePoint, MODULUS};
use std::ops::Add;
use crate::utils::{log2};
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

pub fn folded_reverse_bit_order(values: &Vec<CirclePoint>) -> Vec<CirclePoint>{
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
        
        let domain_tmp :Vec<CirclePoint> = domain.iter().cloned().step_by(2).collect();
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
        let hash_output = hash(hash_input);
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
    let o = fft(&folded_reverse_bit_order(evaluations), Some(&folded_reverse_bit_order(domain)));
    let zero = CirclePoint::zero();
    return o[halflen..].iter().all(|&c| c == zero );
}

//need to_bytes and from_bytes functions for CirclePoint

pub fn chunkify(values: &Vec<CirclePoint>) -> Vec<Vec<u8>>{
    let mut v = Vec::new();
    for i in 0..values.len(){
        let element = values[i..i + (FOLD_SIZE_RATIO as usize)]
        .iter()
        .flat_map(|&x| x.to_bytes())
        .collect();
        v.push(element);
    }
    v
}

pub fn unchunkify(field: &FieldElement<M31>, data: &[u8]) -> Vec<FieldElement<M31>>
{
    data.chunks(16).map(|chunk| field.from_bytes(chunk)).collect()
}

pub struct Proof {
    // trees : Vec<Vec<Option<Vec<u8>>>>
    roots: Vec<Vec<u8>>,
    branches: Vec<Vec<Vec<Vec<u8>>>>,
    leaf_values: Vec<Vec<Vec<CirclePoint>>>,
    final_values: Vec<CirclePoint>
}

pub fn prove_low_degree(evaluations: &Vec<CirclePoint>) -> Proof {
    let domain = folded_reverse_bit_order(&get_initial_domain_of_size(MODULUS, evaluations.len()));
    let values = folded_reverse_bit_order(evaluations);
    let leaves = Vec::new();
    let trees = Vec::new();
    let roots = Vec::new();
    let rounds = log2((evaluations.len() as u32)/(BASE_CASE_SIZE as u32)) /FOLDS_PER_ROUND as usize;
    print!("generating proof");
    for i in 0..rounds{
        leaves.push(values);
        trees.push(merkelize(chunkify(&values)));
        roots.push(trees[trees.len()-1][1].unwrap());
        let fold_factor : u32 = 1; //wasn't able to translate fold_factor = E(get_challenges(b''.join(roots), M, 4))
        (domain,values) = fold(values, fold_factor, domain);
    }
    let entropy= roots
        .iter()
        .chain(values.iter().map(|x| x.to_bytes()))
        .flatten()
        .cloned()
        .collect();
    let challenges = get_challenges(entropy, evaluations.len() as u32 >> FOLDS_PER_ROUND , NUM_CHALLENGES as usize);
    let round_challenges = (0..rounds).map(|i| {
        challenges.iter().map(|&c| c >> ((i as u32)*FOLDS_PER_ROUND)).collect()}).collect();
    let mut branches = Vec::new();
    for (r_challenges, tree) in round_challenges.iter().zip(trees.iter()) {
        let mut round_branches = Vec::new();
        for &c in r_challenges {
            round_branches.push(get_branch(tree, c));
        }
        branches.push(round_branches);
    }
    let mut leaf_values = Vec::new();
    for (leaf_list, r_challenges) in leaves.iter().zip(round_challenges.iter()) {
        let mut round_leaf_values = Vec::new();
        for &c in r_challenges {
            let start = c * FOLD_SIZE_RATIO;
            let end = (c + 1) * FOLD_SIZE_RATIO;
            if start < leaf_list.len() && end <= leaf_list.len() {
                round_leaf_values.push(leaf_list[start..end].to_vec()); 
            } else {
                round_leaf_values.push(Vec::new());
            }
        }
        leaf_values.push(round_leaf_values);
    }
    Proof { roots, branches, leaf_values, final_values: values }
}

pub fn verify_low_degree(proof: Proof) -> bool{
    let roots = proof.roots;
    let branches = proof.branches;
    let leaf_values = proof.leaf_values;
    let final_values = proof.final_values;
    let M = MODULUS;
    let len_evaluations = final_values.len() << (FOLDS_PER_ROUND as usize * roots.len());

    let mut entropy = Vec::new();
    for root in roots {
        entropy.extend(root);
    }
    for &x in &final_values {
        entropy.extend(x.to_bytes());
    }
    let challenges = get_challenges(&entropy, len_evaluations as u32>> FOLDS_PER_ROUND, NUM_CHALLENGES as usize);

    for i in 0..roots.len(){
        let fold_factor : u32 = 1; //wasn't able to translate fold_factor = E(get_challenges(b''.join(roots), M, 4))
        let evaluation_size = len_evaluations >> (i * FOLDS_PER_ROUND as usize);
        let positions: Vec<usize> = challenges.iter().flat_map(|&c| (c * FOLD_SIZE_RATIO..(c + 1) * FOLD_SIZE_RATIO).map(|x| x as usize)).collect();
        let mut domain = Vec::new();
        if i==0{
            let domain_tmp = Vec::new();
            for pos in positions{
                domain_tmp.push(get_single_domain_value(M, evaluation_size, rbo_index_to_original(evaluation_size, pos)));
            }
            domain = domain_tmp;
        } else {
            let domain_tmp = Vec::new();
            for pos in positions{
                domain_tmp.push(halve_single_domain_value(&get_single_domain_value(M, evaluation_size*2, rbo_index_to_original(evaluation_size*2, pos*2))));
            }
            domain = domain_tmp;
        }
        let flattened_values = leaf_values[i].iter().flat_map(|v| v.iter().cloned()).collect();
        let (folded_values, _) = fold(flattened_values, fold_factor, domain);
        let mut expected_values = Vec::new();
        if i < roots.len()-1{
            let mut expected_values_list = Vec::new();
            for (j, &c) in challenges.iter().enumerate(){
                let value = leaf_values[i+1][j][(c % FOLD_SIZE_RATIO) as usize];
                expected_values_list.push(value);
            }
            expected_values = expected_values_list;
        } else {
            let mut expected_values_list = Vec::new();
            for c in challenges{
                expected_values_list.push(final_values[c as usize]);
            }
            expected_values = expected_values_list;
        }
        assert!(folded_values == expected_values);
        for (j, &c) in challenges.iter().enumerate(){
            assert!(verify_branch(roots[i], c as i32, chunkify(&leaf_values[i][j].to_vec())[0], branches[i][j]));
        }
        let mut challenges_new = Vec::new();
        for c in challenges{
            challenges_new.push(c >> FOLDS_PER_ROUND);
        }
        challenges = challenges_new;
    }
    let final_domain = folded_reverse_bit_order(&halve_domain(&get_initial_domain_of_size(M, final_values.len()*2), true));
    assert!(is_rbo_low_degree(&final_values, &final_domain));
    return true;
}