use ndarray::{Array, Dimension};

const BASE_CASE_SIZE : u32= 64;
const FOLDS_PER_ROUND : u32 = 3;
const BASE : u32 = 2;
const FOLD_SIZE_RATIO : u32= BASE.pow(FOLDS_PER_ROUND);
const NUM_CHALLENGES : u32 = 80;

pub fn fold<D>(values : &Array<u32, D>, coeff: u32, first_round: bool) where D : Dimension, {

    for i in 0..FOLDS_PER_ROUND{
        let full_len = values.shape()[1];
        let half_len = full_len/2;
    }
}
