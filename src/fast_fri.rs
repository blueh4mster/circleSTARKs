use std::collections::HashMap;

use ndarray::{s, Array, ArrayView, ArrayBase, OwnedRepr, Dimension};

use crate::{precomputes::{compute_bit_orders, get_subdomains, inverse_x, inverse_y}, utils::HALF};

const BASE_CASE_SIZE : u32= 64;
const FOLDS_PER_ROUND : u32 = 3;
const BASE : u32 = 2;
const FOLD_SIZE_RATIO : u32= BASE.pow(FOLDS_PER_ROUND);
const NUM_CHALLENGES : u32 = 80;

pub fn fold<D>(mut values : &Array<u32, D>, coeff: u32, first_round: bool) -> &ArrayBase<OwnedRepr<u32>, D>
where D : Dimension {

    let (rbos, folded_rbos) = compute_bit_orders();
    let sub_domains = get_subdomains();
    let invx = inverse_x(sub_domains);
    let invy = inverse_y(sub_domains);

    for i in 0..FOLDS_PER_ROUND{
        let full_len = values.shape()[1];
        let half_len = full_len/2;
        let new_shape =(half_len,);
        let left = values.slice(s![..;2]);
        let right = values.slice(s![1..;2]);
        let mut f0 = Array::zeros(new_shape);
        for j in 0..half_len{
            let value = (left[j]+right[j])*HALF;
            f0[j] = value;
        }
        let folded_rbos_slice = folded_rbos[full_len..full_len*2].iter().step_by(2).collect();
        let invy_slice = invy[full_len..full_len*2].to_vec();
        let invx_slice = invx[full_len*2..full_len*3].to_vec();
        let twiddle = if i == 0 && first_round {
            ArrayView::from(folded_rbos_slice.iter().map(|&idx| invy[idx]).collect())
        } else {
            ArrayView::from(folded_rbos_slice.iter().map(|&idx| invx[idx]).collect())
        };
        let left_ndim = left.ndim();
        let twiddle_shape: Vec<usize> = std::iter::once(half_len) 
        .chain(std::iter::repeat(1).take(left_ndim - 1)) 
        .collect();
        let twiddle_box = twiddle.to_owned().into_shape_clone(twiddle_shape).unwrap();
        let mut f1;
        for k in 0..half_len{
            let value = (left[k] - right[k])*HALF*twiddle_box[k];
        }
        let mut new_values;
        for l in 0..f0.len(){
            new_values.push(f0[l] + f1[l]*coeff);
        }
        values = new_values;
    }
    values
}

