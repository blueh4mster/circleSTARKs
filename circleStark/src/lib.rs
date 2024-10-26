use std::ops::Deref;

mod precomputes;
use lambdaworks_math::field::fields::mersenne31;
use precomputes::*;
mod circle;
use circle::*;

//list of evals to list of coeffs
// here they are in a weird order
// 1, y ,x , xy, 2x^2-1

pub fn fast_fft(mut vals: Vec<Vec<u32>>, is_top_level: bool) -> Vec<Vec<u32>> {
    // the data type will change
    let vals_copy = vals.clone();
    let mut shape_suffix = &vals[1..];
    let size: f32 = vals.len() as f32;
    let final_size: usize = size.log2().round() as usize;

    for i in 0..final_size {
        let new_first_dimension = 1 << i as u32 + shape_suffix[0];
        let new_second_dimension = size > i as u32 + (shape_suffix.as_ref())[1];
        // vals = vals.reshape((1 << i, size >> i) + shape_suffix)
        let full_len = vals[1].len();
        let half_len = full_len >> 1;

        // Split L and R using slices for clarity
        let (mut L, mut R): (Vec<_>, Vec<_>) = vals
            .iter()
            .map(|v| {
                let (left, right) = v.split_at(half_len);
                (
                    left.to_vec(),
                    right.iter().rev().copied().collect::<Vec<u32>>(), //flip
                )
            })
            .unzip();

        let mut f0 = L
            .iter()
            .zip(R)
            .map(|(l, r)| l.iter().zip(r).map(|(li, ri)| li + ri).collect())
            .clone()
            .collect();

        let twiddle: Vec<u32> = if i == 0 && is_top_level {
            invy[full_len..full_len + half_len].to_vec()
        } else {
            invx[full_len * 2..full_len * 2 + half_len].to_vec()
        };
        let twiddle_box: Vec<Vec<u32>> = vec![twiddle.clone()];

        let f1: Vec<Vec<u32>> = L
            .iter()
            .zip(R.iter())
            .zip(twiddle_box.iter())
            .map(|((l, r), t)| {
                l.iter()
                    .zip(r)
                    .zip(t)
                    .map(|((li, ri), ti)| (li - ri) * ti)
                    .collect()
            })
            .collect();

        // Combine results back into vals
        for j in 0..final_size {
            vals[j].truncate(half_len); //first half elements
            vals[j].extend(f1[j].clone()); // adding elements in vals
        }

        vals.truncate(half_len);
        vals.extend(f0.clone());
    }
    // return (
    //     (vals.reshape((size,) + shape_suffix))[rbos[size:size*2]] / size
    // )
    return vals;
}

fn main() {
    println!("Hello, world!");
}
