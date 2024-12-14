use ndarray::{s, Array, Axis, Dimension};
use std::any::Any;
use crate::circle::{CirclePoint}

fn is_tuple(value: &dyn Any) -> bool {
    value.is::<(CirclePoint, CirclePoint)>()
}
pub fn log2(x: u32) -> usize {
    (x as f64).log2() as usize
}
pub fn reverse_bit_order<D>(vals: &Array<u32, D>) -> Array<u32, D>
where
    D: Dimension,
{
    let size = vals.len_of(Axis(0));
    let log_size = (size as u32).ilog2() as usize;
    let mut current = vals.to_owned();

    for i in 0..log_size {
        let sub_blocks = 1 << i;
        let block_size = size >> i;

        // Reshape the array
        let reshaped = current
            .to_shape((sub_blocks, block_size))
            .expect("Reshape failed");

        // Create output array
        let mut output = Array::zeros(reshaped.raw_dim());

        // Split each block into left and right halves
        for (block_idx, block) in reshaped.rows().into_iter().enumerate() {
            let half_len = block.len() / 2;

            // Left half goes in first half of output
            output
                .slice_mut(s![block_idx, ..half_len])
                .assign(&block.slice(s![..;2]));

            // Right half goes in second half of output
            output
                .slice_mut(s![block_idx, half_len..])
                .assign(&block.slice(s![1..;2]));
        }

        current = output; //yet to check
    }

    // Reshape back to original shape
    current
}

pub fn folded_reverse_bit_order<D>(vals: &Array<u32, D>) -> Array<u32, D>
where
    D: Dimension,
{
    // Create an output array of the same shape as input
    let mut output = Array::zeros(vals.raw_dim());

    // Handle even indices (0, 2, 4, ...)
    let even_indices = vals.slice(s![..;2]);
    let reversed_even = reverse_bit_order(&even_indices);
    output.slice_mut(s![..;2]).assign(&reversed_even);

    // Handle odd indices (1, 3, 5, ...), but reverse the order before and after bit reversal
    let odd_indices = vals.slice(s![1..;2]);
    let reversed_odd = reverse_bit_order(&odd_indices.slice(s![..;-1])).slice(s![..;-1]);
    output.slice_mut(s![1..;2]).assign(&reversed_odd);

    output
}
