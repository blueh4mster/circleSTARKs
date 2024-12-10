use ndarray::{Array, ArrayView, Axis, Dimension, ShapeBuilder};

pub fn reverse_bit_order<D>(vals: &Array<u32, D>) -> Array<u32, D>
{
    let size = vals.len_of(Axis(0));
    let log_size = (size as u32).log2() as usize;
    let mut current = vals.to_owned();

    for i in 0..log_size {
        let sub_blocks = 1 << i;
        let block_size = size >> i;

        // Reshape the array
        let reshaped = current.into_shape((sub_blocks, block_size))
            .expect("Reshape failed");

        // Create output array
        let mut output = Array::zeros(reshaped.raw_dim());

        // Split each block into left and right halves
        for (block_idx, block) in reshaped.rows().into_iter().enumerate() {
            let half_len = block.len() / 2;
            
            // Left half goes in first half of output
            output.slice_mut(s![block_idx, ..half_len]).assign(&block.slice(s![..;2]));
            
            // Right half goes in second half of output
            output.slice_mut(s![block_idx, half_len..]).assign(&block.slice(s![1..;2]));
        }

        current = output;
    }

    // Reshape back to original shape
    current
}