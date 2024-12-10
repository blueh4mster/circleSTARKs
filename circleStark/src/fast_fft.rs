use crate::circle::CirclePoint;
use crate::utils::*;
use ndarray::s;
use ndarray::{Array, ArrayView, Axis, Dimension};

// yet to fix new functions with ndarray and reverse bit order

//list of evals to list of coeffs
// here they are in a weird order
// 1, y ,x , xy, 2x^2-1
fn fft<D, T>(
    vals: &Array<T, D>,
    is_top_level: bool,
    rbos: &[usize],
    invx: &[T],
    invy: &[T],
) -> Array<T, D>
where
    D : Dimension
    T: Clone,
{
    let mut vals = vals.to_owned();
    let size = vals.len_of(Axis(0));
    let log_size = log2(size as u32);

    for i in 0..log_size {
        let sub_blocks = 1 << i;
        let block_size = size >> i;

        vals.to_shape((sub_blocks, block_size)).unwrap();
        let full_len = block_size;
        let half_len = full_len >> 1;

        let l = vals.slice(s![.., ..half_len]).to_owned();
        let mut r = vals.slice(s![.., half_len..]).to_owned();
        r.flip_axis(Axis(1));

        let f0 = &l + &r;

        // Select twiddle factors based on level and top-level flag
        let twiddle = if i == 0 && is_top_level {
            ArrayView::from(&invy[full_len..full_len + half_len])
        } else {
            ArrayView::from(&invx[full_len * 2..full_len * 2 + half_len])
        };

        // Reshape twiddle for broadcasting
        let twiddle_box = twiddle.into_shape((1, half_len)).unwrap();

        let f1 = (&l - &r) * &twiddle_box;

        // Update the original array
        vals.slice_mut(s![.., ..half_len]).assign(&f0);
        vals.slice_mut(s![.., half_len..]).assign(&f1);
    }

    // Reshape and apply reverse bit order, then divide by size
    let reshaped = vals.into_shape(size).unwrap();
    let reordered = Array::from_shape_vec(
        size,
        rbos[size..size * 2]
            .iter()
            .map(|&idx| reshaped[idx])
            .collect(),
    )
    .unwrap();

    reordered / T::from(size).unwrap()
}

fn inv_fft<D, T>(vals: &Array<T, D>, sub_domains: &[CirclePoint], rbos: &[usize]) -> Array<T, D>
where
    D: Dimension,
    T: Clone,
{
    let mut vals = reverse_bit_order(vals);
    let size = vals.len_of(Axis(0));
    let log_size = log2(size as u32);

    for i in (0..log_size).rev() {
        let sub_blocks = 1 << i;
        let block_size = size >> i;

        vals.to_shape((sub_blocks, block_size)).unwrap();
        let full_len = block_size;
        let half_len = full_len >> 1;

        let f0 = vals.slice(s![.., ..half_len]).to_owned();
        let f1 = vals.slice(s![.., half_len..]).to_owned();

        // Select twiddle factors based on level
        let twiddle = if i == 0 {
            ArrayView::from(
                &sub_domains[full_len..full_len + half_len]
                    .iter()
                    .map(|p| p.y)
                    .collect::<Vec<_>>(),
            )
        } else {
            ArrayView::from(
                &sub_domains[full_len * 2..full_len * 2 + half_len]
                    .iter()
                    .map(|p| p.x)
                    .collect::<Vec<_>>(),
            )
        };

        let f1_times_twiddle = &f1 * &twiddle.to_shape((1, half_len)).unwrap();

        let l = &f0 + &f1_times_twiddle;
        let r = &f0 - &f1_times_twiddle;

        vals.slice_mut(s![.., ..half_len]).assign(&l);

        let mut r_flipped = r.to_owned();
        r_flipped.flip_axis(Axis(1));
        vals.slice_mut(s![.., half_len..]).assign(&r_flipped);
    }

    vals.to_shape(size).unwrap()
}

fn bary_eval<T>(vals: &Array<T, Ix1>, pt: &CirclePoint, invx: &[T], invy: &[T]) -> T
where
    T: Clone,
{
    let size = vals.len();
    let log_size = log2(size);
    let mut vals = vals.to_owned();

    let mut baryfac = pt.y;
    let mut current_level = 0;

    for i in 0..log_size {
        let full_len = vals.len();
        let half_len = full_len >> 1;

        let l = vals[..half_len].to_owned();
        let mut r = vals[half_len..].to_owned();
        r.flip();

        let f0 = &l + &r;

        let twiddle = if current_level == 0 {
            ArrayView::from(&invy[full_len..full_len + half_len])
        } else {
            ArrayView::from(&invx[full_len * 2..full_len * 2 + half_len])
        };

        if current_level == 1 {
            baryfac = pt.x;
        } else if current_level > 1 {
            baryfac = baryfac * baryfac * T::from(2).unwrap() - T::from(1).unwrap();
        }

        let twiddle_box = twiddle.into_shape((half_len, 1)).unwrap();
        let f1 = (&l - &r) * &twiddle_box;

        vals = f0 + baryfac * &f1;
        current_level += 1;
    }

    vals[0] / T::from(size).unwrap()
}
