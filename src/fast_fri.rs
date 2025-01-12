use ndarray::{s, Array, ArrayView, ArrayBase, OwnedRepr, Dimension};

use crate::{precomputes::{compute_bit_orders, get_subdomains, inverse_x, inverse_y}, utils::{log2, rbo_index_to_original, HALF}};

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

pub fn fold_with_positions<D>(mut values : &Array<u32, D>, mut domain_size: usize,mut positions: &Array<usize,D>, coeff: u32, first_round: bool) -> &ArrayBase<OwnedRepr<u32>, D>
where D : Dimension {
    let sub_domains = get_subdomains();
    let invx = inverse_x(sub_domains);
    let invy = inverse_y(sub_domains);
    positions = &positions.slice(s![..;2]).to_owned();
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
        let twiddle;
        if i == 0 && first_round{
            let unrbo_positions = rbo_index_to_original(domain_size, positions, true);
            let indices = &unrbo_positions+domain_size;
            twiddle = ArrayView::from(indices.iter().map(|&idx| invy.get(idx)).collect());
        } else {
            let shifted_positions = (positions << 1) >> i;
            let unrbo_positions = rbo_index_to_original(domain_size*2, &shifted_positions, true);
            let indices = &unrbo_positions + (domain_size*2);
            twiddle = ArrayView::from(indices.iter().map(|&idx| invx.get(idx)).collect());
        }
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
        positions = &positions.slice(s![..;2]).to_owned();
        domain_size = domain_size/2;
    }
    values
}

pub fn prove_low_degree<D>(evaluations: &Array<usize, D>, extra_entropy: &Vec<u8>) where D:Dimension{
    //commit merkle root 
    let (rbos, folded_rbos) = compute_bit_orders();
    let indices = &folded_rbos[evaluations.len()..evaluations.len()*2];
    let values = ArrayView::from(indices.iter().map(|&x| evaluations.get(x)).collect());
    let leaves;
    let roots;
    let trees;
    //prove descent
    let rounds = (log2(evaluations.len() as u32) / BASE_CASE_SIZE as usize) / FOLDS_PER_ROUND as usize;
    for i in 0..rounds{
        
    }

}

//     for i in range(rounds):
//         leaves.append(values)
//         trees.append(merkelize_top_dimension(values.reshape(
//             (len(values) // FOLD_SIZE_RATIO, FOLD_SIZE_RATIO)
//             + values.shape[1:]
//         )))
//         roots.append(trees[-1][1])
//         
//         fold_factor = ExtendedM31(get_challenges(b''.join(roots), modulus, 4))
//         values = fold(values, fold_factor, i==0)
//     entropy = extra_entropy + b''.join(roots) + values.tobytes()
//     challenges = get_challenges(
//         entropy, len(evaluations) >> FOLDS_PER_ROUND, NUM_CHALLENGES
//     )
//     round_challenges = (
//         challenges.reshape((1,)+challenges.shape)
//         >> cp.arange(0, rounds * FOLDS_PER_ROUND, FOLDS_PER_ROUND)
//         .reshape((rounds,) + (1,) * challenges.ndim)
//     )

//     branches = [
//         [get_branch(tree, c) for c in r_challenges]
//         for i, (r_challenges, tree) in enumerate(zip(round_challenges, trees))
//     ]
//     round_challenges_xfold = (
//         round_challenges.reshape(round_challenges.shape + (1,)) * 8
//         + cp.arange(FOLD_SIZE_RATIO).reshape(1, 1, FOLD_SIZE_RATIO)
//     )

//     leaf_values = [
//         leaves[i][round_challenges_xfold[i]]
//         for i in range(rounds)
//     ]
//     return {
//         "roots": roots,
//         "branches": branches,
//         "leaf_values": leaf_values,
//         "final_values": values
//     }