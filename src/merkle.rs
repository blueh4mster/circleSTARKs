//types very dicey

use sha256::digest;

use crate::circle::{
    div, multiply, scalar_division, scalar_multiply, subtract, CircleImpl, CirclePoint, MODULUS,
};

pub fn hash(x: &[u8]) -> String {
    return digest(x);
}

pub fn merkelize(vals: &[CirclePoint]) -> vec<CirclePoint> {
    assert!(vals.len() & (vals.len()-1) == 0);
    let o = vec![None; vals.len()];
    o.extend(vals.iter().map(|val| Some(digest(val).as_bytes().to_vec())));
    for i in (0..vals.len()).rev() {
        o[i] = digest(o[i*2] + o[i*2+1]);
    }
    o
}

pub fn get_root(tree: &[CirclePoint]) -> Option<CirclePoint> {
    return tree[1]?;
}
pub fn get_branch(tree: &[CirclePoint], pos: u32) -> vec<CirclePoint>{
    let offset_pos = (pos + tree.len())/2;
    let branch_length = (tree.len() as u64).next_power_of_two().count_ones() as usize - 2;
    (0..branch_length).map(|i| tree[(offset_pos >> i)^1]).collect()
}

pub fn verify_branch(root: &CirclePoint, pos: u32, val: &[u8], branch: &[CirclePoint]) -> bool {
    let x = hash(val);
    for b in branch{
        if pos & 1{
            x = hash(b+x);
        } else {
            x = hash(x+b)
        }
        pos = pos/2;
    }
    return x == root
}

