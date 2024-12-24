use sha256::digest;

pub fn hash(x: Vec<u8>) -> Vec<u8> {
    return digest(x).as_bytes().to_vec();
}

pub fn merkelize(vals: Vec<Vec<u8>>) -> Vec<Option<Vec<u8>>> {
    assert!(vals.len() & (vals.len()-1) == 0);
    let mut o = vec![None; vals.len()];
    o.extend(vals.iter().map(|val| Some(digest(val).as_bytes().to_vec())));
    for i in (0..vals.len()).rev() {
        let o1 = o[i*2].clone().unwrap();
        let s1 = o1.as_slice();
        let o2 = o[i*2+1].clone().unwrap();
        let s2 = o2.as_slice();
        let result: Vec<u8> = s1.iter()
        .zip(s2.iter())
        .map(|(a, b)| a + b) // Add corresponding elements
        .collect();
        o[i] = Some(digest(result.as_slice()).as_bytes().to_vec());
    }
    o
}

pub fn get_root(tree: Vec<&[u8]>) -> Option<&[u8]> {
    return Some(tree[1]);
}
pub fn get_branch(tree: Vec<Vec<u8>>, pos: usize) -> Vec<Vec<u8>>{
    let offset_pos = (pos + tree.len())/2;
    let branch_length = (tree.len() as u64).next_power_of_two().count_ones() as usize - 2;
    (0..branch_length).map(|i| tree[(offset_pos >> i)^1].clone()).collect()
}

pub fn verify_branch(root: Vec<u8>, mut pos: i32, val: Vec<u8>, branch: Vec<Vec<u8>>) -> bool {
    let mut x = hash(val);
    for b in branch{
        if pos != 0{
            let result: Vec<u8> = b.iter().zip(x.iter()).map(|(first, second)| first + second).collect();
            x = hash(result);
        } else {
            let result: Vec<u8> = x.iter().zip(b.iter()).map(|(f, s)| f + s).collect();
            x = hash(result);
        }
        pos = pos/2;
    }
    return x == root
}

