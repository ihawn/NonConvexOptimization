extern crate ndarray;
extern crate ndarray_linalg;

use ndarray::prelude::*;
use ndarray_linalg::Solve;
use crate::structs::{Hyperplane, Pyramid};

//Function to generate all valid intersections given 3 pyramids. Should only be used in initial iteration
pub fn intersect_pyramids(p1: Pyramid, p2: Pyramid, p3: Pyramid) -> Vec<[f64; 3]>
{
    let default_hyp = Hyperplane { parent_id: 0, direction: 0, coeff: [0.0, 0.0, 0.0, 0.0] };
    let mut hyp_arr: [Hyperplane; 12] = [default_hyp; 12];
    let mut int_points: Vec<[f64; 3]> = vec![];


    //Build list of hyperplanes to intersect
    for i in 0..4
    {
        hyp_arr[3*i] = p1.hyperplanes[i];
        hyp_arr[3*i+1] = p2.hyperplanes[i];
        hyp_arr[3*i+2] = p3.hyperplanes[i];
    }

    for i in 0..hyp_arr.len()
    {
        for j in i+1..hyp_arr.len()
        {
            for k in j+1..hyp_arr.len()
            {
                if hyp_arr[i].direction != hyp_arr[j].direction &&
                   hyp_arr[j].direction != hyp_arr[k].direction &&
                   hyp_arr[i].direction != hyp_arr[k].direction
                {
                    //let pt = intersect_hyperplanes(hyp_arr[i], hyp_arr[j], hyp_arr[k]);
                    
                    //if valid_intersection(p1, p2, p3, pt) { int_points.push(pt); }
                }
            }
        }
    }

    int_points
}

pub fn intersect_hyperplanes(h1: Hyperplane, h2: Hyperplane, h3: Hyperplane) -> Array1<f64>
{
    let a: Array2<f64> = array![
        [h1.coeff[0], h1.coeff[1], h1.coeff[2]],
        [h2.coeff[0], h2.coeff[1], h2.coeff[2]],
        [h3.coeff[0], h3.coeff[1], h3.coeff[2]]];


    let b: Array1<f64> = array![
        -h1.coeff[3],
        -h2.coeff[3],
        -h3.coeff[3]];

    a.solve_into(b).unwrap()
}