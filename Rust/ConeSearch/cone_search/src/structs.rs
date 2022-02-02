//Structs to use throughout the program

#[derive(Copy, Clone)]
pub struct Hyperplane
{
    pub parent_id: usize,
    pub direction: u8, //0-3
    pub coeff: [f64; 4]
}

#[derive(Copy, Clone)]
pub struct Pyramid
{
    pub id: usize,
    pub peak: [f64; 3],
    pub ell: f64,
    pub dist: f64, //distance to another pyramid. Sort of a temp variable but it's more convenient to store it here
    pub hyperplanes: [Hyperplane; 4]
}