//Structs to use throughout the program

pub struct Hyperplane
{
    pub parent_id: usize,
    pub direction: u8, //0-3
    pub coeff: Vec<f64>
}

pub struct Pyramid
{
    pub id: usize,
    pub peak: Vec<f64>,
    pub ell: f64,
    pub dist: f64, //distance to another pyramid. Sort of a temp variable but it's more convenient to store it here
    pub hyperplanes: Vec<Hyperplane>
}