use structs::Hyperplane;

mod structs;
mod intersections;
mod cone_search;

fn main()
{
    cone_search::solve();
}
