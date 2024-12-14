use grid::Grid;

pub mod grid;
pub mod spin;

fn main() {
    // Defining initial values.
    let size_of_the_square_matrix = 100;
    let coupling_between_neighboring_spins = 0.44;
    let applied_field = 0.02;

    let mut grid = Grid::new_random(size_of_the_square_matrix, size_of_the_square_matrix);

    let number_of_sweeps = 7000;
    for step in 0..number_of_sweeps {
        if step % 100 == 0 {
            println!("Sweep number: {}", step);
        }
        grid.step(coupling_between_neighboring_spins, applied_field);
    }

    println!("Final configuration (sample element): {:?}", grid.get(1, 1));
}
