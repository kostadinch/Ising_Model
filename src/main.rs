use rand::prelude::*;

fn main() {
    // Defining initial values.
    let size_of_the_square_matrix = 100;
    let coupling_between_neighboring_spins = 0.44;
    let applied_field = 0.02;

    // Defining initial configuration.
    let mut initial_configuration: Vec<Vec<i32>> = (0..size_of_the_square_matrix)
        .map(|_| {
            (0..size_of_the_square_matrix)
                .map(|_| if rand::random::<bool>() { 1 } else { -1 })
                .collect()
        })
        .collect();

    // Going through sweeps to get a more energetically favorable configuration.
    let number_of_sweeps = 7000;
    let mut final_configuration = initial_configuration.clone();
    for _ in 0..number_of_sweeps {
        initial_configuration = sweep(
            &mut initial_configuration,
            size_of_the_square_matrix,
            coupling_between_neighboring_spins,
            applied_field,
        );
    }

    println!(
        "Final configuration (sample element): {}",
        initial_configuration[0][0]
    );
}

/// # Sweep function
/// Sweep function that implements the Monte Carlo method using the Metropolis-Hastings algorithm.
fn sweep(
    initial_configuration: &mut Vec<Vec<i32>>,
    size: usize,
    coupling_between_neighboring_spins: f64,
    applied_field: f64,
) -> Vec<Vec<i32>> {
    let mut rng = thread_rng();
    let mut updated_configuration = initial_configuration.clone();

    // The loop goes through every element of the initial configuration and checks if changing the
    // element's spin will result in a higher energy configuration.
    for i in 0..size {
        for j in 0..size {
            let spin_up = initial_configuration[i][if j == size - 1 { 0 } else { j + 1 }];
            let spin_down = initial_configuration[i][if j == 0 { size - 1 } else { j - 1 }];
            let spin_left = initial_configuration[if i == 0 { size - 1 } else { i - 1 }][j];
            let spin_right = initial_configuration[if i == size - 1 { 0 } else { i + 1 }][j];

            let initial_spin = initial_configuration[i][j];
            let final_spin = -initial_spin;

            let initial_energy = -coupling_between_neighboring_spins
                * (initial_spin * (spin_up + spin_down + spin_left + spin_right)) as f64
                - applied_field * (initial_spin as f64);

            let final_energy = -coupling_between_neighboring_spins
                * (final_spin * (spin_up + spin_down + spin_left + spin_right)) as f64
                - applied_field * (final_spin as f64);

            let change_in_energy = final_energy - initial_energy;
            let probability_of_change = (-change_in_energy).exp();

            if change_in_energy < 0.0 || rng.gen::<f64>() < probability_of_change {
                updated_configuration[i][j] = final_spin;
            }
        }
    }

    updated_configuration
}
