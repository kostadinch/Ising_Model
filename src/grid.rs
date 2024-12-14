use crate::spin::Spin;

/// # Grid
/// This is a struct that represents a grid of spins.
#[derive(Debug)]
pub struct Grid {
    spins: Vec<Spin>,
    width: usize,
    height: usize,
}

impl Grid {
    /// # New random grid
    /// This function creates a new grid of spins, where each spin has a random orientation.
    pub fn new_random(width: usize, height: usize) -> Self {
        let spins = (0..width * height)
            .map(|_| {
                if rand::random::<bool>() {
                    Spin::Up
                } else {
                    Spin::Down
                }
            })
            .collect();

        Self {
            spins,
            width,
            height,
        }
    }

    /// # New constant grid
    /// This function creates a new grid of spins, where each spin has the same orientation.
    pub fn new_constant(width: usize, height: usize, spin: Spin) -> Self {
        let spins = vec![spin; width * height];

        Self {
            spins,
            width,
            height,
        }
    }

    /// # Get index
    /// This function gets the index of a spin at the given coordinates. It applies periodic
    /// boundary conditions to the input coordinates.
    fn get_index(&self, x: i64, y: i64) -> usize {
        // The modulo operator can return negative values, so we add the width/height after the
        // first modulo to ensure that the result is always positive.
        let x_periodic = ((x % self.width as i64) + self.width as i64) % self.width as i64;
        let y_periodic = ((y % self.height as i64) + self.height as i64) % self.height as i64;
        (y_periodic * self.width as i64 + x_periodic) as usize
    }

    /// # Get a spin
    /// This retrieves the spin at the given coordinates, also accounting for periodic boundary
    /// conditions.
    pub fn get(&self, x: i64, y: i64) -> Spin {
        let index = self.get_index(x, y);
        self.spins[index]
    }

    /// # Get a spin as a plus/minus one
    /// This retrieves the spin at the given coordinates as a plus/minus one, also accounting for
    /// periodic boundary conditions.
    pub fn get_spin_as_float(&self, x: i64, y: i64) -> f64 {
        match self.get(x, y) {
            Spin::Up => 1.0,
            Spin::Down => -1.0,
        }
    }

    /// # Set a spin
    /// This sets the spin at the given coordinates, also accounting for periodic boundary
    /// conditions.
    pub fn set(&mut self, x: i64, y: i64, spin: Spin) {
        let index = self.get_index(x, y);
        self.spins[index] = spin;
    }

    /// # Get field energy
    /// Gets the magnetic field energy at a site.
    fn field_energy(&self, x: i64, y: i64, field: f64) -> f64 {
        // Get the nearest neighbours and the spin at the site.
        let our_spin = self.get_spin_as_float(x, y);
        let upper_neighbor = self.get_spin_as_float(x, y + 1);
        let lower_neighbor = self.get_spin_as_float(x, y - 1);
        let left_neighbor = self.get_spin_as_float(x - 1, y);
        let right_neighbor = self.get_spin_as_float(x + 1, y);

        // Calculate the magnetic field energy.
        (our_spin + upper_neighbor + lower_neighbor + left_neighbor + right_neighbor) * field
    }

    /// # Get the interaction energy
    /// Gets the interaction energy at a site.
    fn interaction_energy(&self, x: i64, y: i64, coupling: f64) -> f64 {
        // Get the nearest neighbours and the spin at the site.
        let our_spin = self.get_spin_as_float(x, y);
        let upper_neighbor = self.get_spin_as_float(x, y + 1);
        let lower_neighbor = self.get_spin_as_float(x, y - 1);
        let left_neighbor = self.get_spin_as_float(x - 1, y);
        let right_neighbor = self.get_spin_as_float(x + 1, y);

        // Calculate the interaction energy.
        -coupling * our_spin * (upper_neighbor + lower_neighbor + left_neighbor + right_neighbor)
    }

    /// # Get total energy
    /// Gets the total energy at a site.
    pub fn total_energy(&self, x: i64, y: i64, coupling: f64, field: f64) -> f64 {
        self.interaction_energy(x, y, coupling) + self.field_energy(x, y, field)
    }

    /// # Single site step
    /// This function performs a single Monte Carlo step at a single site.
    pub fn single_site_step(&mut self, x: i64, y: i64, coupling: f64, field: f64) {
        // Get the current energy at the site.
        let current_energy = self.total_energy(x, y, coupling, field);

        // Flip the spin.
        let current_spin = self.get(x, y);
        let new_spin = current_spin.flip();
        self.set(x, y, new_spin);

        // Get the new energy at the site.
        let new_energy = self.total_energy(x, y, coupling, field);

        // Calculate exp(-ΔE); this is the probability of accepting the new configuration.
        let probability_of_acceptance = (-(new_energy - current_energy).exp()).min(1.0);

        // Create a random number between 0 and 1.
        let random_number = rand::random::<f64>();

        // If the random number is less than the probability of accepting the new
        // configuration, accept the new configuration.
        if random_number > probability_of_acceptance {
            self.set(x, y, current_spin);
        }
    }

    /// # Step
    /// This function performs a single Monte Carlo step.
    pub fn step(&mut self, coupling: f64, field: f64) {
        // Iterate over all the spins.
        for y in 0..self.height {
            for x in 0..self.width {
                self.single_site_step(x as i64, y as i64, coupling, field);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use std::collections::HashSet;

    use super::*;

    #[test]
    fn test_new_random() {
        let width = 50;
        let height = 50;
        let grid = Grid::new_random(width, height);
        assert_eq!(grid.spins.len(), width * height);

        // Make sure that not all the spins are the same.
        let unique_spins = grid.spins.iter().collect::<HashSet<_>>();
        assert!(unique_spins.len() > 1);
    }

    #[test]
    fn test_new_constant() {
        let width = 50;
        let height = 50;
        let grid = Grid::new_constant(width, height, Spin::Up);
        assert_eq!(grid.spins.len(), width * height);

        // Make sure that all the spins are the same.
        let unique_spins = grid.spins.iter().collect::<HashSet<_>>();
        assert_eq!(unique_spins.len(), 1);

        // Make sure that all the spins are the same as the one we set.
        let spin_value = **unique_spins.iter().next().unwrap();
        assert_eq!(spin_value, Spin::Up);
    }

    #[test]
    fn test_get() {
        let width = 50;
        let height = 50;
        let mut grid = Grid::new_constant(width, height, Spin::Up);
        grid.set(0, 0, Spin::Down);

        // Test the periodic boundary conditions that get is using.
        assert_eq!(grid.get(0, 0), Spin::Down);
        assert_eq!(grid.get(50, 0), Spin::Down);
        assert_eq!(grid.get(0, 50), Spin::Down);
        assert_eq!(grid.get(500, 500), Spin::Down);
        assert_eq!(grid.get(-50, 0), Spin::Down);
    }

    #[test]
    fn test_get_spin_as_float() {
        let width = 50;
        let height = 50;
        let mut grid = Grid::new_constant(width, height, Spin::Up);
        grid.set(0, 0, Spin::Down);

        // Test the periodic boundary conditions that get_spin_as_float is using.
        assert_eq!(grid.get_spin_as_float(0, 0), -1.0);
        assert_eq!(grid.get_spin_as_float(50, 0), -1.0);
        assert_eq!(grid.get_spin_as_float(0, 50), -1.0);
        assert_eq!(grid.get_spin_as_float(500, 500), -1.0);
        assert_eq!(grid.get_spin_as_float(-50, 0), -1.0);
    }

    #[test]
    fn test_get_index() {
        let width = 50;
        let height = 50;
        let grid = Grid::new_constant(width, height, Spin::Up);

        // Test the periodic boundary conditions that get_index is using.
        assert_eq!(grid.get_index(0, 0), 0);
        assert_eq!(grid.get_index(50, 0), 0);
        assert_eq!(grid.get_index(3, 50), 3);
        assert_eq!(grid.get_index(500, 500), 0);
        assert_eq!(grid.get_index(-50, 0), 0);
        assert_eq!(grid.get_index(-1, 0), 49);
    }

    #[test]
    fn test_set() {
        let width = 50;
        let height = 50;
        let mut grid = Grid::new_constant(width, height, Spin::Up);
        grid.set(65, 14, Spin::Down);
        grid.set(-1, 14, Spin::Down);

        assert_eq!(grid.get(15, 14), Spin::Down);
        assert_eq!(grid.get(49, 14), Spin::Down);
    }

    #[test]
    fn test_field_energy() {
        let width = 50;
        let height = 50;
        let grid = Grid::new_constant(width, height, Spin::Up);
        assert_eq!(grid.field_energy(0, 0, 1.0), 4.0);
    }

    #[test]
    fn test_interaction_energy() {
        let width = 50;
        let height = 50;
        let grid = Grid::new_constant(width, height, Spin::Up);
        assert_eq!(grid.interaction_energy(0, 0, 1.0), -4.0);
    }
}
