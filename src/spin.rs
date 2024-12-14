/// Represents the spin at a site on a lattice.
#[derive(Debug, Clone, Copy, PartialEq, Hash, Eq)]
pub enum Spin {
    Up,
    Down,
}

impl Spin {
    /// # Flip
    /// Returns a new spin that is the opposite of the current spin.
    pub fn flip(&self) -> Spin {
        match self {
            Spin::Up => Spin::Down,
            Spin::Down => Spin::Up,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_flip() {
        assert_eq!(Spin::Up.flip(), Spin::Down);
        assert_eq!(Spin::Down.flip(), Spin::Up);
    }
}
