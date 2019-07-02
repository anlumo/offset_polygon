use std::fmt;

pub struct CombinatorialExplosionError;

impl fmt::Display for CombinatorialExplosionError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Combinatorial explosion.")
    }
}

impl fmt::Debug for CombinatorialExplosionError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        <Self as fmt::Display>::fmt(self, f)
    }
}

impl std::error::Error for CombinatorialExplosionError {
}
