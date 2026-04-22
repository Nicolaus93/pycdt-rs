pub const NO_NEIGHBOR: usize = usize::MAX;
pub const EPS: f64 = 1e-12;

pub type Point = [f64; 2];
pub type TriangleVertices = [usize; 3];
pub type TriangleNeighbors = [usize; 3];

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PointLocation {
    Interior(usize),
    OnEdge(usize, usize),
    NotFound,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn no_neighbor_is_usize_max() {
        assert_eq!(NO_NEIGHBOR, usize::MAX);
    }

    #[test]
    fn eps_matches_expected_value() {
        assert_eq!(EPS, 1e-12);
    }
}
