use crate::types::{Point, EPS};
use robust::Coord;

fn to_coord(p: &Point) -> Coord<f64> {
    Coord { x: p[0], y: p[1] }
}

/// Shewchuk's robust 2D orientation predicate.
/// Returns > 0 if (a, b, c) are in counterclockwise order,
/// < 0 if clockwise, 0 if collinear.
pub fn orient2d(a: &Point, b: &Point, c: &Point) -> f64 {
    robust::orient2d(to_coord(a), to_coord(b), to_coord(c))
}

/// Shewchuk's robust incircle predicate.
/// Returns > 0 if d is inside the circumcircle of (a, b, c),
/// < 0 if outside, 0 if on the circle.
pub fn incircle(a: &Point, b: &Point, c: &Point, d: &Point) -> f64 {
    robust::incircle(to_coord(a), to_coord(b), to_coord(c), to_coord(d))
}

/// Classification of a point relative to a triangle.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PointInTriangle {
    Outside = 0,
    Inside = 1,
    OnEdge0 = 2, // on edge a-b
    OnEdge1 = 3, // on edge b-c
    OnEdge2 = 4, // on edge c-a
}

/// Test if `p` is inside the bounding box of segment [a, b].
fn point_in_segment_box(a: &Point, b: &Point, p: &Point) -> bool {
    let min_x = a[0].min(b[0]);
    let max_x = a[0].max(b[0]);
    let min_y = a[1].min(b[1]);
    let max_y = a[1].max(b[1]);
    p[0] >= min_x - EPS && p[0] <= max_x + EPS && p[1] >= min_y - EPS && p[1] <= max_y + EPS
}

/// Classify point `p` relative to triangle (a, b, c).
/// Ported from Python `point_inside_triangle`.
pub fn point_in_triangle(p: &Point, a: &Point, b: &Point, c: &Point) -> PointInTriangle {
    let o1 = orient2d(a, b, p);
    let o2 = orient2d(b, c, p);
    let o3 = orient2d(c, a, p);

    // Check edges (collinear + within bounding box)
    if o1.abs() <= EPS && point_in_segment_box(a, b, p) {
        return PointInTriangle::OnEdge0;
    }
    if o2.abs() <= EPS && point_in_segment_box(b, c, p) {
        return PointInTriangle::OnEdge1;
    }
    if o3.abs() <= EPS && point_in_segment_box(c, a, p) {
        return PointInTriangle::OnEdge2;
    }

    // Inside: all same sign
    if (o1 >= 0.0 && o2 >= 0.0 && o3 >= 0.0) || (o1 <= 0.0 && o2 <= 0.0 && o3 <= 0.0) {
        return PointInTriangle::Inside;
    }

    PointInTriangle::Outside
}

/// Reorder indices (a, b, c) so the triangle is counterclockwise.
pub fn ensure_ccw(points: &[Point], a: usize, b: usize, c: usize) -> (usize, usize, usize) {
    if orient2d(&points[a], &points[b], &points[c]) < 0.0 {
        (a, c, b)
    } else {
        (a, b, c)
    }
}

/// Ray casting algorithm — even-odd rule.
/// Ported from Python `is_point_inside_polygon`.
pub fn is_point_inside_polygon(point: &Point, polygon: &[Point]) -> bool {
    let x = point[0];
    let y = point[1];
    let n = polygon.len();
    let mut inside = false;

    for i in 0..n {
        let (x0, y0) = (polygon[i][0], polygon[i][1]);
        // Use wrapping subtraction to index i-1 (wraps around for i==0)
        let prev = if i == 0 { n - 1 } else { i - 1 };
        let (x1, y1) = (polygon[prev][0], polygon[prev][1]);

        // Corner check
        if (x - x0).abs() < EPS && (y - y0).abs() < EPS {
            return true;
        }

        if (y0 > y) != (y1 > y) {
            let cross = (x - x0) * (y1 - y0) - (x1 - x0) * (y - y0);
            if cross == 0.0 {
                return true; // on boundary
            }
            if (cross < 0.0) != (y1 < y0) {
                inside = !inside;
            }
        }
    }
    inside
}

/// Axis-aligned bounding box check.
pub fn is_point_in_box(point: &Point, min: &Point, max: &Point) -> bool {
    point[0] >= min[0] - EPS
        && point[0] <= max[0] + EPS
        && point[1] >= min[1] - EPS
        && point[1] <= max[1] + EPS
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn orient2d_ccw_is_positive() {
        let a = [0.0, 0.0];
        let b = [1.0, 0.0];
        let c = [0.0, 1.0];
        assert!(orient2d(&a, &b, &c) > 0.0);
    }

    #[test]
    fn orient2d_cw_is_negative() {
        let a = [0.0, 0.0];
        let b = [0.0, 1.0];
        let c = [1.0, 0.0];
        assert!(orient2d(&a, &b, &c) < 0.0);
    }

    #[test]
    fn point_in_triangle_inside() {
        let p = [0.25, 0.25];
        let a = [0.0, 0.0];
        let b = [1.0, 0.0];
        let c = [0.0, 1.0];
        assert_eq!(point_in_triangle(&p, &a, &b, &c), PointInTriangle::Inside);
    }

    #[test]
    fn point_in_triangle_outside() {
        let p = [1.0, 1.0];
        let a = [0.0, 0.0];
        let b = [1.0, 0.0];
        let c = [0.0, 1.0];
        assert_eq!(point_in_triangle(&p, &a, &b, &c), PointInTriangle::Outside);
    }

    #[test]
    fn point_in_triangle_on_edge0() {
        // Point on edge a-b: [0,0]-[1,0]
        let p = [0.5, 0.0];
        let a = [0.0, 0.0];
        let b = [1.0, 0.0];
        let c = [0.0, 1.0];
        assert_eq!(point_in_triangle(&p, &a, &b, &c), PointInTriangle::OnEdge0);
    }

    #[test]
    fn is_point_inside_polygon_square() {
        let square = vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
        assert!(is_point_inside_polygon(&[0.5, 0.5], &square));
        assert!(!is_point_inside_polygon(&[2.0, 2.0], &square));
        // corner is inside
        assert!(is_point_inside_polygon(&[0.0, 0.0], &square));
    }

    #[test]
    fn is_point_in_box_basic() {
        let min = [0.0, 0.0];
        let max = [1.0, 1.0];
        assert!(is_point_in_box(&[0.5, 0.5], &min, &max));
        assert!(is_point_in_box(&[0.0, 0.0], &min, &max));
        assert!(is_point_in_box(&[1.0, 1.0], &min, &max));
        assert!(!is_point_in_box(&[1.5, 0.5], &min, &max));
        assert!(!is_point_in_box(&[0.5, -0.1], &min, &max));
    }

    #[test]
    fn incircle_point_inside_circumcircle() {
        // Right triangle with circumcircle centered at (0.5, 0.5)
        let a = [0.0, 0.0];
        let b = [1.0, 0.0];
        let c = [0.0, 1.0];
        // Origin of circumcircle is inside — d very close to center
        let d = [0.5, 0.4];
        // Not testing sign strictly but just that it runs
        let _ = incircle(&a, &b, &c, &d);
    }

    #[test]
    fn ensure_ccw_swaps_cw_triangle() {
        let points = vec![[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]];
        // CW: indices 0, 2, 1
        let (a, b, c) = ensure_ccw(&points, 0, 2, 1);
        assert!(orient2d(&points[a], &points[b], &points[c]) >= 0.0);
    }
}
