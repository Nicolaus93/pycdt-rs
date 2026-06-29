use numpy::ndarray::Array2;
use numpy::{IntoPyArray, PyArray2, PyReadonlyArray2};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

#[pyclass(name = "PyTriangulation")]
pub struct PyTriangulation {
    inner: pycdt_core::Triangulation,
}

#[pymethods]
impl PyTriangulation {
    #[getter]
    fn points<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<f64>>> {
        let n = self.inner.points.len();
        let flat: Vec<f64> = self.inner.points.iter().flat_map(|&p| p).collect();
        Ok(Array2::from_shape_vec((n, 2), flat)
            .map_err(|e| PyValueError::new_err(e.to_string()))?
            .into_pyarray(py))
    }

    #[getter]
    fn triangle_vertices<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<u64>>> {
        let m = self.inner.triangle_vertices.len();
        let flat: Vec<u64> = self
            .inner
            .triangle_vertices
            .iter()
            .flat_map(|&t| t.map(|x| x as u64))
            .collect();
        Ok(Array2::from_shape_vec((m, 3), flat)
            .map_err(|e| PyValueError::new_err(e.to_string()))?
            .into_pyarray(py))
    }

    #[getter]
    fn triangle_neighbors<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<i64>>> {
        let m = self.inner.triangle_neighbors.len();
        let flat: Vec<i64> = self
            .inner
            .triangle_neighbors
            .iter()
            .flat_map(|&n| {
                n.map(|x| {
                    if x == pycdt_core::types::NO_NEIGHBOR {
                        -1i64
                    } else {
                        x as i64
                    }
                })
            })
            .collect();
        Ok(Array2::from_shape_vec((m, 3), flat)
            .map_err(|e| PyValueError::new_err(e.to_string()))?
            .into_pyarray(py))
    }

    #[getter]
    fn constrained_edges(&self) -> Vec<(usize, usize)> {
        self.inner.constrained_edges.iter().copied().collect()
    }

    #[getter]
    fn num_triangles(&self) -> usize {
        self.inner.num_triangles()
    }

    #[getter]
    fn num_points(&self) -> usize {
        self.inner.num_points()
    }

    fn _set_triangles(
        &mut self,
        vertices: PyReadonlyArray2<u64>,
        neighbors: PyReadonlyArray2<i64>,
    ) -> PyResult<()> {
        let verts = vertices.as_array();
        let neighs = neighbors.as_array();
        self.inner.triangle_vertices = verts
            .rows()
            .into_iter()
            .map(|r| [r[0] as usize, r[1] as usize, r[2] as usize])
            .collect();
        self.inner.triangle_neighbors = neighs
            .rows()
            .into_iter()
            .map(|r| {
                r.iter()
                    .map(|&x| {
                        if x < 0 {
                            pycdt_core::types::NO_NEIGHBOR
                        } else {
                            x as usize
                        }
                    })
                    .collect::<Vec<_>>()
                    .try_into()
                    .unwrap()
            })
            .collect();
        Ok(())
    }
}

#[pyfunction]
fn triangulate(points: PyReadonlyArray2<f64>) -> PyResult<PyTriangulation> {
    let arr = points.as_array();
    let pts: Vec<[f64; 2]> = arr.rows().into_iter().map(|r| [r[0], r[1]]).collect();
    let inner = pycdt_core::build::triangulate(&pts);
    Ok(PyTriangulation { inner })
}

#[pyfunction]
fn add_constraints(t: &mut PyTriangulation, edges: Vec<(usize, usize)>) -> PyResult<()> {
    pycdt_core::constrained::add_constraints(&mut t.inner, &edges);
    Ok(())
}

#[pyfunction]
fn remove_holes(t: &mut PyTriangulation, holes: PyReadonlyArray2<f64>) -> PyResult<()> {
    let arr = holes.as_array();
    let pts: Vec<[f64; 2]> = arr.rows().into_iter().map(|r| [r[0], r[1]]).collect();
    pycdt_core::build::remove_holes(&mut t.inner, &pts);
    Ok(())
}

#[pyfunction]
fn remove_holes_by_edges(t: &mut PyTriangulation, edges: Vec<(usize, usize)>) -> PyResult<()> {
    pycdt_core::build::remove_holes_by_edges(&mut t.inner, &edges);
    Ok(())
}

#[pyfunction]
fn remove_super_triangle(t: &mut PyTriangulation) -> PyResult<()> {
    pycdt_core::build::remove_super_triangle(&mut t.inner);
    Ok(())
}

#[pyfunction]
fn update_triangulation(
    t: &mut PyTriangulation,
    new_points: PyReadonlyArray2<f64>,
) -> PyResult<()> {
    let arr = new_points.as_array();
    let pts: Vec<[f64; 2]> = arr.rows().into_iter().map(|r| [r[0], r[1]]).collect();
    pycdt_core::build::update_triangulation(&mut t.inner, &pts);
    Ok(())
}

#[pyfunction]
fn build_polygons_from_edges(_t: &PyTriangulation, edges: Vec<(usize, usize)>) -> Vec<Vec<usize>> {
    pycdt_core::build::build_polygons_from_edges(&edges)
}

#[pymodule]
fn pycdt_rs(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyTriangulation>()?;
    m.add_function(wrap_pyfunction!(triangulate, m)?)?;
    m.add_function(wrap_pyfunction!(add_constraints, m)?)?;
    m.add_function(wrap_pyfunction!(remove_holes, m)?)?;
    m.add_function(wrap_pyfunction!(remove_holes_by_edges, m)?)?;
    m.add_function(wrap_pyfunction!(remove_super_triangle, m)?)?;
    m.add_function(wrap_pyfunction!(update_triangulation, m)?)?;
    m.add_function(wrap_pyfunction!(build_polygons_from_edges, m)?)?;
    Ok(())
}
