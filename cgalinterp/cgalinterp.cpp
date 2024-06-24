/*
  cgalinterp.cpp
 
  Copyright (C) 2023 Ben Vitale
 
  Based on cgalinterp.cpp from PyLidar
  Copyright (C) 2015 John Armston, Pete Bunting, Neil Flood, Sam Gillingham
  https://github.com/ubarsc/cgalinterp/blob/master/src/cgalinterp.cpp 

  Description

  This version differs from the PyLidar original by having a much narrower API,
  much less generally useful, but simpler and more CPU-efficient for the
  purposes of the UCS-flood analysis pipeline.

  In particular, the inputs to our cgalinterp.naturalneighbor function,
  exported to Python and implemented here in C++, are three 1-d numpy arrays,
  'xin', 'yin', 'zin', and a scalar float 'step'. Together, the arrays define a
  single sparse input array of coordinates (x,y) of known 2d points, and the Z
  value associated with each point. The output array is derived implicitly as a
  dense grid of points, sweeping from xin_min to xin_max, and yin_min to
  yin_max, spaced regularly at distance 'step'. The 'naturalneighbor' function
  computes an interpolated Z value for each point in the output grid. The
  return value is a 1-d numpy array of Z values associated with the output
  points, in row-major order, and suitable for writing directly into a
  rasterio dataset.

  License

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
 
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
 
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

// Following needed or CGAL templates misbehave.
#ifndef NDEBUG
#  define NDEBUG
#endif

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION 


#include <math.h>
#include <thread>
#include <vector>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Interpolation_traits_2.h>
#include <CGAL/natural_neighbor_coordinates_2.h>
#include <CGAL/interpolation_functions.h>
#include <CGAL/algorithm.h>
#include <CGAL/Origin.h>
#include <CGAL/squared_distance_2.h>

#include <Python.h>
#include "numpy/arrayobject.h"

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::FT CGALCoordType;
typedef K::Vector_2 CGALVector;
typedef K::Point_2 CGALPoint;
    
typedef CGAL::Delaunay_triangulation_2<K> DelaunayTriangulation;
typedef CGAL::Interpolation_traits_2<K> InterpTraits;
typedef DelaunayTriangulation::Vertex_handle Vertex_handle;
typedef DelaunayTriangulation::Face_handle Face_handle;
    
typedef std::vector<std::pair<CGALPoint, CGALCoordType>> CoordinateVector;
typedef std::map<CGALPoint, CGALCoordType, K::Less_xy_2> PointValueMap;

/* An exception object for this module. */
struct CGALInterpState {
  PyObject *error;
};

static inline CGALInterpState* GetState(PyObject *obj) {
  return (CGALInterpState*) PyModule_GetState(obj);
}

#ifndef __cpp_lib_ssize

template <class T>
constexpr auto ssize(const T& c)
  -> std::common_type_t<std::ptrdiff_t,
			std::make_signed_t<decltype(c.size())>>
{
  using R = std::common_type_t<std::ptrdiff_t,
			       std::make_signed_t<decltype(c.size())>>;
  return static_cast<R>(c.size());
}

#endif  // __cpp_lib_ssize

npy_intp y_start_, y_end_;
int thread_num_;

class Interpolator {
public:
  static std::unique_ptr<Interpolator> Create(
    PyArrayObject *xin, PyArrayObject *yin, PyArrayObject *zin, float step,
    float no_data);

  PyArrayObject* Interpolate();
  
private:
  Interpolator() {}
  DelaunayTriangulation dt_;
  PointValueMap values_;
  std::vector<float> xout_, yout_;
  float no_data_;
  PyArrayObject *out_;

  void InterpolateBlock(npy_intp y_start, npy_intp y_end);
};

// Create output grid points from input points.
static std::vector<float> make_coords(const PyArrayObject *arr, float step) {
  assert(PyArray_TYPE(arr) == NPY_FLOAT);
  assert(PyArray_STRIDE(arr, 0) == sizeof(float));
  assert(PyArray_ITEMSIZE(arr, 0) == sizeof(float));
  
  const float* const arr_data =
    static_cast<float*>(PyArray_DATA(const_cast<PyArrayObject*>(arr)));
  const npy_intp arr_len = PyArray_DIM(arr, 0);
  const auto minmax = std::minmax_element(arr_data, arr_data + arr_len);
  const float min = *minmax.first;
  const float max = *minmax.second;

  const int out_len = 1 + static_cast<int>(std::ceil((max - min) / step));
  std::vector<float> out;
  out.reserve(out_len);
  for (float x = min; x <= max; x += step) {
    out.push_back(x);
  }
  return out;
}

std::unique_ptr<Interpolator> Interpolator::Create(
    PyArrayObject *xin, PyArrayObject *yin, PyArrayObject *zin, float step,
    float no_data) {
  auto interp = std::unique_ptr<Interpolator>(new Interpolator);
  interp->no_data_ = no_data;

  // Setup Delaunay Triangulation and Natural Neighbor context;
  const npy_intp in_len = PyArray_DIM(xin, 0);
  for (npy_intp i = 0; i < in_len; ++i) {
    const float x = *(static_cast<float*>(PyArray_GETPTR1(xin, i)));
    const float y = *(static_cast<float*>(PyArray_GETPTR1(yin, i)));
    const CGALCoordType z = *(static_cast<float*>(PyArray_GETPTR1(zin, i)));
    const K::Point_2 point(x, y);
    interp->dt_.insert(point);
    interp->values_.insert(std::make_pair(point, z));
  }

  // Setup output grid points.
  interp->xout_ = make_coords(xin, step);
  interp->yout_ = make_coords(yin, step);
  
  // Allocate output array.
  npy_intp dims[2] = {
    (npy_intp) interp->yout_.size(),
    (npy_intp) interp->xout_.size()
  };
  interp->out_ = (PyArrayObject*) PyArray_EMPTY(2, dims, NPY_FLOAT, 0);
  if (!interp->out_) return nullptr;
  
  return interp;
}

PyArrayObject *Interpolator::Interpolate() {
  const int num_blocks =
    yout_.size() > 200 ? std::thread::hardware_concurrency() : 1;
                   
  const int rows_per_block = yout_.size() / num_blocks;

  std::vector<std::thread> threads;
  for (npy_intp yi = 0; yi < ssize(yout_); yi += rows_per_block) {
    const int yi_end = std::min(ssize(yout_), yi + rows_per_block);
    threads.push_back(std::thread(
	  [this, yi, yi_end] { InterpolateBlock(yi, yi_end); }));
  }
  for (auto& thread : threads) {
    thread.join();
  }
  return out_;
}

void Interpolator::InterpolateBlock(npy_intp y_start, npy_intp y_end) {
  using NNResultType = CGAL::Triple<
    std::back_insert_iterator<CoordinateVector>, CGALCoordType, bool>;

  DelaunayTriangulation::Face_handle hint;
  for (npy_intp yi = y_start; yi < y_end; ++yi) {
    // Reverse y array to place largest UTM northings in first row of GeoTIFF.
    const float y = yout_[yout_.size() - 1 - yi];
    
    for (npy_intp xi = 0; xi < ssize(xout_); ++xi) {
      K::Point_2 point(xout_[xi], y);
      CoordinateVector coords;
      hint = dt_.locate(point, hint);
      const NNResultType result = CGAL::natural_neighbor_coordinates_2(
        dt_, point, std::back_inserter(coords), hint);
      float *const outp = static_cast<float*>(PyArray_GETPTR2(out_, yi, xi));
      if (!result.third) {
	// point not within convex hull of dataset.
	*outp = no_data_;
      } else {
	const CGALCoordType norm = result.second;
	const CGALCoordType outValue = CGAL::linear_interpolation(
	    coords.begin(), coords.end(), norm,
            CGAL::Data_access<PointValueMap>(values_));
	*outp = outValue;
      }
    }
  }
}

static PyObject *cgalinterp_naturalneighbour(PyObject *self, PyObject *args)
{
  PyArrayObject *xin, *yin, *zin;
  float step, noData = std::nan("");
  
  if (!PyArg_ParseTuple(args, "OOOf|f:NaturalNeighbour",
			&xin, &yin, &zin, &step, &noData))
    return nullptr;
  
  if (!PyArray_Check(xin) || !PyArray_Check(yin) || !PyArray_Check(zin)) {
      PyErr_SetString(GetState(self)->error,
		      "First 3 arguments must be numpy arrays");
      return nullptr;
    }

  // Check inputs are 1-d.
  if (PyArray_NDIM(xin) != 1 || PyArray_NDIM(yin) != 1
      || PyArray_NDIM(zin) != 1) {
      PyErr_SetString(GetState(self)->error,
		      "Input X, Y, and Z arrays must be 1-dimensional");
      return nullptr;    
  }
  
  // Check dimensions equal.
  if ((PyArray_DIM(xin, 0) != PyArray_DIM(yin, 0))
      || (PyArray_DIM(xin, 0) != PyArray_DIM(zin, 0))) {
      PyErr_SetString(GetState(self)->error,
		      "Input X, Y and Z arrays must be the same length");
      return nullptr;
    }
  
  // Check types.
  if ((PyArray_TYPE(xin) != NPY_FLOAT) || (PyArray_TYPE(yin) != NPY_FLOAT)
      || (PyArray_TYPE(zin) != NPY_FLOAT)) {
      PyErr_SetString(GetState(self)->error,
		      "X, Y, and Z input arrays must have type float32");
      return nullptr;
    }
  
  if (PyArray_DIM(xin, 0) < 3) {
      PyErr_SetString(GetState(self)->error,
		      "Not enough points in X, Y, Z arrays; need at least 3.");
      return nullptr;
    }

  try {
    auto interp = Interpolator::Create(xin, yin, zin, step, noData);

    if (!interp) {
      PyErr_SetString(GetState(self)->error, "Failed to create interpolator");
      return nullptr;
    }

    PyArrayObject *out = interp->Interpolate();
    return (PyObject*) out;
  }
  catch (std::exception &e) {
    PyErr_SetString(GetState(self)->error, e.what());
    return nullptr;
  }
}

/* Our list of functions in this module*/
static PyMethodDef CGALInterpMethods[] = {
    {"NaturalNeighbour", cgalinterp_naturalneighbour, METH_VARARGS,
"Perform Natural Neighbour Interpolation\n"
"  arr = NaturalNeighbour(xvals, yvals, zvals, xgrid, ygrid, no_data=NaN)\n"
"\n"
"  xvals, yvals, zvals:\n"
"    1d arrays of x, y, and z values of input points. Must have same length.\n"
"  no_data: \n"
"    Optional float. Points outside convex hull of input will have output z\n"
"    value set to this no_data value.\n"
    }, 
    {nullptr}
};

static int cgalinterp_traverse(PyObject *m, visitproc visit, void *arg) {
  Py_VISIT(GetState(m)->error);
  return 0;
}

static int cgalinterp_clear(PyObject *m)
{
  Py_CLEAR(GetState(m)->error);
  return 0;
}

static struct PyModuleDef moduledef = {
  PyModuleDef_HEAD_INIT,
  "cgalinterp",
  nullptr,
  sizeof(struct CGALInterpState),
  CGALInterpMethods,
  nullptr,
  cgalinterp_traverse,
  cgalinterp_clear,
  nullptr
};

PyMODINIT_FUNC 
PyInit_cgalinterp(void)
{
  PyObject *pModule;
  struct CGALInterpState *state;

  // Initialize numpy.
  import_array();
  
  pModule = PyModule_Create(&moduledef);
  if (!pModule) return nullptr;
  
  state = GetState(pModule);
  
  /* Create and add our exception type */
  state->error = PyErr_NewException("cgal.error", nullptr, nullptr);
  if (!state->error) {
      Py_DECREF(pModule);
      return nullptr;
    }
  
  return pModule;
}
