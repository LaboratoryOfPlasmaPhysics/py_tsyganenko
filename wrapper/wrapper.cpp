#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/chrono.h>

namespace py = pybind11;

extern "C" {
void t96_01_(int IOPT, double PARMOD[10], double *PS, double *X, double *Y,
             double *Z, double *BX, double *BY, double *BZ);
void recalc_08_(int *IYEAR, int *IDAY, int *IHOUR, int *MIN, int *ISEC,
                double *VGSEX, double *VGSEY, double *VGSEZ);

void igrf_geo_08_(double *R, double *THETA, double *PHI, double *BR,
                  double *BTHETA, double *BPHI);
};

auto T96(double PDYN, double DST, double BYIMF, double BZIMF, double PS,
         double X, double Y, double Z) {
  double BX, BY, BZ;
  double PARMOD[10] = {PDYN, DST, BYIMF, BZIMF, 0.};
  t96_01_(1, PARMOD, &PS, &X, &Y, &Z, &BX, &BY, &BZ);
  return std::tuple{BX, BY, BZ};
}

void recalc(int year, int day, int hour, int minute, int second, double VGSEX,
            double VGSEY, double VGSEZ) {
  recalc_08_(&year, &day, &hour, &minute, &second, &VGSEX, &VGSEY, &VGSEZ);
}

auto igrf_geo(double R, double THETA, double PHI) {
  double BR, BTHETA, BPHI;
  igrf_geo_08_(&R, &THETA, &PHI, &BR, &BTHETA, &BPHI);
  return std::tuple{BR, BTHETA, BPHI};
}

PYBIND11_MODULE(py_tsyganenko, m) {
  m.doc() = "py_tsyganenko module";
  m.def("T96", T96, py::arg("PDYN"), py::arg("DST"), py::arg("BYIMF"),
        py::arg("BZIMF"), py::arg("PS"), py::arg("X"), py::arg("Y"),
        py::arg("Z"));
  m.def("recalc", recalc, py::arg("year"), py::arg("day"), py::arg("hour"),
        py::arg("minute"), py::arg("second"), py::arg("VGSEX"),
        py::arg("VGSEY"), py::arg("VGSEZ"));
  m.def("igrf_geo", igrf_geo, py::arg("R"), py::arg("THETA"), py::arg("PHI"));
}
