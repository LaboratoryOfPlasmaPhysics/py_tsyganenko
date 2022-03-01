#include <pybind11/chrono.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

extern "C"
{
    void t96_01_(int IOPT, double PARMOD[10], double* PS, double* X, double* Y, double* Z,
        double* BX, double* BY, double* BZ);
    void recalc_08_(int* IYEAR, int* IDAY, int* IHOUR, int* MIN, int* ISEC, double* VGSEX,
        double* VGSEY, double* VGSEZ);

    void igrf_geo_08_(
        double* R, double* THETA, double* PHI, double* BR, double* BTHETA, double* BPHI);

    void igrf_gsw_08_(
        double* XGSW, double* YGSW, double* ZGSW, double* HXGSW, double* HYGSW, double* HZGSW);

    void geomag_08_(
        double* XGEO, double* YGEO, double* ZGEO, double* XMAG, double* YMAG, double* ZMAG, int* J);

    void gswgse_08_(
        double* XGSW, double* YGSW, double* ZGSW, double* XGSE, double* YGSE, double* ZGSE, int* J);

    void sphcar_08_(double* R, double* THETA, double* PHI, double* X, double* Y, double* Z, int* J);
    void bcarsp_08_(double* X, double* Y, double* Z, double* BX, double* BY, double* BZ, double* BR,
        double* BTHETA, double* BPHI);

    void bspcar_08_(double* THETA, double* PHI, double* BR, double* BTHETA, double* BPHI,
        double* BX, double* BY, double* BZ);
};

template <typename T>
struct np_coord_change
{
    T f;
    np_coord_change(T f) : f { f } {};

    auto operator()(const py::array_t<double>& COORDS)
    {
        py::buffer_info in_buff = COORDS.request();
        if (in_buff.ndim != 2)
            throw std::runtime_error("Number of dimensions must be 2");
        if (in_buff.shape[1] != 3)
            throw std::runtime_error("POS vector shape must be [3:x]");
        auto result = py::array_t<double>(in_buff.shape);
        auto in = COORDS.unchecked<2>();
        auto out = result.mutable_unchecked<2>();
        for (py::ssize_t i = 0; i < in_buff.shape[0]; i++)
        {
            auto [X, Y, Z] = f(in(i, 0), in(i, 1), in(i, 2));
            out(i, 0) = X;
            out(i, 1) = Y;
            out(i, 2) = Z;
        }
        return result;
    }
};
auto T96(
    double PDYN, double DST, double BYIMF, double BZIMF, double PS, double X, double Y, double Z)
{
    double BX, BY, BZ;
    double PARMOD[10] = { PDYN, DST, BYIMF, BZIMF, 0. };
    t96_01_(1, PARMOD, &PS, &X, &Y, &Z, &BX, &BY, &BZ);
    return std::tuple { BX, BY, BZ };
}

py::array_t<double> T96_v(
    double PDYN, double DST, double BYIMF, double BZIMF, double PS, py::array_t<double> POS)
{
    py::buffer_info in_buff = POS.request();
    if (in_buff.ndim != 2)
        throw std::runtime_error("Number of dimensions must be 2");
    if (in_buff.shape[1] != 3)
        throw std::runtime_error("POS vector shape must be [3:x]");
    auto result = py::array_t<double>(in_buff.shape);
    auto in = POS.unchecked<2>();
    auto out = result.mutable_unchecked<2>();
    for (py::ssize_t i = 0; i < in_buff.shape[0]; i++)
    {
        auto [BX, BY, BZ] = T96(PDYN, DST, BYIMF, BZIMF, PS, in(i, 0), in(i, 1), in(i, 2));
        out(i, 0) = BX;
        out(i, 1) = BY;
        out(i, 2) = BZ;
    }

    return result;
}

void recalc(
    int year, int day, int hour, int minute, int second, double VGSEX, double VGSEY, double VGSEZ)
{
    recalc_08_(&year, &day, &hour, &minute, &second, &VGSEX, &VGSEY, &VGSEZ);
}

auto igrf_geo(double R, double THETA, double PHI)
{
    double BR, BTHETA, BPHI;
    igrf_geo_08_(&R, &THETA, &PHI, &BR, &BTHETA, &BPHI);
    return std::tuple { BR, BTHETA, BPHI };
}

auto igrf_gsw(double XGSW, double YGSW, double ZGSW)
{
    double HXGSW, HYGSW, HZGSW;
    igrf_gsw_08_(&XGSW, &YGSW, &ZGSW, &HXGSW, &HYGSW, &HZGSW);
    return std::tuple { HXGSW, HYGSW, HZGSW };
}

auto geo_to_mag(double XGEO, double YGEO, double ZGEO)
{
    int J = 1;
    double XMAG, YMAG, ZMAG;
    geomag_08_(&XGEO, &YGEO, &ZGEO, &XMAG, &YMAG, &ZMAG, &J);
    return std::tuple { XMAG, YMAG, ZMAG };
}

auto mag_to_geo(double XMAG, double YMAG, double ZMAG)
{
    int J = -1;
    double XGEO, YGEO, ZGEO;
    geomag_08_(&XGEO, &YGEO, &ZGEO, &XMAG, &YMAG, &ZMAG, &J);
    return std::tuple { XGEO, YGEO, ZGEO };
}

auto sph_to_cart(double R, double THETA, double PHI)
{
    int J = 1;
    double X, Y, Z;
    sphcar_08_(&R, &THETA, &PHI, &X, &Y, &Z, &J);
    return std::tuple { X, Y, Z };
}

auto cart_to_sph(double X, double Y, double Z)
{
    int J = -1;
    double R, THETA, PHI;
    sphcar_08_(&R, &THETA, &PHI, &X, &Y, &Z, &J);
    return std::tuple { R, THETA, PHI };
}

auto B_cart_to_sph(double X, double Y, double Z, double BX, double BY, double BZ)
{
    double BR, BTHETA, BPHI;
    bcarsp_08_(&X, &Y, &Z, &BX, &BY, &BZ, &BR, &BTHETA, &BPHI);
    return std::tuple { BR, BTHETA, BPHI };
}

auto B_sph_to_cart(double THETA, double PHI, double BR, double BTHETA, double BPHI)
{
    double BX, BY, BZ;
    bspcar_08_(&THETA, &PHI, &BR, &BTHETA, &BPHI, &BX, &BY, &BZ);
    return std::tuple { BX, BY, BZ };
}

py::array_t<double> B_cart_to_sph_v(py::array_t<double> COORDS, py::array_t<double> B)
{
    py::buffer_info in_coords_buff = COORDS.request();
    py::buffer_info in_b_buff = B.request();
    if (in_coords_buff.ndim != 2 && in_b_buff.ndim != 2)
        throw std::runtime_error("Number of dimensions must be 2");
    if (in_coords_buff.shape[1] != 3 && in_b_buff.shape[1] != 3)
        throw std::runtime_error("COORDS and B vectors shape must be [3:x]");
    if (in_coords_buff.shape[0] == in_b_buff.shape[0])
        throw std::runtime_error("COORDS and B vectors shape must have same length");
    auto result = py::array_t<double>(in_coords_buff.shape);
    auto in_coords = COORDS.unchecked<2>();
    auto in_b = B.unchecked<2>();
    auto out = result.mutable_unchecked<2>();
    for (py::ssize_t i = 0; i < in_coords_buff.shape[0]; i++)
    {
        auto [BR, BTHETA, BPHI] = B_cart_to_sph(
            in_coords(i, 0), in_coords(i, 1), in_coords(i, 2), in_b(i, 0), in_b(i, 1), in_b(i, 2));
        out(i, 0) = BR;
        out(i, 1) = BTHETA;
        out(i, 2) = BPHI;
    }

    return result;
}

py::array_t<double> B_sph_to_cart_v(py::array_t<double> COORDS, py::array_t<double> B)
{
    py::buffer_info in_coords_buff = COORDS.request();
    py::buffer_info in_b_buff = B.request();
    if (in_coords_buff.ndim != 2 && in_b_buff.ndim != 2)
        throw std::runtime_error("Number of dimensions must be 2");
    if (in_coords_buff.shape[1] != 2)
        throw std::runtime_error("COORDS vector shape must be [2:x]");
    if (in_b_buff.shape[1] != 3)
        throw std::runtime_error("B vector shape must be [3:x]");
    if (in_coords_buff.shape[0] == in_b_buff.shape[0])
        throw std::runtime_error("COORDS and B vectors shape must have same length");
    auto result = py::array_t<double>(in_b_buff.shape);
    auto in_coords = COORDS.unchecked<2>();
    auto in_b = B.unchecked<2>();
    auto out = result.mutable_unchecked<2>();
    for (py::ssize_t i = 0; i < in_coords_buff.shape[0]; i++)
    {
        auto [BX, BY, BZ]
            = B_sph_to_cart(in_coords(i, 0), in_coords(i, 1), in_b(i, 0), in_b(i, 1), in_b(i, 2));
        out(i, 0) = BX;
        out(i, 1) = BY;
        out(i, 2) = BZ;
    }
    return result;
}

auto gsw_to_gse(double XGSW, double YGSW, double ZGSW)
{
    int J = 1;
    double XGSE, YGSE, ZGSE;
    gswgse_08_(&XGSW, &YGSW, &ZGSW, &XGSE, &YGSE, &ZGSE, &J);
    return std::tuple { XGSE, YGSE, ZGSE };
}

auto gse_to_gsw(double XGSE, double YGSE, double ZGSE)
{
    int J = -1;
    double XGSW, YGSW, ZGSW;
    gswgse_08_(&XGSW, &YGSW, &ZGSW, &XGSE, &YGSE, &ZGSE, &J);
    return std::tuple { XGSW, YGSW, ZGSW };
}

PYBIND11_MODULE(py_tsyganenko, m)
{
    m.doc()
        = "py_tsyganenko module wrapper for N. A. Tsyganenko Fortran codes from "
          "https://geo.phys.spbu.ru/~tsyganenko/empirical-models/";
    auto geopack = m.def_submodule("Geopack",
        "Geopack module wrapper, see "
        "https://geo.phys.spbu.ru/~tsyganenko/empirical-models/coordinate_systems/geopack");

    geopack.def("recalc", recalc, py::arg("year"), py::arg("day"), py::arg("hour"),
        py::arg("minute"), py::arg("second"), py::arg("VGSEX"), py::arg("VGSEY"), py::arg("VGSEZ"));

    geopack.def("igrf_geo", igrf_geo, py::arg("R"), py::arg("THETA"), py::arg("PHI"));
    geopack.def("igrf_geo", np_coord_change(igrf_geo), py::arg("GEO_COORDS"));

    geopack.def("igrf_gsw", igrf_gsw, py::arg("XGSW"), py::arg("YGSW"), py::arg("ZGSW"));
    geopack.def("igrf_gsw", np_coord_change(igrf_gsw), py::arg("GSW_COORDS"));


    geopack.def("geo_to_mag", geo_to_mag, py::arg("XGEO"), py::arg("YGEO"), py::arg("ZGEO"));
    geopack.def("mag_to_geo", mag_to_geo, py::arg("XMAG"), py::arg("YMAG"), py::arg("ZMAG"));
    geopack.def("geo_to_mag", np_coord_change(geo_to_mag), py::arg("GEO_COORDS"));
    geopack.def("mag_to_geo", np_coord_change(mag_to_geo), py::arg("MAG_COORDS"));

    geopack.def("gsw_to_gse", gsw_to_gse, py::arg("XGSW"), py::arg("YGSW"), py::arg("ZGSW"));
    geopack.def("gse_to_gsw", gse_to_gsw, py::arg("XGSE"), py::arg("YGSE"), py::arg("ZGSE"));
    geopack.def("gsw_to_gse", np_coord_change(gsw_to_gse), py::arg("GSW_COORDS"));
    geopack.def("gse_to_gsw", np_coord_change(gse_to_gsw), py::arg("GSE_COORDS"));

    geopack.def("sph_to_cart", sph_to_cart, py::arg("R"), py::arg("THETA"), py::arg("PHI"));
    geopack.def("cart_to_sph", cart_to_sph, py::arg("X"), py::arg("Y"), py::arg("Z"));
    geopack.def("sph_to_cart", np_coord_change(sph_to_cart), py::arg("SPH_COORDS"));
    geopack.def("cart_to_sph", np_coord_change(cart_to_sph), py::arg("CART_COORDS"));

    geopack.def("B_cart_to_sph", B_cart_to_sph, py::arg("X"), py::arg("Y"), py::arg("Z"),
        py::arg("BX"), py::arg("BY"), py::arg("BZ"));
    geopack.def("B_sph_to_cart", B_sph_to_cart, py::arg("THETA"), py::arg("PHI"), py::arg("BR"),
        py::arg("BTHETA"), py::arg("BPHI"));
    geopack.def("B_cart_to_sph", B_cart_to_sph_v, py::arg("COORDS"), py::arg("B"));
    geopack.def("B_sph_to_cart", B_sph_to_cart_v, py::arg("COORDS"), py::arg("B"));

    auto models = m.def_submodule("Models",
        "Magnetic field models wrapper, see "
        "https://geo.phys.spbu.ru/~tsyganenko/empirical-models/magnetic_field/t96");

    models.def("T96", T96, py::arg("PDYN"), py::arg("DST"), py::arg("BYIMF"), py::arg("BZIMF"),
        py::arg("PS"), py::arg("X"), py::arg("Y"), py::arg("Z"));
    models.def("T96", T96_v, py::arg("PDYN"), py::arg("DST"), py::arg("BYIMF"), py::arg("BZIMF"),
        py::arg("PS"), py::arg("POS"));
}
