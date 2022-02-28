project(
  'py_tsyganenko',
  'cpp', 'c', 'fortran',
  version : '0.1.3',
  default_options : ['warning_level=2', 'cpp_std=c++17', 'fortran_std=legacy'],
  license : 'GPL3'
)

pybind11_dep = dependency('pybind11', required : true)

pymod = import('python')
python3 = pymod.find_installation('python3')

tsyganenko = library('tsyganenko','src/Geopack-2008_dp.for','src/T96.for',
                        fortran_args:['-fdollar-ok', '-ffixed-form', '-fdefault-double-8', '-fdefault-real-8'])

srcs = [
  'wrapper/wrapper.cpp'
  ]

python3.extension_module('py_tsyganenko', srcs,
                         link_with:tsyganenko,
                         include_directories:'src',
                         dependencies: [ pybind11_dep, python3.dependency() ],
                         install: true
                        )

