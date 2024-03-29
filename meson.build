project(
  'py_tsyganenko',
  'cpp', 'c', 'fortran',
  version : '0.1.6',
  default_options : ['warning_level=2', 'cpp_std=c++17', 'fortran_std=legacy'],
  license : 'GPL3'
)

pybind11_dep = dependency('pybind11', required : true)

pymod = import('python')
python3 = pymod.find_installation('python3')

extra_files = [ 'pyproject.toml', 'README.md',
                '.github/workflows/ccpp.yml','.github/workflows/tests-win.yml', '.github/workflows/pythonpublish-linux.yml',
                '.github/workflows/pythonpublish-osx.yml', '.github/workflows/pythonpublish-win.yml']

tsyganenko = static_library('tsyganenko','src/Geopack-2008_dp.for','src/T96.for',
                        extra_files: extra_files,
                        fortran_args:['-fdollar-ok', '-ffixed-form', '-fdefault-double-8', '-fdefault-real-8', '-static-libgfortran'],
                        link_args: ['-static-libgfortran'])

libquadmath = meson.get_compiler('fortran').find_library('quadmath', static:true)

srcs = [
  'wrapper/wrapper.cpp'
  ]

if build_machine.system() == 'windows'
    link_args = ['-static-libgfortran', '-static-libstdc++','-static-libgcc','-static']
elif build_machine.system() == 'darwin'
    link_args = ['-static-libgfortran', '-static-libstdc++','-static-libgcc']
else
    link_args = ['-static-libgfortran']
endif

python3.extension_module('py_tsyganenko', srcs,
                         link_with:[tsyganenko],
                         include_directories:'src',
                         dependencies: [ pybind11_dep, python3.dependency(),libquadmath],
                         link_args: link_args,
                         install: true
                        )

