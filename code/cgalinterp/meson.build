project('cgalinterp', 'cpp')

py = import('python').find_installation(pure: false)

lib_gmp = dependency('gmp')
lib_mpfr = dependency('mpfr')
# numpy_dep = dependency('numpy')

incdir_numpy = run_command(py, ['-c',
    'import os; os.chdir(".."); import numpy; print(numpy.get_include())'
  ],
  check : true
).stdout().strip()

inc_np = include_directories(incdir_numpy)

py.extension_module(
    'cgalinterp',
    'cgalinterp.cpp',
    include_directories: inc_np,
    dependencies: [lib_mpfr],
    install: true
)
