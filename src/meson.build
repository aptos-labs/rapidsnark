cmake = import('cmake')

# Required on macos because somehow the homebrew version of gmp has a library
# binary that supports 64-bit limbs, but gmp.h only allows this if you set this 
# flag
add_project_arguments('-D_LONG_LONG_LIMB', language : 'cpp')

gmp_dep = dependency('gmp')




# oneTBB subproject/dependency
opt_var = cmake.subproject_options()
opt_var.add_cmake_defines({'TBB_TEST': false})
tbb_dep = cmake.subproject('oneTBB').dependency('tbb')

deps = [gmp_dep, tbb_dep]


# I'm using `includes` for header-only libs
incdirs = include_directories('../includes')

src = [ 
  'alt_bn128.cpp',
  #'alt_bn128_test.cpp',
  'binfile_utils.cpp',
  'curve.cpp',
  'f2field.cpp',
  'fft.cpp',
  'fq.cpp',
  'fq_generic.cpp',
  'fq_raw_generic.cpp',
  'fr.cpp',
  'fr_generic.cpp',
  'fr_raw_generic.cpp',
  'fullprover.cpp',
  'groth16.cpp',
  'logger.cpp',
  'main.cpp',
  'misc.cpp',
  'multiexp.cpp',
  'naf.cpp',
  'scalar.cpp',
  'splitparstr.cpp',
  #'splitparstr_test.cpp',
  #'test_prover.cpp',
  ]

rapidsnark_lib = static_library(
  'rapidsnark', 
  src, 
  dependencies: deps,
  include_directories: incdirs
)

