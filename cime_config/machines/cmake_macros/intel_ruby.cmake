string(APPEND CPPDEFS " -DNO_SHR_VMATH -DCNL")
if (DEBUG)
  string(APPEND FFLAGS " -check all -ftrapuv")
endif()
string(APPEND SLIBS " -llapack -lblas")
string(APPEND LDFLAGS " -L/usr/tce/packages/gcc/gcc-10.3.1-magic/lib/gcc/x86_64-redhat-linux/10/")
set(KOKKOS_OPTIONS "--with-serial --ldflags='-L/usr/tce/packages/gcc/gcc-10.3.1-magic/lib/gcc/x86_64-redhat-linux/10/'")
set(MPI_LIB_NAME "mpich")
set(MPI_PATH "/usr/tce/packages/mvapich2/mvapich2-2.3.7-intel-classic-2021.6.0/")