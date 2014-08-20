## Makefile.am -- Process this file with automake to produce Makefile.in
INCLUDES     =  -I$(srcdir)/				\
		-I$(srcdir)/../				\
                -I$(srcdir)/../Body/	                \
                -I$(srcdir)/../IonicModel/              \
		-I$(srcdir)/../Mesh			\
		-I$(srcdir)/../Material			\
		-I$(srcdir)/../VoomMath			\
                -I$(blitz_includes)                     \
                -I$(tvmet_includes)                     \
                -I/u/local/apps/trilinos/current/include
AM_LDFLAGS    = -L$(blitz_libraries)                    \
	        -L../	                                \
	        -L../Mesh/                              \
	        -L/u/local/apps/mumps/4.9.2/lib                 \
	        -L/u/local/apps/trilinos/current/lib     \
	        -L/u/local/apps/scalapack/current               \
	        -L/u/local/apps/superlu_dist/2.4/lib            \
	        -L/u/local/apps/parmetis/3.1.1                  \
	        -L/u/local/compilers/intel/current/current/mkl/lib/em64t
LDADD         = -lblitz                                    \
	        -lQuadrature                                    \
	        -lShape                                         \
	        -lVoomMath -lMaterial

lib_LIBRARIES=libEPSolver.a
libEPSolver_a_SOURCES = EPSolver.cc EP_ExplicitSolver.cc EP_BackwardEuler.cc \
			EP_CrankNicholson.cc

