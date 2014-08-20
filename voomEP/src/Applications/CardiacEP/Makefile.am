bin_PROGRAMS    = EPSimulator ComputeECG
voomdir = ../../
INCLUDES        =-I$(voomdir)/       			 \
		-I$(voomdir)/Body/			 \
		-I$(voomdir)/Shape/          		 \
		-I$(voomdir)/VoomMath/  	    	 \
		-I$(voomdir)/Quadrature/		 \
		-I$(voomdir)/Mesh/			 \
		-I$(voomdir)/IonicModel/		 \
		-I$(voomdir)/Solvers			 \
		-I$(voomdir)/Material			 \
		-I$(blitz_includes) 			 \
		-I$(tvmet_includes) 			 \
		-I/u/local/apps/trilinos/current/include
ComputeECG_SOURCES      = ComputeECG.cc
EPSimulator_SOURCES         = EPSimulator.cc
AM_LDFLAGS    = -L$(blitz_libraries) -L$(xerces_libraries)		\
        -L$(voomdir)/ -L$(voomdir)/Quadrature/      			\
        -L$(voomdir)/Shape/						\
	-L$(voomdir)/VoomMath/						\
	-L$(voomdir)/Mesh/ -L$(voomdir)/Body/ -L$(voomdir)/IonicModel/ 	\
	-L$(voomdir)/Solvers -L$(voomdir)/Material			\
	-L/u/local/apps/mumps/4.10.0/intel-11.1/openmpi-1.4/lib 	\
	-L/u/local/apps/trilinos/current/lib 				\
	-L/u/local/compilers/intel/current/current/mkl/lib/em64t -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_sequential	-lmkl_core -lguide -lmkl_blacs_openmpi_lp64 -lpthread 						\
	-L/u/local/apps/superlu_dist/3.1/intel-11.1/openmpi-1.4/lib 	\
	-L/u/local/apps/parmetis/4.0.2/intel-11.1/openmpi-1.4/lib 	\
	-L/u/local/compilers/intel/11.1/080/mkl/lib/em64t
LDADD = -lblitz -lEPSolver -lCardiacBody -lIonicModel 	\
	-lMesh -lShape -lQuadrature  -lVoomMath		\
	-lMaterial -lbelos -lifpack -lamesos -laztecoo 	\
	-lteuchos 					\
	-lepetra					\
	-lepetraext 					\
	-ltriutils					\
	-ldmumps -lmumps_common -lpord -lbelosepetra	\
	-lsuperlu_dist_3.1				\
	-lparmetis -lmetis				\
	-lmkl_intel_lp64 -lmkl_sequential 		\
	-lmkl_core				\
	-lpthread -lmpi_f90 -lmpi_f77 -lgfortranbegin 	\
	-lgfortran 					\
	-lmkl_scalapack_lp64 -lmkl_blacs_openmpi_lp64 	\
	-lmkl_lapack  					\
	-lifcore                        		\
        -limf

