## Makefile.am -- Process this file with automake to produce Makefile.in
AM_CPPFLAGS= -I$(srcdir)/../VoomMath -I$(srcdir)/.. -I$(blitz_includes) -I$(tvmet_includes) -I$(srcdir)/../Node -I$(srcdir)/../Quadrature
lib_LIBRARIES=libShape.a
libShape_a_SOURCES = LoopShellShape.cc ShapeTri3.cc ShapeTri6.cc ShapeQ4.cc \
		     ShapeTet4CP.cc ShapeTet10.cc ShapeHex8.cc Hermite.cc \
		     ShapeBar.cc

