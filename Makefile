#
# $Id: Makefile,v 3.31 2016/05/09 by RY$
#

## default options
AUX= ./Tools
GOURMET_HOME_PATH = /usr/local/OCTA8/GOURMET
ARCH   = linux_64
CC     = gcc
CXX    = g++
CCOPT  = -O
LINKS  = -lm -lplatform -lstdc++
GOURMET_LIB_PATH = $(GOURMET_HOME_PATH)/lib/$(ARCH)
GOURMET_INCLUDE_PATH = $(GOURMET_HOME_PATH)/include

OSTYPE = $(shell uname)
#ENV = GCC
#ENV = ICC
#ENV = ICC_MKL_OMP
#ENV = MINGW
#ENV = CYGWIN

## options for GCC/CYGWIN
ifeq ($(ENV), CYGWIN)
#ifneq (,$(findstring CYGWIN,$(OSTYPE)))
      ARCH   = cygwin
      CC     = gcc 
      CXX    = g++ 
      CCOPT  = -O3 -fno-inline
      LINKS  = -static -lm -lplatform 
endif

## options for MINGW32
ifeq ($(ENV), MINGW)
      ARCH   = win32
      CC     = i686-w64-mingw32-gcc
      CXX    = i686-w64-mingw32-g++
      CCOPT  = -O3 -fno-inline
      LINKS  = -static -lm -lplatform 
endif

## options for MINGW64
ifeq ($(ENV), MINGW64)
      ARCH   = win64
      CC     = x86_64-w64-mingw32-gcc
      CXX    = x86_64-w64-mingw32-g++
      CCOPT  = -O3 -fno-inline
      LINKS  = -static -lm -lplatform 
endif

ifeq ($(ENV), CLANG)
     ARCH    = macosx
     CC	     = clang
     CXX     = clang++
     LINKS   = -L/usr/local/lib -lm -lplatform -stdlib=libc++
     CCOPT   = -I/usr/local/include -g -fcolor-diagnostics -stdlib=libc++ 
endif

## options for GCC/LINUX
ifeq ($(ENV), GCC)
      ARCH   = linux_64
      CC     = gcc
      CXX    = g++
      CCOPT  = -O3 
      LINKS  = -lm -lplatform -lstdc++
endif

## options for ICC/LINUX
ifeq ($(ENV), ICC)
      ARCH   = linux_64
      CC     = icc 
      CXX    = icpc 
      CCOPT  = -O3 -xSSSE3 -axAVX,SSE4.2,SSE4.1,SSSE3,SSE3,SSE2 -w0
#      LINKS  = -lm -lplatform -lcxaguard -lstdc++
      LINKS  = -lm -lplatform -lstdc++
endif

## options for GCC+MKL+OMP/LINUX
ifeq ($(ENV), ICC_MKL_OMP)
      ARCH   = linux_64
      CC     = icc 
      CXX    = icpc 
      CCOPT  = -O3 -xSSSE3 -axAVX,SSE4.2,SSE4.1,SSSE3,SSE3,SSE2\
	-ip -qopenmp -parallel -w0
#      LINKS  = -lplatform -lcxaguard -lstdc++\
#	-lmkl_intel_lp64 -lmkl_intel_thread  -lmkl_core -lm
      LINKS  = -lplatform -lstdc++\
	-lmkl_intel_lp64 -lmkl_intel_thread  -lmkl_core -lm
endif

OBJS  	= mt19937ar.o\
	operate_electrolyte.o\
	fluct.o\
	alloc.o\
	solute_rhs.o\
	fftsg.o\
	fftsg3d.o\
	avs_output.o\
	avs_output_p.o\
	output.o\
	resume.o\
	make_phi.o\
	fluid_solver.o\
	particle_solver.o\
	md_force.o\
	profile.o\
	interaction.o\
	operate_omega.o\
	fft_wrapper.o\
	f_particle.o\
	init_fluid.o\
	init_particle.o\
	input.o\
	rigid_body.o\
	operate_surface.o\
	matrix_diagonal.o\
	periodic_spline.o\
	sp_3d_ns.o

ifeq ($(HDF5), ON)
      LINKS  += -lhdf5 -lhdf5_hl
      CCOPT  += -DWITH_EXTOUT -L/opt/hdf5.1.8/lib -I/opt/hdf5.1.8/include
      OBJS   += output_writer.o
endif
CFLAGS 	= $(CCOPT) -I$(GOURMET_INCLUDE_PATH)
LINKS  += -L$(GOURMET_LIB_PATH) 

XYZ_OBJS= alloc.o\
	rigid_body.o\
	$(AUX)/udf2xyz.o

TARGET 	= kapsel 
XYZ	= udf2xyz

## Implicit rules

.SUFFIXES: .c .cxx .o .out

## Build rules

all: $(TARGET) $(XYZ)

$(TARGET): $(OBJS)
	$(CXX) $(OBJS) -o $(TARGET) $(CFLAGS) $(LINKS)

$(XYZ): $(XYZ_OBJS)
	$(CXX) $(XYZ_OBJS) -o $(XYZ) $(CFLAGS) $(LINKS)

## Compile

.cxx.o: 
	$(CXX) -c $< $(CFLAGS) -o $@

.c.o: 
	$(CC) -c $< $(CFLAGS) -o $@

## Clean

clean:
	rm -f $(OBJS) $(AUX)/$(XYZ_OBJS) $(TARGET) $(XYZ)
	rm -f *~ *.bak

depend:
	makedepend -- $(CFLAGS) -- *.cxx *.c *.h
