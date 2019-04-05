#
# $Id: Makefile,v 3.31 2016/05/09 by RY$
#
### Command line options for make ###
#
### FOR LINUX ##
#ENV = GCC
#ENV = ICC
#ENV = ICC_OMP
#ENV = ICC_MKL_OMP
ENV = ICC_MKL_MPI
#ENV = ICC_MKL_OMP_MPI

#
### FOR WINDOWS ###
#ENV = CYGWIN
#ENV = MINGW64
#ENV = MINGW
#
### FOR MAC ###
#ENV = GCC_MAC
#ENV = CLANG
#
### FFT LIB ###
#FFT = FFTW
#FFT = IMKL
#FFT = OOURA
### HDF5 SUPPORT ###
#HDF5 = ON

## default options
# Use OCTA environment variables
#GOURMET_HOME_PATH = $(PF_FILES)
#ENGINE_HOME_PATH  = $(PF_ENGINE)
#ARCH              = $(PF_ENGINEARCH)
# OR
# Define environment variables explicitly here
GOURMET_HOME_PATH  = /home/OCTA81/GOURMET
ENGINE_HOME_PATH   = /home/OCTA81/ENGINES
ARCH               = linux_64
#
## FFTW options
FFTW_LIB_PATH = /opt/fftw/3.3.7/lib
FFTW_INCLUDE_PATH = /opt/fftw/3.3.7/include

AUX= ./Tools
CC     = gcc
CXX    = g++
CCOPT  = -O
LINKS  = -lm -lplatform -lstdc++
GOURMET_LIB_PATH = $(GOURMET_HOME_PATH)/lib/$(ARCH)
GOURMET_INCLUDE_PATH = $(GOURMET_HOME_PATH)/include
TARGET_DIR=$(ENGINE_HOME_PATH)/bin/$(ARCH)
OSTYPE = $(shell uname)

ifeq ($(ENV),)
  ifneq (,$(findstring CYGWIN,$(OSTYPE)))
    ifeq ($(ARCH),win32)
      CC     = i686-w64-mingw32-gcc
      CXX    = i686-w64-mingw32-g++
      CCOPT  = -O3 -fno-inline
      LINKS  = -lm -lplatform -static
    else
      ifeq ($(ARCH),win64)
        CC     = x86_64-w64-mingw32-gcc
        CXX    = x86_64-w64-mingw32-g++
        CCOPT  = -O3 -fno-inline
        LINKS  = -lm -lplatform -static
      else
        ARCH   = cygwin
        CCOPT  = -O3 -fno-inline
        LINKS  = -lm -lplatform -static
      endif
    endif
  endif
  ifeq ($(OSTYPE), Linux)
    ARCH   = linux_64
    CCOPT  = -O3 
    LINKS  = -lm -lplatform -lstdc++ -static
  endif
endif

## options for GCC/CYGWIN/WINDOWS (not supported)
ifeq ($(ENV), CYGWIN)
#ifneq (,$(findstring CYGWIN,$(OSTYPE)))
      ARCH   = cygwin
      CC     = gcc 
      CXX    = g++ 
      CCOPT  = -O3 -fno-inline
      LINKS  = -lm -lplatform 
endif

## options for MINGW32/CYGWIN/WINDOWS (not supported)
ifeq ($(ENV), MINGW)
      ARCH   = win32
      CC     = i686-w64-mingw32-gcc
      CXX    = i686-w64-mingw32-g++
      CCOPT  = -O3 -fno-inline
      LINKS  = -static -lm -lplatform 
endif

## options for MINGW64/CYGWIN/WINDOWS (not supported)
ifeq ($(ENV), MINGW64)
      ARCH   = win64
      CC     = x86_64-w64-mingw32-gcc
      CXX    = x86_64-w64-mingw32-g++
      CCOPT  = -O3 -fno-inline
      LINKS  = -static -lm -lplatform 
endif

## options for CLANG/MAC (not supported)
ifeq ($(ENV), CLANG)
     ARCH    = macosx
     CC	     = clang
     CXX     = clang++
     LINKS   = -L/usr/local/lib -lm -lplatform -stdlib=libc++
     CCOPT   = -I/usr/local/include -g -fcolor-diagnostics -stdlib=libc++ 
endif

## options for GCC/MAC (not supported)
ifeq ($(ENV), GCC_MAC)
     ARCH    = macosx
     CC	     = gcc-5
     CXX     = g++-5
     CCOPT  = -I/usr/local/include -O3 -fno-inline -std=c++11
     LINKS  = -L/usr/local/lib -lm -lplatform_gcc-5 
endif

## options for GCC/LINUX (not supported)
ifeq ($(ENV), GCC)
      ARCH   = linux_64
      CC     = gcc
      CXX    = g++
      CCOPT  = -O3 -std=c++11
      LINKS  = -lm -lplatform -lstdc++ -static
	ifeq ($(HDF5), ON)
		LINKS  += -L/usr/local/hdf5/lib
		CCOPT  += -I/usr/local/hdf5/include
	endif
    ifeq ($(FFT), OOURA)
		CCOPT += -D_FFT_OOURA
    endif
endif

## options for ICC/LINUX (not supported)
ifeq ($(ENV), ICC)
      ARCH   = linux_64
      CC     = icc 
      CXX    = icpc 
      CCOPT  = -std=c++11 -O3 -xSSSE3 -axAVX,SSE4.2,SSE4.1,SSSE3,SSE3,SSE2 -Wall -D__INTEL_COMPILER
#      LINKS  = -lm -lplatform -lcxaguard -lstdc++
      LINKS  = -lm -lplatform -lstdc++ -static-intel
	ifeq ($(HDF5), ON)
		LINKS  += -L/usr/local/hdf5/lib
		CCOPT  += -I/usr/local/hdf5/include
	endif
    ifeq ($(FFT), FFTW)
		CCOPT += -I$(FFTW_INCLUDE_PATH) -D_FFT_FFTW
		LINKS += -L$(FFTW_LIB_PATH) -lfftw3
    endif
    ifeq ($(FFT), IMKL)
		CCOPT += -D_FFT_IMKL
		LINKS += -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread
    endif
    ifeq ($(FFT), OOURA)
		CCOPT += -D_FFT_OOURA
    endif
endif

## options for ICC/LINUX with OMP (not supported)
ifeq ($(ENV), ICC_OMP)
     ARCH   = linux_64
     CC     = icc 
     CXX    = icpc 
     CCOPT  = -std=c++11 -O3 -xSSSE3 -axCOMMON-AVX512,CORE-AVX512,CORE-AVX2,CORE-AVX-I,AVX,SSE4.2,SSE4.1,SSSE3,SSE3,SSE2 -ip -qopenmp -parallel -w0 -D_OPENMP -DTIME_MEASURE
     LINKS  = -lm -lplatform -lstdc++
     ifeq ($(HDF5), ON)
	LINKS  += -L/opt/hdf5.1.8/lib
	CCOPT  += -I/opt/hdf5.1.8/include
     endif
     ifeq ($(FFT), FFTW)
	CCOPT += -I$(FFTW_INCLUDE_PATH) -D_FFT_FFTW
	#LINKS += -L$(FFTW_LIB_PATH) -lfftw3_threads -lfftw3
	LINKS += -L$(FFTW_LIB_PATH) -lfftw3_omp -lfftw3 
     endif
     ifeq ($(FFT), IMKL)
	CCOPT += -D_FFT_IMKL
	LINKS += -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lpthread -static_intel 
     endif
endif

## options for icc+MKL+OMP/LINUX (not supported)
ifeq ($(ENV), ICC_MKL_OMP)
      ARCH   = linux_64
      CC     = mpiicc 
      CXX    = mpiicpc 
      CCOPT  = -O3 -xSSSE3 -axAVX,SSE4.2,SSE4.1,SSSE3,SSE3,SSE2 -ip -parallel -w0
#      LINKS  = -lplatform -lcxaguard -lstdc++ -lmkl_intel_lp64 -lmkl_intel_thread  -lmkl_core -lm
      LINKS  = -lplatform -lstdc++ -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lm -static-intel
	ifeq ($(HDF5), ON)
		LINKS  += -L/usr/local/hdf5/lib
		CCOPT  += -I/usr/local/hdf5/include
	endif
endif

## options for mpicc+MKL/LINUX
ifeq ($(ENV), ICC_MKL_MPI)
      ARCH   = linux_64
      CC     = mpiicc 
      CXX    = mpiicpc 
      CCOPT  = -std=c++11 -O3 -xSSSE3 -axAVX,SSE4.2,SSE4.1,SSSE3,SSE3,SSE2 -ip -parallel -w0 -D_FFT_MPI_IMKL -D_MPI
      LINKS  = -lplatform -lstdc++ -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lm -static-intel
	ifeq ($(HDF5), ON)
		LINKS  += -L/usr/local/hdf5/lib
		CCOPT  += -I/usr/local/hdf5/include
	endif
endif

## options for mpicc+MKL+OMP/LINUX (not supported)
ifeq ($(ENV), ICC_MKL_OMP_MPI)
      ARCH   = linux_64
      CC     = mpiicc 
      CXX    = mpiicpc 
      CCOPT  = -std=c++11 -O3 -xSSSE3 -axAVX,SSE4.2,SSE4.1,SSSE3,SSE3,SSE2 -ip -qopenmp -parallel -w0 -D_FFT_MPI_IMKL -D_MPI -D_OPENMP -DTIME_MEASURE
      LINKS  = -lplatform -lstdc++ -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lm -static-intel
	ifeq ($(HDF5), ON)
		LINKS  += -L/usr/local/hdf5/lib
		CCOPT  += -I/usr/local/hdf5/include
	endif
endif
OBJS  	= \
	operate_mpi.o\
	operate_mpi_particle.o\
	fft_wrapper_base.o\
	memory_model.o\
	variable.o\
	operate_electrolyte.o\
	fluct.o\
	alloc.o\
	solute_rhs.o\
	fftsg.o\
	fftsg3d.o\
	avs_output.o\
	avs_output_p.o\
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
	matrix_diagonal.o \
	operate_surface.o \
	output.o \
	periodic_spline.o \
	rigid_body.o \
	sp_3d_ns.o \
	mt19937ar.o

## options for HDF5 support
ifeq ($(HDF5), ON)
      LINKS  += -lhdf5 -lhdf5_hl
      CCOPT  += -DWITH_EXTOUT
      OBJS   += output_writer.o
endif

#CXXFLAGS = $(CCOPT) -I$(GOURMET_INCLUDE_PATH) -D_MPI -D_OPENMP -DNDEBUG -D__INTEL_COMPILER
CXXFLAGS = $(CCOPT) -I$(GOURMET_INCLUDE_PATH) -DNDEBUG


LINKS  += -L$(GOURMET_LIB_PATH) 

XYZ_OBJS= alloc.o\
	rigid_body.o\
	$(AUX)/udf2xyz.o

TARGET 	= kapsel_u2m_mpi
XYZ	= udf2xyz

ENGINE = $(TARGET)
CONVERTER =  $(XYZ)
ifeq ($(ARCH), cygwin)
  ENGINE = $(TARGET).exe
  CONVERTER = $(XYZ).exe
endif
ifeq ($(ARCH), win32)
  ENGINE = $(TARGET).exe
  CONVERTER = $(XYZ).exe
endif
ifeq ($(ARCH), win64)
  ENGINE = $(TARGET).exe
  CONVERTER = $(XYZ).exe
endif

## Implicit rules

.SUFFIXES: .c .cxx .o .out

## Build rules

all: $(TARGET) $(XYZ)

$(TARGET): $(OBJS)
	$(CXX) $(OBJS) -o $(TARGET) $(CXXFLAGS) $(LINKS)

$(XYZ): $(XYZ_OBJS)
	$(CXX) $(XYZ_OBJS) -o $(XYZ) $(CXXFLAGS) $(LINKS)

## Compile

.cxx.o: 
	$(CXX) -c $< $(CXXFLAGS) -o $@

.c.o: 
	$(CC) -c $< $(CXXFLAGS) -o $@

## Clean

clean:
	rm -f $(OBJS) $(AUX)/$(XYZ_OBJS) $(TARGET) $(XYZ)
	rm -f *~ *.bak

## Install

install:
	if [ x$(ENGINE_HOME_PATH) = x ]; then \
		echo " Environment variable PF_ENIGNE must be defined to install engine."; \
		exit 1;\
	fi
	if [ ! -d $(TARGET_DIR) ]; then \
		mkdir -p $(TARGET_DIR) ; \
	fi;\
	if [ -f $(ENGINE) ]; then \
		cp -f $(ENGINE) $(TARGET_DIR)/ ;\
	fi
	if [ -f $(CONVERTER) ]; then \
		cp -f $(CONVERTER) $(TARGET_DIR)/ ;\
	fi

depend:
	makedepend -- $(CXXFLAGS) -- *.cxx *.c *.h
