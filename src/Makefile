# Compiler
FC = gfortran
#FC = ifort

# FFT Option : Give option of which FFT to use
USE_FFT=FFTW

# FFT Option : Give option of which FFT library to use
USE_LAPACKLIB=LAPACK

# Flags for debugging 
#FCFLAGS = -g -fbounds-check -fcheck=all -ffpe-trap=invalid,zero,overflow -Wall \
          -fcheck=all -ffpe-trap=invalid,zero,overflow -D$(USE_FFT) -fopenmp 
#FCFLAGS = -O3 -xHost -ipo#-march=native 
#FCFLAGS = -O3 -axAVX 
#FCFLAGS = -O0 -D$(USE_FFT) -fopenmp # checked - Ok   
FCFLAGS = -O3 -D$(USE_FFT) -fopenmp  
FCFLAGS += -I/usr/local/include
FCFLAGS += -I//usr/local/Cellar/fftw/3.3.8/include
LDFLAGS = -llapack -lblas -lfftw3 -lm 

ifeq ($(USE_LAPACKLIB),LAPACK)
        #FCFLAGS += -I${FFTWINCLUDE}
        #         #FCFLAGS += -I/cm/shared/appsLocal/fftw3_scalar/gcc/3.3.2/include
        #LDFLAGS += -L${FFTWDIR} -lfftw3 -L${LAPACK_DIR} -llapack -L${BLASDIR} -lblas

        RED_SOURCES += algebra_lapack.f90
        ALGEBRA_OBJS += algebra_lapack.o

endif

ifeq ($(USE_LAPACKLIB),MKL)
        FCFLAGS += -I${FFTWINCLUDE} -I$(MKLROOT)/include -I$(MKLROOT)/include/intel64/lp64
        LDFLAGS += -L${FFTWDIR} -lfftw3 -L$(MKLROOT)/lib/intel64 -lmkl_intel_lp64 \
                                        -lmkl_lapack95_lp64 \
                                        -lmkl_sequential -lmkl_core -lpthread -lm \
                                        -Wl,-rpath,$(MKLROOT)/lib/intel64
        RED_SOURCES += algebra_mkl.f90
        ALGEBRA_OBJS += algebra_mkl.o
endif

# Selects the FFT way according to given option
ifeq ($(USE_FFT),GG)                               # GG -> Transforms and Inverse transforms as explained by Gary Glatzmaier
	RED_SOURCES += fourier_GG.f90
        RED_SOURCES += chebyshev_GG.f90
	FFT_OBJS_F += fourier_GG.o
        FFT_OBJS_C += chebyshev_GG.o
else ifeq ($(USE_FFT),FFTW)                        # FFTW -> Full usage of FFTW routines for both Chebyshev and Fourier directions
	RED_SOURCES += fourier_fftw.f90
        RED_SOURCES += chebyshev_fftw.f90
	FFT_OBJS_F += fourier_fftw.o
        FFT_OBJS_C += chebyshev_fftw.o
else ifeq ($(USE_FFT),FFTWF)                       # FFTWF -> FFTW routines used in Fourier direction only and Chebyshev is as in GAG
	RED_SOURCES += fourier_fftw.f90
        RED_SOURCES += chebyshev_GG.f90
	FFT_OBJS_F += fourier_fftw.o
        FFT_OBJS_C += chebyshev_GG.o
else ifeq ($(USE_FFT),FFTWC)                       # FFTWC -> FFTW routines used in Chebyshev direction only and Fourier is as in GAG 
	RED_SOURCES += fourier_GG.f90
        RED_SOURCES += chebyshev_fftw.f90
	FFT_OBJS_F += fourier_GG.o
        FFT_OBJS_C += chebyshev_fftw.o
endif

# List of executables to be built within the package
PROGRAMS = main

DEPS := $(addsuffix .o, $(basename $(RED_SOURCES)))

# "make" builds all
all:  $(PROGRAMS) 

namelists.o: dble.o

constants.o: dble.o

$(ALGEBRA_OBJS): dble.o 

$(FFT_OBJS_C): dble.o constants.o 

$(FFT_OBJS_F): dble.o constants.o 

timeschemes.o: dble.o constants.o 

init.o: dble.o constants.o $(FFT_OBJS_C) $(FFT_OBJS_F) timeschemes.o 

mat_assembly.o: dble.o constants.o $(FFT_OBJS_C) init.o $(ALGEBRA_OBJS) timeschemes.o

output.o: dble.o integ.o init.o $(FFT_OBJS_C) $(FFT_OBJS_F) timeschemes.o 

integ.o: dble.o $(FFT_OBJS_C) init.o

nonlin.o: dble.o constants.o $(FFT_OBJS_F)

rhs_create.o: dble.o constants.o $(FFT_OBJS_C) $(FFT_OBJS_F) timeschemes.o init.o 

fourierloop.o: dble.o constants.o $(FFT_OBJS_C) timeschemes.o init.o mat_assembly.o $(ALGEBRA_OBJS) rhs_create.o

steptime.o: dble.o constants.o $(FFT_OBJS_C) $(FFT_OBJS_F) timeschemes.o init.o mat_assembly.o nonlin.o rhs_create.o fourierloop.o

rhs_create_rk.o: dble.o constants.o $(FFT_OBJS_C) $(FFT_OBJS_F) timeschemes.o init.o

fourierloop_rk.o: dble.o constants.o $(FFT_OBJS_C) timeschemes.o init.o mat_assembly.o $(ALGEBRA_OBJS) rhs_create_rk.o

steptime_rk.o: dble.o constants.o $(FFT_OBJS_C) $(FFT_OBJS_F) timeschemes.o init.o mat_assembly.o nonlin.o fourierloop_rk.o 

rhs_create_imexrk.o: dble.o constants.o $(FFT_OBJS_C) timeschemes.o init.o

fourierloop_imexrk.o: dble.o constants.o $(FFT_OBJS_C) timeschemes.o init.o mat_assembly.o $(ALGEBRA_OBJS) rhs_create_imexrk.o

steptime_imexrk.o: dble.o constants.o $(FFT_OBJS_C) $(FFT_OBJS_F) timeschemes.o init.o mat_assembly.o nonlin.o fourierloop_imexrk.o  

solver.o: dble.o constants.o namelists.o $(FFT_OBJS_C) $(FFT_OBJS_F) timeschemes.o init.o mat_assembly.o integ.o output.o nonlin.o rhs_create.o fourierloop.o steptime.o rhs_create_rk.o fourierloop_rk.o steptime_rk.o rhs_create_imexrk.o fourierloop_imexrk.o steptime_imexrk.o $(ALGEBRA_OBJS)

main.o: dble.o constants.o namelists.o $(FFT_OBJS_C) $(FFT_OBJS_F) timeschemes.o init.o mat_assembly.o integ.o output.o nonlin.o rhs_create.o fourierloop.o steptime.o rhs_create_rk.o fourierloop_rk.o steptime_rk.o rhs_create_imexrk.o fourierloop_imexrk.o steptime_imexrk.o $(ALGEBRA_OBJS) solver.o

main: dble.o constants.o namelists.o $(FFT_OBJS_C) $(FFT_OBJS_F) timeschemes.o init.o mat_assembly.o integ.o output.o nonlin.o rhs_create.o fourierloop.o steptime.o rhs_create_rk.o fourierloop_rk.o steptime_rk.o rhs_create_imexrk.o fourierloop_imexrk.o steptime_imexrk.o $(ALGEBRA_OBJS) solver.o

%: %.o
	$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.F90 
	$(FC) $(FCFLAGS) -c $< 

%.o: %.F90 $(DEPS)
	$(FC) -c -o $@ $<

.PHONY: clean veryclean

clean:
	rm -f *.o *.mod *.MOD

veryclean: clean
	rm -f *~ $(PROGRAMS)

