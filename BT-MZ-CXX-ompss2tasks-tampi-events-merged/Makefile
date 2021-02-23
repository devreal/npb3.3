SHELL=/bin/sh
BENCHMARK=bt-mz
BENCHMARKU=BT-MZ
VEC=
IMPL=cxx-ompss2tasks-tampi-events-merged

include ../config/make.def

HEADER := $(wildcard *.h)

OBJS = bt.o  mpi_setup.o  set_constants.o  zone_setup.o \
        ${COMMON}/c_print_results.o ${COMMON}/c_timers.o add_impl.o  bt_omp.o \
        exact_rhs_impl.o exact_rhs_omp.o exch_qbc_impl.o exch_qbc_omp.o x_solve_impl.o \
        y_solve_impl.o z_solve_impl.o error_omp.o error_impl.o initialize_impl.o initialize_omp.o \
        rhs_impl.o rhs_omp.o print_zone.o


include ../sys/make.common

# Flags used for your MPI implementation
MPICXXFLAGS=$(shell mpicxx -showme:compile)
#MPICXXFLAGS=$(shell CC --cray-print-opts=cflags)
MPILINKFLAGS=$(shell mpicxx -showme:link)
#MPILINKFLAGS=$(shell CC --cray-print-opts=libs)

OMPSSROOT=${HOME}/opt/ompss-2-2019.11.2/include/

# The OmpSs-2 compiler requires a different set of arguments, e.g., C++11
OMPCC = mcxx --cxx=clang++ -std=c++11 --ompss-2 -k -g -O2 -march=haswell -DNDEBUG -mcmodel=medium

TAMPI_INC=-I$(TAMPIROOT)/include
TAMPI_LIB=-L$(TAMPIROOT)/lib -ltampi

CXXFLAGS+=$(TAMPI_INC)

# npbparams.h is included by header.h
# The following rule should do the trick but many make programs (not gmake)
# will do the wrong thing and rebuild the world every time (because the
# mod time on header.h is not changed. One solution would be to
# touch header.h but this might cause confusion if someone has
# accidentally deleted it. Instead, make the dependency on npbparams.h
# explicit in all the lines below (even though dependence is indirect).

# header.h: npbparams.h

${PROGRAM}: config
	@if [ x$(VERSION) = xvec ] ; then	\
		${MAKE} VEC=_vec exec;		\
	elif [ x$(VERSION) = xVEC ] ; then	\
		${MAKE} VEC=_vec exec;		\
	else					\
		${MAKE} exec;			\
	fi

%_omp.o: %_omp.cc npbparams.h
	${OMPCC} $(MPICXXFLAGS) $(TAMPI_INC) -c $< -o $@

exec: $(OBJS)
	${OMPCC} ${CXXLINKFLAGS} ${TAMPI_LIB} ${MPILINKFLAGS} -o ${PROGRAM} ${OBJS} ${F_LIB}

.f.o:
	${FCOMPILE} $<

bt.o:             bt.cc $(HEADER)
initialize.o:     initialize.cc $(HEADER)
set_constants.o:  set_constants.cc $(HEADER)
zone_setup.o:     zone_setup.cc $(HEADER)
verify.o:         verify.cc $(HEADER)
exch_qbc.o:       exch_qbc.cc $(HEADER)
mpi_setup.o:      mpi_setup.cc $(HEADER)
add_impl.o:       add_impl.cc $(HEADER)
error_impl.o:     error_impl.cc $(HEADER)
exact_rhs_impl.o: exact_rhs_impl.cc $(HEADER)
exch_qbc_impl.o:  exch_qbc_impl.cc $(HEADER)
x_solve_impl.o:   x_solve_impl.cc $(HEADER)
y_solve_impl.o:   y_solve_impl.cc $(HEADER)
z_solve_impl.o:   z_solve_impl.cc $(HEADER)

clean:
	- rm -f *.o *~ mputil*
	- rm -f npbparams*.h core
