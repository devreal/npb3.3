SHELL=/bin/sh
BENCHMARK=bt-mz
BENCHMARKU=BT-MZ
VEC=
IMPL=cxx-omptasks

include ../config/make.def


OBJS = bt.o exact_rhs.o exch_qbc.o  mpi_setup.o  set_constants.o  zone_setup.o \
        ${COMMON}/c_print_results.o ${COMMON}/c_timers.o

CXX=mpicxx

include ../sys/make.common

CXXLINKFLAGS += -fopenmp
CXXFLAGS += -fopenmp


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

exec: $(OBJS)
	${CXXLINK} ${CXXLINKFLAGS} -o ${PROGRAM} ${OBJS} ${F_LIB}

.f.o:
	${FCOMPILE} $<

bt.o:             bt.cc  mpi_stuff.h header.h npbparams.h
initialize.o:     initialize.cc  header.h npbparams.h
set_constants.o:  set_constants.cc  header.h npbparams.h
zone_setup.o:     zone_setup.cc mpi_stuff.h header.h npbparams.h
verify.o:         verify.cc  mpi_stuff.h header.h npbparams.h
exch_qbc.o:       exch_qbc.cc mpi_stuff.h header.h npbparams.h
mpi_setup.o:      mpi_setup.cc  mpi_stuff.h header.h npbparams.h

clean:
	- rm -f *.o *~ mputil*
	- rm -f npbparams*.h core
