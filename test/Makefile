PROG=test_copy test_asum test_iamax test_iamin test_nrm2 test_swap test_axpy test_gemv
all: $(PROG)

INC_DIR=../include

LIB_DIR=../lib

LIBST=-lmnblas
LIBDYN=-lmnblasdyn

OPTIONS_COMPIL  =-Wall -O2 -fPIC -I$(INC_DIR)

OPTIONS_LINK_STATIC  =   -fopenmp -L$(LIB_DIR) $(LIBST)
OPTIONS_LINK_DYN  = -fopenmp -L$(LIB_DIR) $(LIBDYN)


test_copy: test_copy.o flop.o $(LIB_DIR)/libmnblas.a
	gcc -o test_copy test_copy.o flop.o $(OPTIONS_LINK_STATIC)

test_copy.o: test_copy.c $(INC_DIR)/mnblas.h
	gcc $(OPTIONS_COMPIL) -c test_copy.c

test_asum: test_asum.o flop.o $(LIB_DIR)/libmnblas.a
	gcc -o test_asum test_asum.o flop.o $(OPTIONS_LINK_STATIC)

test_asum.o: test_asum.c $(INC_DIR)/mnblas.h
	gcc $(OPTIONS_COMPIL) -c test_asum.c

test_iamax: test_iamax.o flop.o $(LIB_DIR)/libmnblas.a
	gcc -o test_iamax test_iamax.o flop.o $(OPTIONS_LINK_STATIC)

test_iamax.o: test_iamax.c $(INC_DIR)/mnblas.h
	gcc $(OPTIONS_COMPIL) -c test_iamax.c

test_iamin: test_iamin.o flop.o $(LIB_DIR)/libmnblas.a
	gcc -o test_iamin test_iamin.o flop.o $(OPTIONS_LINK_STATIC)

test_iamin.o: test_iamin.c $(INC_DIR)/mnblas.h
	gcc $(OPTIONS_COMPIL) -c test_iamin.c

test_nrm2: test_nrm2.o flop.o $(LIB_DIR)/libmnblas.a
	gcc -o test_nrm2 test_nrm2.o flop.o $(OPTIONS_LINK_STATIC) -lm

test_nrm2.o: test_nrm2.c $(INC_DIR)/mnblas.h
	gcc -lm $(OPTIONS_COMPIL) -c test_nrm2.c

test_axpy: test_axpy.o flop.o $(LIB_DIR)/libmnblas.a
	gcc -o test_axpy test_axpy.o flop.o $(OPTIONS_LINK_STATIC)

test_axpy.o: test_axpy.c $(INC_DIR)/mnblas.h
	gcc -lm $(OPTIONS_COMPIL) -c test_axpy.c

test_swap: test_swap.o flop.o $(LIB_DIR)/libmnblas.a
	gcc -o test_swap test_swap.o flop.o $(OPTIONS_LINK_STATIC)

test_swap.o: test_swap.c $(INC_DIR)/mnblas.h
	gcc $(OPTIONS_COMPIL) -c test_swap.c

test_gemv: test_gemv.o flop.o $(LIB_DIR)/libmnblas.a
	gcc -o test_gemv test_gemv.o flop.o $(OPTIONS_LINK_STATIC)

test_gemv.o: test_gemv.c $(INC_DIR)/mnblas.h
	gcc $(OPTIONS_COMPIL) -c test_gemv.c

test: $(PROG)
	@echo;echo "*** TESTS ***";echo;\
	for file in $^;\
		do echo ;\
		./$$file 2>/dev/null || break;\
	done
clean:
	rm -f *.o $(PROG) *~
