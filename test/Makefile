all: 

INC_DIR=../include

LIB_DIR=../lib

LIBST=-lmnblas
LIBDYN=-lmnblasdyn

OPTIONS_COMPIL  =-Wall -O2 -fPIC -I$(INC_DIR)

OPTIONS_LINK_STATIC  =   -fopenmp -L$(LIB_DIR) $(LIBST)
OPTIONS_LINK_DYN  = -fopenmp -L$(LIB_DIR) $(LIBDYN)


#test_dot: test_dot.o flop.o $(LIB_DIR)/libmnblas.a
#	gcc -o test_dot test_dot.o flop.o $(OPTIONS_LINK_STATIC)

#test_dot.o: test_dot.c $(INC_DIR)/mnblas.h
#	gcc $(OPTIONS_COMPIL) -c test_dot.c


clean:
	rm -f *.o *~
