EXECS=prog1
MPICC?=mpicxx

all: ${EXECS}

prog1: prog1.cpp
	${MPICC} -std=c++14 -o prog1 prog1.cpp

clean:
	rm -f ${EXECS}
