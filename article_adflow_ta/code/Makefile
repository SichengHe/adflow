PETSC_DIR = /opt/homebrew/opt/petsc

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

ex20adj: ex20adj.o
	-${PCC_LINKER} -o $@ $< ${PETSC_TS_LIB}
	${RM} ex20adj.o

.PHONY: clean
clean::
	${RM} ex20adj *.o
