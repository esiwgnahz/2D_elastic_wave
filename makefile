all: main

FLAGS  = -g -O0 -I${PETSC_DIR}/include -I$(PETSC_DIR)/$(PETSC_ARCH)/include

OBJECT = output.o shape_functions.o matrix_elemental.o matrix_assemble.o body_force.o model.o solver.o main.o

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

main: $(OBJECT) chkopts
	-$(CLINKER) $(FLAGS) -o main $(OBJECT) $(PETSC_KSP_LIB)

# %.o: %.c
# 	gcc $(FLAGS) $< -c

include ${PETSC_DIR}/lib/petsc/conf/test
