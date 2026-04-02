SRC = MD.c control.c util.c 
OBJ = $(SRC:.c=.o)
CC = cc

# OPTIMIZED FLAGS for ARCHER2 - GNU Compiler
# -O3: Maximum optimization level
# -march=znver2: Optimize for AMD EPYC 7742 (ARCHER2 compute nodes)
# -ffast-math: Allow aggressive floating-point optimizations
# -funroll-loops: Unroll loops for better instruction pipelining
# -ftree-vectorize: Enable auto-vectorization (default with -O3)
CFLAGS = -O3 -march=znver2 -ffast-math -funroll-loops -flto -Wall

# Linker flags: math library + link-time optimization
# -flto enables cross-module inlining and optimisation at link time
LDFLAGS = -flto -lm

all: MD

MD: $(OBJ)
	$(CC) $(CFLAGS) -o $@ $(OBJ) $(LDFLAGS)

clean:
	rm -f MD $(OBJ) output.dat*

# Dependency tracking
$(OBJ): coord.h Makefile

.PHONY: all clean
