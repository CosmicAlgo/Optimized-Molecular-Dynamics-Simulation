#!/bin/bash
#SBATCH --job-name=md_optimized
#SBATCH --time=2:00:00
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=serial
#SBATCH --qos=serial
#SBATCH --account=

# -------------------------------------------------------
# Optimized benchmark: optimized code with production flags
# -------------------------------------------------------

module load PrgEnv-gnu
CW2=${SLURM_SUBMIT_DIR}
RESULTS=${CW2}/results
cd ${CW2}

# Force-create results directory
[ -f "${RESULTS}" ] && rm -f "${RESULTS}"
mkdir -p ${RESULTS}

echo "=== Optimized MD Benchmark ==="
echo "Date: $(date)"
echo "Node: $(hostname)"
echo "Compiler: $(cc --version | head -1)"
echo "Working dir: $(pwd)"
echo ""

# Use a temporary build directory
BUILDDIR=$(mktemp -d ${CW2}/build_opt_XXXXXX)
cd ${BUILDDIR}

# Copy optimized source + input data
cp ${CW2}/optimized/MD.c ${CW2}/optimized/control.c \
   ${CW2}/optimized/util.c ${CW2}/optimized/coord.h \
   ${CW2}/optimized/Makefile .
cp ${CW2}/input.dat .

# Compile with production flags
make clean 2>/dev/null
make
echo "CFLAGS: $(grep ^CFLAGS Makefile)"
echo ""

# Build correctness checker
cc -O2 -o diff-output ${CW2}/Test/diff-output.c -lm

# Run 5 times
for run in 1 2 3 4 5; do
    echo "--- Run $run ---"
    rm -f output.dat*
    ./MD > ${RESULTS}/run_opt_${run}.log 2>&1
    grep "timesteps took" ${RESULTS}/run_opt_${run}.log | tail -1
    if [ $run -eq 5 ]; then
        cp output.dat00500 ${RESULTS}/output_optimized.dat
    fi
done

# Clean up temp build directory
rm -rf ${BUILDDIR}

# Correctness check against baseline (if available)
if [ -f ${RESULTS}/output_baseline.dat ]; then
    echo ""
    echo "=== Correctness Check ==="
    ${RESULTS}/diff-output ${RESULTS}/output_optimized.dat ${RESULTS}/output_baseline.dat
fi

echo ""
echo "=== Done ==="
echo "Timing logs: results/run_opt_*.log"
echo "Output:      results/output_optimized.dat"
