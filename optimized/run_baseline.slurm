#!/bin/bash
#SBATCH --job-name=md_baseline
#SBATCH --time=4:00:00
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=serial
#SBATCH --qos=serial
#SBATCH --account=

# -------------------------------------------------------
# Baseline benchmark: original code compiled with -O0 -g
# -------------------------------------------------------

module load PrgEnv-gnu
CW2=${SLURM_SUBMIT_DIR}
RESULTS=${CW2}/results
cd ${CW2}

# Force-create results directory
[ -f "${RESULTS}" ] && rm -f "${RESULTS}"
mkdir -p ${RESULTS}

echo "=== Baseline MD Benchmark ==="
echo "Date: $(date)"
echo "Node: $(hostname)"
echo "Compiler: $(cc --version | head -1)"
echo "Working dir: $(pwd)"
echo ""

# Use a temporary build directory
BUILDDIR=$(mktemp -d ${CW2}/build_baseline_XXXXXX)
cd ${BUILDDIR}

# Copy original source + input data
cp ${CW2}/original/MD.c ${CW2}/original/control.c \
   ${CW2}/original/util.c ${CW2}/original/coord.h \
   ${CW2}/original/Makefile .
cp ${CW2}/input.dat .

# Compile with -O0 -g (no optimisation)
make clean 2>/dev/null
make CFLAGS="-O0 -g" LDFLAGS="-lm"
echo "Compiled with: -O0 -g"
echo ""

# Run 5 times
for run in 1 2 3 4 5; do
    echo "--- Run $run ---"
    rm -f output.dat*
    ./MD > ${RESULTS}/run_baseline_${run}.log 2>&1
    grep "timesteps took" ${RESULTS}/run_baseline_${run}.log | tail -1
    if [ $run -eq 1 ]; then
        cp output.dat00500 ${RESULTS}/output_baseline.dat
    fi
done

# Clean up temp build directory
rm -rf ${BUILDDIR}

echo ""
echo "=== Baseline complete ==="
echo "Reference output: results/output_baseline.dat"
echo "Timing logs:      results/run_baseline_*.log"
