#!/bin/bash
#SBATCH --job-name=md_comparison
#SBATCH --time=6:00:00
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=serial
#SBATCH --qos=serial
#SBATCH --account=

# -------------------------------------------------------
# Full comparison benchmark:
#   1. Baseline (original, -O0 -g)
#   2. gprof profile of baseline
#   3. Optimised code with -O0 -g  (isolate source changes)
#   4. Optimised code with production flags
#   5. Original code with production flags (isolate flag effect)
#   6. Correctness verification for all variants
# -------------------------------------------------------

module load PrgEnv-gnu

CW2=${SLURM_SUBMIT_DIR}
ORIG_DIR=${CW2}/original
OPT_DIR=${CW2}/optimized
RESULTS=${CW2}/results
NRUNS=5

cd ${CW2}

# Force-create results directory
[ -f "${RESULTS}" ] && rm -f "${RESULTS}"
mkdir -p ${RESULTS}

echo "=== MD Full Comparison Benchmark ==="
echo "Date: $(date)"
echo "Node: $(hostname)"
echo "Compiler: $(cc --version | head -1)"
echo "Runs per variant: $NRUNS"
echo "CW2 root: ${CW2}"
echo ""

# Build correctness checker
cc -O2 -o ${RESULTS}/diff-output ${CW2}/Test/diff-output.c -lm

# Use a temporary build directory to avoid polluting source dirs
BUILDDIR=$(mktemp -d ${CW2}/build_XXXXXX)
cd ${BUILDDIR}
cp ${CW2}/input.dat .

# -------------------------------------------------------
# Helper function
# -------------------------------------------------------
run_variant() {
    local name="$1"
    local src_dir="$2"
    local cflags="$3"
    local ldflags="$4"

    echo "========================================"
    echo "Variant: $name"
    echo "Source:  $src_dir"
    echo "CFLAGS:  $cflags"
    echo "========================================"

    rm -f MD *.o output.dat*
    cp ${src_dir}/MD.c ${src_dir}/control.c ${src_dir}/util.c ${src_dir}/coord.h .

    cc $cflags -c MD.c control.c util.c 2> ${RESULTS}/compile_${name}.err
    cc $cflags -o MD MD.o control.o util.o $ldflags 2>> ${RESULTS}/compile_${name}.err
    if [ $? -ne 0 ]; then
        echo "  COMPILATION FAILED -- see results/compile_${name}.err"
        cat ${RESULTS}/compile_${name}.err
        echo ""
        return 1
    fi

    for run in $(seq 1 $NRUNS); do
        rm -f output.dat*
        ./MD > ${RESULTS}/run_${name}_${run}.log 2>&1
        echo -n "  Run $run: "
        grep "timesteps took" ${RESULTS}/run_${name}_${run}.log | tail -1
    done

    cp output.dat00500 ${RESULTS}/out_${name}.dat 2>/dev/null
    echo ""
}

# -------------------------------------------------------
# 1. BASELINE: original code, -O0 -g
# -------------------------------------------------------
run_variant "baseline" "$ORIG_DIR" "-O0 -g" "-lm"

# -------------------------------------------------------
# 2. gprof profile of baseline
# -------------------------------------------------------
echo "=== gprof profile of baseline ==="
rm -f MD *.o output.dat* gmon.out
cp ${ORIG_DIR}/MD.c ${ORIG_DIR}/control.c ${ORIG_DIR}/util.c ${ORIG_DIR}/coord.h .
cc -O0 -g -pg -c MD.c control.c util.c
cc -O0 -g -pg -o MD MD.o control.o util.o -lm
./MD > /dev/null 2>&1
gprof MD gmon.out > ${RESULTS}/gprof_baseline.txt 2>&1
echo "  Saved: results/gprof_baseline.txt"
head -30 ${RESULTS}/gprof_baseline.txt
echo ""

# -------------------------------------------------------
# 3. OPTIMISED CODE: -O0 -g (isolate source-level changes)
# -------------------------------------------------------
run_variant "opt_O0" "$OPT_DIR" "-O0 -g" "-lm"

# -------------------------------------------------------
# 4. OPTIMISED CODE: production flags
# -------------------------------------------------------
run_variant "opt_prod" "$OPT_DIR" \
    "-O3 -march=znver2 -ffast-math -funroll-loops -flto -Wall" \
    "-flto -lm"

# -------------------------------------------------------
# 5. ORIGINAL CODE: production flags (isolate flag effect)
# -------------------------------------------------------
run_variant "orig_prod" "$ORIG_DIR" \
    "-O3 -march=znver2 -ffast-math -funroll-loops -flto -Wall" \
    "-flto -lm"

# -------------------------------------------------------
# 6. Correctness verification
# -------------------------------------------------------
echo "========================================"
echo "Correctness Verification"
echo "========================================"
CHECKER=${RESULTS}/diff-output
for v in opt_O0 opt_prod orig_prod; do
    if [ -f ${RESULTS}/out_${v}.dat ] && [ -f ${RESULTS}/out_baseline.dat ]; then
        echo -n "  $v vs baseline: "
        ${CHECKER} ${RESULTS}/out_${v}.dat ${RESULTS}/out_baseline.dat 2>&1 | tail -1
    fi
done

# Clean up temp build directory
rm -rf ${BUILDDIR}

echo ""
echo "=== All benchmarks complete ==="
echo "All results in: ${RESULTS}/"
echo ""
echo "Key files to copy back:"
echo "  results/gprof_baseline.txt   - for report Table 2"
echo "  results/run_*_*.log          - timing data for report tables"
echo "  results/out_*.dat            - correctness verification"
echo "  slurm-*.out                  - full job output"
