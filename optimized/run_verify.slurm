#!/bin/bash
#SBATCH --job-name=md_verify
#SBATCH --time=00:10:00
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=serial
#SBATCH --qos=serial
#SBATCH --account=

# Fast correctness verification (~5 min total):
# 1. Run original code + production flags  -> reference output (~68s, known correct)
# 2. Run optimised code + production flags -> test output    (~27s)
# 3. Run optimised code at -O0            -> check collision count only (~87s)
# 4. diff-output comparison of 1 vs 2
# No -O0 baseline re-run needed (would take ~875s).

module load PrgEnv-gnu
CW2=${SLURM_SUBMIT_DIR}
RESULTS=${CW2}/results
[ -f "${RESULTS}" ] && rm -f "${RESULTS}"
mkdir -p ${RESULTS}

# Explicit flag sets - override whatever the Makefile defaults are
PROD_CFLAGS="-O3 -march=znver2 -ffast-math -funroll-loops -flto"
PROD_LDFLAGS="-flto -lm"
O0_CFLAGS="-O0 -g"
O0_LDFLAGS="-lm"

echo "=== Correctness Verification (fast) ==="
echo "Date: $(date)"
echo "Node: $(hostname)"
echo ""

# Build diff-output tool
cc -O2 -o ${RESULTS}/diff-output ${CW2}/Test/diff-output.c -lm

# -------------------------------------------------------
# Step 1: Original code + production flags (reference)
# Override Makefile CFLAGS: original Makefile has -fsanitize=address
# -------------------------------------------------------
REF=$(mktemp -d ${CW2}/verify_ref_XXXXXX)
cp ${CW2}/original/MD.c ${CW2}/original/control.c \
   ${CW2}/original/util.c ${CW2}/original/coord.h \
   ${CW2}/original/Makefile ${CW2}/input.dat ${REF}/
cd ${REF}
make clean 2>/dev/null
make CFLAGS="${PROD_CFLAGS}" LDFLAGS="${PROD_LDFLAGS}" 2>&1
echo "=== Original + prod flags ==="
./MD > ${RESULTS}/verify_orig_prod.log 2>&1
REF_COLL=$(grep "^collisions" ${RESULTS}/verify_orig_prod.log | tail -1)
REF_TIME=$(grep "timesteps took" ${RESULTS}/verify_orig_prod.log | tail -1)
echo "  Collisions: ${REF_COLL}"
echo "  Time:       ${REF_TIME}"
for f in output.dat500 output.dat00500; do
    [ -f "$f" ] && cp "$f" ${RESULTS}/verify_ref.dat && echo "  Output saved: $f"
done
cd ${CW2} && rm -rf ${REF}

# -------------------------------------------------------
# Step 2: Optimised code + production flags (already correct Makefile,
# but pass flags explicitly for safety)
# -------------------------------------------------------
OPT=$(mktemp -d ${CW2}/verify_opt_XXXXXX)
cp ${CW2}/optimized/MD.c ${CW2}/optimized/control.c \
   ${CW2}/optimized/util.c ${CW2}/optimized/coord.h \
   ${CW2}/optimized/Makefile ${CW2}/input.dat ${OPT}/
cd ${OPT}
make clean 2>/dev/null
make CFLAGS="${PROD_CFLAGS}" LDFLAGS="${PROD_LDFLAGS}" 2>&1
echo ""
echo "=== Optimised + prod flags ==="
./MD > ${RESULTS}/verify_opt_prod.log 2>&1
OPT_COLL=$(grep "^collisions" ${RESULTS}/verify_opt_prod.log | tail -1)
OPT_TIME=$(grep "timesteps took" ${RESULTS}/verify_opt_prod.log | tail -1)
echo "  Collisions: ${OPT_COLL}"
echo "  Time:       ${OPT_TIME}"
for f in output.dat500 output.dat00500; do
    [ -f "$f" ] && cp "$f" ${RESULTS}/verify_opt_prod.dat && echo "  Output saved: $f"
done
cd ${CW2} && rm -rf ${OPT}

# -------------------------------------------------------
# Step 3: Optimised code at -O0 (collision count check, ~87s)
# -------------------------------------------------------
OPT0=$(mktemp -d ${CW2}/verify_opt0_XXXXXX)
cp ${CW2}/optimized/MD.c ${CW2}/optimized/control.c \
   ${CW2}/optimized/util.c ${CW2}/optimized/coord.h \
   ${CW2}/optimized/Makefile ${CW2}/input.dat ${OPT0}/
cd ${OPT0}
make clean 2>/dev/null
make CFLAGS="${O0_CFLAGS}" LDFLAGS="${O0_LDFLAGS}" 2>&1
echo ""
echo "=== Optimised -O0 (collision check only) ==="
./MD > ${RESULTS}/verify_opt_O0.log 2>&1
OPT0_COLL=$(grep "^collisions" ${RESULTS}/verify_opt_O0.log | tail -1)
OPT0_TIME=$(grep "timesteps took" ${RESULTS}/verify_opt_O0.log | tail -1)
echo "  Collisions: ${OPT0_COLL}"
echo "  Time:       ${OPT0_TIME}"
cd ${CW2} && rm -rf ${OPT0}

# -------------------------------------------------------
# Summary and diff-output comparison
# -------------------------------------------------------
echo ""
echo "========================================="
echo "COLLISION COUNT SUMMARY (expected: 68827)"
echo "  Original + prod:  ${REF_COLL}"
echo "  Optimised + prod: ${OPT_COLL}"
echo "  Optimised -O0:    ${OPT0_COLL}"
echo "========================================="

if [ -f ${RESULTS}/verify_ref.dat ] && [ -f ${RESULTS}/verify_opt_prod.dat ]; then
    echo ""
    echo "=== diff-output: optimised+prod vs original+prod ==="
    ${RESULTS}/diff-output ${RESULTS}/verify_opt_prod.dat ${RESULTS}/verify_ref.dat
else
    echo "WARNING: output .dat files not found for diff-output comparison"
    ls -la ${RESULTS}/verify_*.dat 2>/dev/null || echo "No verify_*.dat files in results/"
fi

echo ""
echo "=== Verification complete: $(date) ==="
