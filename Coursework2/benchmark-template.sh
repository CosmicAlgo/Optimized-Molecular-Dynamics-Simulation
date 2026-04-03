#!/bin/bash
# Sanitized benchmark script for portfolio/GitHub
# (No ARCHER2-specific account ID)
# 
# For use in public repositories and documentation.
# Original scripts with actual account IDs are excluded from git.

#SBATCH --job-name=md_benchmark
#SBATCH --time=2:00:00
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=serial
#SBATCH --qos=serial
# NOTE: --account line removed for portfolio version
# In actual ARCHER2 execution, add: #SBATCH --account=<YOUR_ACCOUNT_ID>

module load PrgEnv-gnu

cd $SLURM_SUBMIT_DIR

# Compile correctness checker
cc -O3 -o Test/diff-output Test/diff-output.c -lm
chmod +x Test/diff-output

# Create results directory
mkdir -p results_benchmark
cd results_benchmark

# Copy source files
cp ../MD.c ../control.c ../util.c ../coord.h ../Makefile ../input.dat .

run_test() {
    local flag_name=$1
    local cflags=$2
    
    echo "running ${flag_name}..."
    
    make clean > /dev/null
    make CFLAGS="${cflags}" > /dev/null
    
    if [ $? -ne 0 ]; then
        echo "Compilation failed for ${flag_name}"
        return 1
    fi
    
    ./MD > "log_${flag_name}.txt"
    
    # Save output for comparison
    cp output.dat00500 "out_${flag_name}.dat"
    
    # Extract timing
    grep "timesteps took" "log_${flag_name}.txt"
    
    # Compare against baseline for correctness
    if [ "${flag_name}" != "baseline" ]; then
        ../Test/diff-output "out_${flag_name}.dat" out_baseline.dat > "diff_${flag_name}.txt"
        if [ $? -eq 0 ]; then
            echo "Correctness check passed for ${flag_name}"
        else
            echo "Correctness check FAILED for ${flag_name}"
        fi
    fi
}

# Run baseline first
run_test "baseline" "-O0 -g"

# Test optimization levels
run_test "O1" "-O1"
run_test "O2" "-O2"
run_test "O3" "-O3"
run_test "Ofast" "-Ofast"

# Test advanced optimization combinations
run_test "O3_native" "-O3 -march=native"
run_test "Ofast_native_LTO" "-Ofast -march=native -flto -funroll-loops"

echo "Benchmark complete. Results in results_benchmark/"
