#!/bin/bash
# Check for sensitive information in portfolio files

echo "Checking for potential sensitive information..."

check_file() {
    local file=$1
    if [ -f "$file" ]; then
        if grep -q "m25oc-s" "$file" 2>/dev/null; then
            echo "WARNING: Found account ID in $file"
            grep -n "m25oc-s" "$file"
        fi
    fi
}

check_file "run_gnu_experiments.slurm"
check_file "run_aocc_experiments.slurm"
check_file "run_cray_experiments.slurm"
check_file "bench_c.slurm"

echo "Check complete."
