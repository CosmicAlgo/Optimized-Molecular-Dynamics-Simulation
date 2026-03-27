#!/bin/bash
# Verify redaction of sensitive data

echo "Verifying redaction..."

if grep -r "m25oc-s" . --include="*.sh" --include="*.slurm" 2>/dev/null; then
    echo "ERROR: Found unredacted account IDs"
    exit 1
else
    echo "OK: No account IDs found in scripts"
fi

if grep -r "B288585" . --include="*.c" --include="*.h" 2>/dev/null; then
    echo "WARNING: Found exam number in source files"
fi

echo "Verification complete."
