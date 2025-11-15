#!/bin/bash
# Script to run CFD comparison tests

echo "=========================================="
echo "ADFLOW TORUS VS 1D TIME SPECTRAL TESTS"
echo "=========================================="
echo ""

cd /home/sicheng/repo/adflow/adflow/tests

echo "1. UNIT TEST (fast - just matrix comparison)"
echo "   Compares spectral operators for n2=1 case"
echo "   Command: testflo unit_tests/test_torus_degenerate_matches_ts.py -v"
echo ""
read -p "Run this test? (y/n) " -n 1 -r
echo
if [[ $REPLY =~ ^[Yy]$ ]]; then
    testflo unit_tests/test_torus_degenerate_matches_ts.py -v
fi
echo ""

echo "2. OPERATOR COMPARISON (fast - no CFD solve)"
echo "   Compares dscalar matrices for torus n2=1 vs 1D TS"
echo "   Command: python compare_torus_n2_1_vs_1d.py"
echo ""
read -p "Run this test? (y/n) " -n 1 -r
echo
if [[ $REPLY =~ ^[Yy]$ ]]; then
    python compare_torus_n2_1_vs_1d.py
fi
echo ""

echo "3. FULL CFD TEST (slow ~1 min - solves both cases)"
echo "   Runs 3×3 torus and 3-pt 1D TS, compares CL/CD"
echo "   Command: python compare_degenerate_torus_vs_1d_ts.py"
echo ""
read -p "Run this test? (y/n) " -n 1 -r
echo
if [[ $REPLY =~ ^[Yy]$ ]]; then
    python compare_degenerate_torus_vs_1d_ts.py
fi
echo ""

echo "4. TRUE DEGENERATE TEST (ω2=0 case)"
echo "   Verifies torus with ω2=0 gives consistent results"
echo "   Command: python test_true_degenerate_torus.py"
echo ""
read -p "Run this test? (y/n) " -n 1 -r
echo
if [[ $REPLY =~ ^[Yy]$ ]]; then
    python test_true_degenerate_torus.py
fi

echo ""
echo "=========================================="
echo "DONE"
echo "=========================================="
