#!/bin/bash
echo     "-----------------------------------------------------------"
echo "RUNNING TEST CASE 1: 2D-ADVECTION"
echo     "-----------------------------------------------------------"
../../bin/Ash3d TC1_XY_MSH.inp > /dev/null 2>&1
echo "   Checking cloud height at 2 hours:"
../../bin/tools/Ash3d_ASCII_check CloudHeight_002.00hrs.dat output/CloudHeight_002.00hrs.dat 1.0e-1
echo "   Checking cloud height at 4 hours:"
../../bin/tools/Ash3d_ASCII_check CloudHeight_004.00hrs.dat output/CloudHeight_004.00hrs.dat 1.0e-1
echo "   Checking cloud load at 2 hours:"
../../bin/tools/Ash3d_ASCII_check CloudLoad_002.00hrs.dat output/CloudLoad_002.00hrs.dat 1.0e-1
echo "   Checking cloud load at 4 hours:"
../../bin/tools/Ash3d_ASCII_check CloudLoad_004.00hrs.dat output/CloudLoad_004.00hrs.dat 1.0e-1

