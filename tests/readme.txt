These tests in directories 1-4 are quick tests that are run via the script
src/run_tests.sh which is launched from 'make check'.  These tests are
designed to run in just a few minutes and to verify that changes to the 
code do not cause a significant change in the output.

All these tests produce ASCII output files which are compared with files
stored in the test_0[1-4]/output folders with the tool, Ash3d_ASCII_check.

Note that these tests will product slightly different output and may produce
'failure' results if a variety of different options are reset, such as:
the CFL condition, Explicit vs Crank-Nicolson, limiters, fast_dt, 
fast_subgrid, useVz_rhoG, etc.

test_01:
This test used a 1-d ASCII windfile and tests fallout of a single grain-
size using a Suzuki source.

test_02:
This test invokes 3d homogeneous diffusion.

test_03:
These 11 subcases test different source types, fall models and vertical
grid specifications.
 Sub-case  0: Suz=4;   WH;          const dz
 Sub-case  1: line;    WH;          const dz
 Sub-case  2: point;   tracer;      const dz
 Sub-case  3: point;   WH;          const dz
 Sub-case  4: point;   Ganser;      const dz
 Sub-case  5: point;   Stokes/slip; const dz
 Sub-case  6: profile; WH;          const dz
 Sub-case  7: line;    WH;          const dz
 Sub-case  8: point;   WH;          piece-wise linear
 Sub-case  9: point;   WH;          constant log
 Sub-case 10: point;   WH;          custom

test_04:
These 5 subcases test more realistic cases with transient, 3d winds (NCEP).
 Sub-case  0: Suz=4
 Sub-case  1: umbrella_air
 Sub-case  2: Suz=4; 12-bin 'deposit' grainsize distribution
 Sub-case  3: umbrella
 Sub-case  4: periodic BC (global)
 Sub-case  5: 1-D wind, Topo with Z_ID=0 (no actual z-grid modification)
 Sub-case  6: 1-D wind, Topo with Z_ID=1 (shifted: effectively same as 5)
 Sub-case  7: 1-D wind, Topo with Z_ID=2 (scaled: same source as 5/6 but with dz_clog)
