#!/bin/bash

resultsdir="./vtune_results$1/"
rm -r $resultsdir
./fix_ptrace_scope.sh vtune -collect hotspots -r $resultsdir -- ./test/integration_tests --gtest_filter=FaceIter.Performance

# to see results in GUI: vtune-gui ./vtune_results/vtune.vtune
