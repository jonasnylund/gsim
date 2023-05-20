#!/bin/bash

name=$1
mkdir -p $name
pushd $name;

../benchmark.sh "exact" "--updateratio 0.0 --theta 0.2 --dtime 0.001" 1
../benchmark.sh "baseline" "--updateratio 0.0 --theta 0.4 --dtime 0.01" 1 #3

# ../benchmark.sh "fast" "--updateratio 0.414 --theta 0.8 --dtime 0.05" 1 #3

# ../benchmark.sh "reduced_phi" "--updateratio 0.414 --theta 0.4 --dtime 0.01" 1 #3
# ../benchmark.sh "reduced_theta" "--updateratio 0.0 --theta 0.8 --dtime 0.01" 1 #3
# ../benchmark.sh "reduced_dt" "--updateratio 0.0 --theta 0.4 --dtime 0.02" 1 #3

# ../benchmark.sh "balanced1" "--updateratio 0.25 --theta 0.4 --dtime 0.0075" 1 #3
# ../benchmark.sh "balanced2" "--updateratio 0.414 --theta 0.4 --dtime 0.007" 1 #3

# ../benchmark.sh "phi_00" "--theta 0.4 --dtime 0.05 --updateratio 0.0" 3
# ../benchmark.sh "phi_01" "--theta 0.4 --dtime 0.05 --updateratio 0.1" 3
# ../benchmark.sh "phi_02" "--theta 0.4 --dtime 0.05 --updateratio 0.2" 3
# ../benchmark.sh "phi_03" "--theta 0.4 --dtime 0.05 --updateratio 0.3" 3
# ../benchmark.sh "phi_04" "--theta 0.4 --dtime 0.05 --updateratio 0.4" 3
# ../benchmark.sh "phi_05" "--theta 0.4 --dtime 0.05 --updateratio 0.5" 3
# ../benchmark.sh "phi_075" "--theta 0.4 --dtime 0.05 --updateratio 0.75" 3
# ../benchmark.sh "phi_10" "--theta 0.4 --dtime 0.05 --updateratio 1.0" 3

popd;