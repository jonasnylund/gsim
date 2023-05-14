#!/bin/bash

./benchmark.sh "exact" "--updateratio 0.0 --theta 0.2 --dtime 0.001" 1

# ./benchmark.sh "fast" "--updateratio 0.414 --theta 1.0 --dtime 0.02" 5
# ./benchmark.sh "baseline" "--updateratio 0.0 --theta 0.5 --dtime 0.01" 10

# ./benchmark.sh "reduced_phi" "--updateratio 0.414 --theta 0.5 --dtime 0.01" 2
# ./benchmark.sh "reduced_theta" "--updateratio 0.0 --theta 1.0 --dtime 0.01" 10
# ./benchmark.sh "reduced_dt" "--updateratio 0.0 --theta 0.5 --dtime 0.02" 10

# ./benchmark.sh "balanced2" "--updateratio 0.25 --theta 0.5 --dtime 0.007" 10

# ./benchmark.sh "phi_00" "--theta 0.25 --dtime 0.02 --updateratio 0.0" 5
# ./benchmark.sh "phi_01" "--theta 0.25 --dtime 0.02 --updateratio 0.1" 5
# ./benchmark.sh "phi_02" "--theta 0.25 --dtime 0.02 --updateratio 0.2" 5
# ./benchmark.sh "phi_03" "--theta 0.25 --dtime 0.02 --updateratio 0.3" 5
# ./benchmark.sh "phi_04" "--theta 0.25 --dtime 0.02 --updateratio 0.4" 5
# ./benchmark.sh "phi_05" "--theta 0.25 --dtime 0.02 --updateratio 0.5" 5
# ./benchmark.sh "phi_075" "--theta 0.25 --dtime 0.02 --updateratio 0.75" 5
# ./benchmark.sh "phi_10" "--theta 0.25 --dtime 0.02 --updateratio 1.0" 5