# readme

file name format is 

`cfg-sim<sim-id>-<scenario-id>-<variant-id>`

use `sed` for efficient updating, e.g. 


+ Change number of simulations from 10 to 500:

`gsed -i 's/nsim:\ 10/nsim:\ 500/' cfg-sim01-sc01-0*.yml`

`gsed -i 's/nsim:\ 10/nsim:\ 500/' cfg-sim02-sc01-v0*.yml`

+ Update scenario label:

`gsed -i 's/sc01/sc02/' cfg-sim01-sc02-0*.yml`
