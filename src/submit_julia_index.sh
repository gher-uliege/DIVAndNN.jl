#!/bin/bash
#SBATCH --job-name=julia
######SBATCH --mail-user=a.barth@ulg.ac.be
######SBATCH --mail-type=ALL
#SBATCH --output=output-%j-%N.out
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4000

# Use as
# sbatch --mem-per-cpu=2000 --time=12:00:00  submit_julia_index.sh  ~/src/DIVAndNN/src/emodnet_bio3.jl 1
#
# where 1 is the index of the specie (1 to 40)

module load  EasyBuild  Python/3.5.1-foss-2016a

export script="$1"
export INDEX="$2"
echo script $script

bt0=$(date +%s)

unset DISPLAY
#printenv

julia <<EOF
include("$script")
EOF

bt1=$(date +%s)

awk  " BEGIN { print \"Run time (hours): \",($bt1 - $bt0)/3600 } "



