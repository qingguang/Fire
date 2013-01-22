#!/bin/bash
#SBATCH -o $HOME/myjob.%j.%N.out
#SBATCH -D $HOME
#SBATCH -J POS_xie
#SBATCH --get-user-env
#SBATCH --clusters=uv2
#SBATCH --ntasks=8
#SBATCH --mail-type=end
#SBATCH --mail-user=sparkallenxie@gmail.com
#SBATCH --export=NONE
#SBATCH --time=00:05:00


source /etc/profile.d/modules.sh


cd $HOME/super/Fire/A2/A2.3/
srun_ps $HOME/super/Fire/A2/A2.3/gccg  ../pent.geo.bin dout dual

