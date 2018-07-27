#!/bin/bash
#!/bin/bash

#SBATCH -p Lewis
#SBATCH --account=biocommunity
#SBATCH -J GVCFgeno
#SBATCH --mem 200G
#SBATCH -N1
#SBATCH -n10
#SBATCH -t 1-00:00
#SBATCH --output=gtGVCFgeno-%A_%a-%j.out

## notifications
#SBATCH --mail-user=buckleyrm@missouri.edu  # email address for notifications
#SBATCH --mail-type=BEGIN,END,FAIL  # which type of notifications to send

#---------------------------------------------------------------------





