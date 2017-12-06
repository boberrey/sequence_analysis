#
#all commands that start with SBATCH contain commands that are just used by SLURM for scheduling  
#################
#set a job name  
#SBATCH --job-name=snakeATAC
#################  
#a file for job output, you can check job progress
#SBATCH --output=quantify_tiles_biochem.out
#################
# a file for errors from the job
#SBATCH --error=slurmSnake.err
#################
#time you think you need; default is one hour
#in minutes in this case, hh:mm:ss
#SBATCH --time=15:00:00
#################
#quality of service; think of it as job priority
#SBATCH --partition=biochem,owners,normal
#SBATCH --qos=normal
#################
#number of nodes you are requesting
#SBATCH --nodes=1
#################
#tasks to run per node; a "task" is usually mapped to a MPI processes.
# for local parallelism (OpenMP or threads), use "--ntasks-per-node=1 --cpus-per-task=16" instead
# Sherlock nodes have 16 cpus. For some reason, you can request more than that on 'owners' partition, but not on others. 
# It can't give you more obviously, but it's interesting that it lets you ask on one partion but not another.
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#################


# Define paths
fastq_dir=$1

# Define other script parameters
num_cores="16"

# Set outputs appropriately 
script=$(basename $0)
script_name=${script%.*}
log_dir=$fastq_dir/$script_name"Logs"
log_file_suffix=".log"
err_file_suffix=".err"

mkdir -p $log_dir


