#!/bin/bash
# example MPI+OpenMP job script for SLURM
#
# Tell SLURM which project's account to use:
#SBATCH -A abnffidp
#
#
# SLURM defaults to the directory you were working in when you submitted the job.
# Output files are also put in this directory. To set a different working directory add:
#
#
#SBATCH --mail-type=ALL
#SBATCH --mail-user=b8048283@ncl.ac.uk
# Tell SLURM if you want to be emailed when your job starts, ends, etc.
# Currently mail can only be sent to addresses @ncl.ac.uk
#
#
# This example has 4 MPI tasks, each with 22 cores
#
# number of tasks
#SBATCH  --ntasks=1
#SBATCH -t 48:00:00
module load GROMACS/2020.5-foss-2020b-Python-3.8.6
#
#
# SLURM recommend using srun instead of mpirun for better job control.
# You need to have a pre-compiled obabel for this script
lig_name=MOL # The ligand residue name, e.g UNK, MOL
skip=1       # The value used in gmx trjconv -skip $skip
hydescore=/mnt/nfs/home/b8055696/workspace/hydescorer-1.5.0-Linux-x64/hydescorer # the excutable file path for hydescore
rm order.dat affinity.dat top1 number.dat output.dat output1.dat path.dat hyde_score.dat
rm -r hyde_analysis
mkdir hyde_analysis
group=`printf "q\n" | gmx make_ndx -f md.gro | grep $lig_name | head -1 | awk '{print $1}'`
printf "1 | $group\n q\n" | gmx make_ndx -f md.gro -o rec_lig.ndx
group=`printf "q\n" | gmx make_ndx -f md.gro -n rec_lig.ndx | grep Protein_$lig_name | tail -1 | awk '{print $1}'`
printf "$group\n $group\n" | gmx trjconv -s md.tpr -f md.xtc -center -pbc mol -n rec_lig.ndx -skip $skip -o ./hyde_analysis/md_noPBC.pdb -sep # total 5000 frames
cd hyde_analysis
num=0
number=`ls | wc -l`

for i in md_noPBC*; do mkdir ${i%.*}; mv $i ${i%.*}; done
echo "#!/bin/bash
# example MPI+OpenMP job script for SLURM
#
# Tell SLURM which project's account to use:
#SBATCH -A abnffidp
#
#
# SLURM defaults to the directory you were working in when you submitted the job.
# Output files are also put in this directory. To set a different working directory add:
#
#
#SBATCH --mail-type=NONE
#SBATCH --mail-user=b8048283@ncl.ac.uk
# Tell SLURM if you want to be emailed when your job starts, ends, etc.
# Currently mail can only be sent to addresses @ncl.ac.uk
#
#
# This example has 4 MPI tasks, each with 22 cores
#
# number of tasks
#SBATCH  --ntasks=1
#SBATCH -t 06:00:00

module load CMake/3.20.1-GCCcore-10.2.0
# module load CMake/3.23.1-GCCcore-11.3.0

more md_noPBC*pdb | grep $lig_name > lig.pdb; more md_noPBC*.pdb  | grep -v $lig_name > rec.pdb; obabel -iPDB lig.pdb -O lig.mol2 > /dev/null 2>&1; $hydescore -i lig.mol2  -o scored_lig.sdf -p rec.pdb -r lig.mol2  -v 4" > submit_mmhyde.sh
for i in `seq 0 $number`; do num=$((num+1)); cp submit_mmhyde.sh md_noPBC$i/${i}_mmhyde.sh; done
for i in `seq 0 $number`; do num=$((num+1)); cd  md_noPBC$i; sbatch ${i}_mmhyde.sh; cd ..; done






