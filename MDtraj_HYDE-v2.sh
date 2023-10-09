#!/bin/bash
# Usage: ./MDtraj_hydescoring.sh <md.tpr>  <md.xtc>  <md.gro> <ligand_name_in_index_eg_UNK>  <skip_number_of_poses> <path_to_hydescore>

###
### MDtraj_hydescoring.sh - calculate the HYDE score for the frames generated during MD simulation
###
### Usage:
### ./MDtraj_hydescoring.sh md.tpr md.xtc md.gro UNK 50 /home/dozeduck/complex/wip1/gsk8-wip1/tc-grp/hydescorer-1.5.0-Linux-x64
###
### Options:
###     <md.tpr>				The tpr file for MD simulation 
###     <md.xtc>				The xtc file generated from MD simulation
###     <md.gro>				The gro file generated from MD simulation
###     <ligand_name_in_index_eg_UNK>	The name for ligand which showed in index.file
###     <skip_number_of_poses>		The number of poses for skipping
###     <path_to_hydescore>			e.g: /home/dozeduck/complex/wip1/gsk8-wip1/tc-grp/hydescorer-1.5.0-Linux-x64
###     -h					Show this message

help() {
	sed -rn 's/^### ?//;T;p;' "$0"
}

if [[ $# == 0 ]]  || [[ "$1" == "-h" ]]; then
	help
	exit 1
fi

rm -r hyde_analysis
mkdir hyde_analysis
group=`printf "q\n" | gmx make_ndx -f $3 | grep $4 | head -1 | awk '{print $1}'`
printf "1 | $group\n q\n" | gmx make_ndx -f $3 -o rec_lig.ndx
group=`printf "q\n" | gmx make_ndx -f $3 -n rec_lig.ndx | grep Protein_$4 | tail -1 | awk '{print $1}'`
printf "$group\n $group\n" | gmx trjconv -s $1 -f $2 -center -pbc nojump -n rec_lig.ndx -o ./hyde_analysis/md_noPBC.pdb -skip $5 -sep # total 5000 frames
cd hyde_analysis
num=0
for i in md_noPBC*; do mkdir ${i%.*}; mv $i ${i%.*}; done 	# move each single file to a specific folder
for i in md_noPBC*; do num=$(($num+1)); cd $i; more md_noPBC*pdb | grep $4 > lig.pdb; more md_noPBC*.pdb  | grep -v $4 > rec.pdb; obabel -iPDB lig.pdb -O lig.mol2 > /dev/null 2>&1; $6/hydescorer -i lig.mol2  -o scored_lig.sdf -p rec.pdb -r lig.mol2  -v 4; echo "This is No.$num, there is still $((`ls .. | grep md_noPBC | wc -l`-$num)) left"; cd ..; done
a=$((`ls | grep md_noPBC | wc -l`-1))
for i in `seq 0 $a`; do echo $i >> order.dat; more md_noPBC$i/scored_lig.sdf  | grep -A1 HYDE_ESTIMATED_AFFINITY_LOWER_BOUNDARY | grep -v nM >> affinity.dat; done
paste order.dat affinity.dat > hyde_score.dat
more hyde_score.dat | grep ` awk 'BEGIN{min = 999999999}{if($2 < min) min = $2}END{print min}' hyde_score.dat` | awk '{print$1}' > top1 # find the binding poses with best estimated affinity.
a=`more top1`; mv md_noPBC$a md_noPBC$a\_top1; cp md_noPBC$a\_top1/md_noPBC*pdb .
rm order.dat affinity.dat top1

