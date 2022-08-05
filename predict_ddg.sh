#!/bin/sh
while getopts o:i:j:p:n: flag
do
    case "${flag}" in
	o) output=${OPTARG};;
	i) input=${OPTARG};;
	j) jump=${OPTARG};;
	p) positions=${OPTARG};;
	n) numtasks=${OPTARG};;
    esac
done

echo "OUTPUT FOLDER: ${output}"
echo "INPUT: ${input}"
echo "JUMP: ${jump}"
echo "POSITIONS: ${positions}"
echo "NTASKS: ${numtasks}"

id=$(sbatch --export=NONE --ntasks=${numtasks} ddgs.sh ${pdbarray[0]} ${output})
echo "submitted initial job $id for pdb ${pdbarray[0]}"
for i in "${pdbarray[@]:1}"; do
  id=$(sbatch --export=NONE --ntasks=${numtasks} ddgs_flex.sh $i ${output});
  echo "submitted job $id for pdb $i"
done
id=$(sbatch --export=NONE --ntasks=1 get_csv.sh ${output} --depend=afterany:$id)
