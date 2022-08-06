#!/bin/sh
while getopts o:i:j:p:n: flag
do
    case "${flag}" in
	o) output=${OPTARG};;
	i) input=${OPTARG};;
	j) jump=${OPTARG};;
	p) positions=(${OPTARG});;
	n) numtasks=${OPTARG};;
    esac
done

echo "OUTPUT FOLDER: ${output}"
echo "INPUT: ${input}"
echo "JUMP: ${jump}"
echo "POSITIONS: ${positions}"
echo "NTASKS: ${numtasks}"

id=$(sbatch --export=NONE --ntasks=${numtasks} ddgs_final_py.sh ${positions[0]} ${output} ${jump} ${input})
echo "submitted initial job $id for position ${positions[0]}"
for i in "${positions[@]:1}"; do
  id=$(sbatch --export=NONE --ntasks=${numtasks} ddgs_final_py.sh $i ${output} ${jump} ${input});
  echo "submitted job $id for position $i"
done
