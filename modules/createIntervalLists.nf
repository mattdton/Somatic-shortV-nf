#!/usr/bin/env nextflow

process createIntervalLists {

// where to publish the outputs
tag "createIntervalLists"
publishDir "${params.outDir}", mode:'copy'

input:
    val(number_of_intervals)

output:
    path("BQSR_intervals")
    path("intervals_list.txt")
    env (intervalList)

shell:
// Create an intervalList based on user specified number_of_intervals
'''
for ((i = 0; i <= !{number_of_intervals-1}; i++)); do
    printf "%04d\n" "$i" >> intervals_list.txt
done


# Read the content of intervals_list.txt into a Bash array named intervalListArray
IFS=$'\n'
intervalList=()
while read line; do
  intervalList+=("$line")
done < "intervals_list.txt"




# Access elements using indexing
echo "First element: ${intervalList[0]}"

gatk SplitIntervals \
        --java-options "-Xmx3g" \
        -R !{params.ref} \
        -scatter-count !{number_of_intervals} \
        -O BQSR_intervals

#-XL !{params.exclude_intervals} \



'''




}