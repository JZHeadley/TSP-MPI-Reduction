#!/bin/bash
gridSize=1000
outFile=results.csv
echo "numCities,numBlocks,numProcs,time,cost" > $outFile
for numCities in `seq 5 10`; #num cities
do
    # numCities=$(echo 2^${i} | bc)
    # echo $numCities
    for numBlocks in `seq 10 10 200`;
    do
        # echo $numBlocks
        for numProcs in `seq 2 2 20`;
        do
            # echo $numProcs
            raw_output=`mpirun -np $numProcs ./tsp $numCities $numBlocks $gridSize $gridSize | tail -n 1`
            cost=`echo $raw_output | grep -o '[0-9]*\.[0-9]*'`
            time=`echo $raw_output | grep -o '[0-9]*' | head -n 1`
            echo "${numCities},${numBlocks},${numProcs},${time},${cost}" | tee -a $outFile

        done
    done
    # output=`./threaded-tsp datasets/${numCities}.cities | tail -n 2`
    # echo ${output} | tee -a results.txt
done
