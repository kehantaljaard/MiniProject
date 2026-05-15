#!/bin/bash

runs=50
mode=$1

times=()

for i in $(seq 1 $runs); do
    t=$(
        /usr/bin/time -f "%e" \
        ./build/miniproject_q3 \
        Datasets/dataset3 \
        output/dataset3 \
        0.9 0.1 "$mode" \
        >/dev/null 2>&1
    )

    times+=("$t")

    printf "Run %d: %s s\n" "$i" "$t"
done


printf "%s\n" "${times[@]}" | awk '
{
    sum += $1
    sumsq += ($1)^2
}
END {
    mean = sum / NR
    std = sqrt((sumsq / NR) - (mean)^2)

    printf("\nMean: %.6f s\n", mean)
    printf("Std Dev: %.6f s\n", std)
}'