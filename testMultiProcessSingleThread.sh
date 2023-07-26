#!/bin/bash -l


NUM_CORES=32

lscpu

module load gcc/12

gcc -mavx2 -O3 SWalgo_V4.c -lgomp -o SWalgo_V4

# Number of cores to utilize


# Path to the binary executable
BINARY_PATH="./SWalgo_V4 100 1000"

# Function to execute the binary on a specific core
run_on_core() {
  core=$1
  taskset -c $core $BINARY_PATH
}

echo "------------------------------------ 32 -----------------------------"

start_time=$(date +%s.%N)

# Loop through the number of cores and run the binary on each core
for ((i = 0; i < NUM_CORES; i++)); do
  run_on_core $i &
done

# Wait for all processes to finish
wait


# End time
end_time=$(date +%s.%N)

# Calculate the total execution time
execution_time=$(echo "$end_time - $start_time" | bc)

echo "Total execution time: $execution_time seconds"



echo "------------------------------------ 16 -----------------------------"

start_time=$(date +%s.%N)

# Loop through the number of cores and run the binary on each core
for ((i = 0; i < 16; i++)); do
  run_on_core $i &
done

# Wait for all processes to finish
wait


# End time
end_time=$(date +%s.%N)

# Calculate the total execution time
execution_time=$(echo "$end_time - $start_time" | bc)

echo "Total execution time: $execution_time seconds"


echo "------------------------------------ 8 -----------------------------"

start_time=$(date +%s.%N)

# Loop through the number of cores and run the binary on each core
for ((i = 0; i < 8; i++)); do
  run_on_core $i &
done

# Wait for all processes to finish
wait


# End time
end_time=$(date +%s.%N)

# Calculate the total execution time
execution_time=$(echo "$end_time - $start_time" | bc)

echo "Total execution time: $execution_time seconds"


echo "------------------------------------ 4 -----------------------------"

start_time=$(date +%s.%N)

# Loop through the number of cores and run the binary on each core
for ((i = 0; i < 4; i++)); do
  run_on_core $i &
done

# Wait for all processes to finish
wait


# End time
end_time=$(date +%s.%N)

# Calculate the total execution time
execution_time=$(echo "$end_time - $start_time" | bc)

echo "Total execution time: $execution_time seconds"


echo "------------------------------------ 2 -----------------------------"

start_time=$(date +%s.%N)

# Loop through the number of cores and run the binary on each core
for ((i = 0; i < 2; i++)); do
  run_on_core $i &
done

# Wait for all processes to finish
wait


# End time
end_time=$(date +%s.%N)

# Calculate the total execution time
execution_time=$(echo "$end_time - $start_time" | bc)

echo "Total execution time: $execution_time seconds"

echo "------------------------------------ 1 -----------------------------"

start_time=$(date +%s.%N)

# Loop through the number of cores and run the binary on each core
for ((i = 0; i < 1; i++)); do
  run_on_core $i &
done

# Wait for all processes to finish
wait


# End time
end_time=$(date +%s.%N)

# Calculate the total execution time
execution_time=$(echo "$end_time - $start_time" | bc)

echo "Total execution time: $execution_time seconds"