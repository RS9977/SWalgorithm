core=32

gcc -mavx2 -O3 SWalgo_V9.c -lgomp -lpthread -o SWalgo_V9

./SWalgo_V9 4096 "${core}0000" "${core}"