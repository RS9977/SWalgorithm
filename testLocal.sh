gcc -O1 SWalgo.c -lgomp -o SWalgo

gcc -O1 SWalgo_V2.c -lgomp -o SWalgo_V2

gcc -O1 -mavx2 SWalgo_V4.c -lgomp -o SWalgo_V4

./SWalgo 50000 60000
./SWalgo_V2 50000 60000
./SWalgo_V4 50000 60000