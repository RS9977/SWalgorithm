#include <immintrin.h>
#include <stdio.h>

int main() {
    // Input 16-bit integers
    int i;
    for(i=0; i<10; i++){
        if(i%4==0)
            i++;
        printf("%d\n", i);
    }

    return 0;
}
