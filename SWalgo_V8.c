/***********************************************************************
 * Smith–Waterman algorithm
 * Purpose:     Local alignment of nucleotide or protein sequences
 * Authors:     Daniel Holanda, Hanoch Griner, Taynara Pinheiro
 ***********************************************************************/
//This is the working AVX256 for 4B integer

//gcc -mavx2 SWalgo_V8.c -lgomp -lpthread -o SWalgo_V8
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include <omp.h>
#include <xmmintrin.h>
#include <smmintrin.h>
#include <immintrin.h>
#include <stdbool.h>
#include <sys/stat.h>
#include <pthread.h>
#include <time.h>
//#include <zmmintrin.h>

/*--------------------------------------------------------------------
 * Text Tweaks
 */
#define RESET   "\033[0m"
#define BOLDRED "\033[1m\033[31m"      /* Bold Red */
/* End of text tweaks */

/*--------------------------------------------------------------------
 * Constants
 */
#define PATH -1
#define NONE 0
#define UP 1
#define LEFT 2
#define DIAGONAL 3

#define Vsize 8

#define CPNS 5.0 
#define NumOfTest 1e2//1e4
#define NumOfThs 20
//#define DEBUG
//#define pragmas
/* End of constants */


/*--------------------------------------------------------------------
 * Functions Prototypes
 */
void similarityScore(long long int ind, long long int ind_u, long long int ind_d, long long int ind_l, long long int ii, long long int jj, int* H, int* P, long long int max_len, long long int* maxPos, long long int *maxPos_max_len);
void similarityScoreIntrinsic(__m256i* H,__m256i* Hu,__m256i* Hd,__m256i* Hl,__m256i* P,__m256i ii,__m256i jj, int* H_main, long long int ind, long long int max_len, long long int* maxPos, long long int *maxPos_max_len);
int  matchMissmatchScore(long long int i, long long int j);
void* enitre_simiraity_pth_worker(void* in);
void backtrack(int* P, long long int maxPos, long long int maxPos_max_len);
void printMatrix(int* matrix);
void printPredecessorMatrix(int* matrix);
void generate(void);
/* End of prototypes */
double interval(struct timespec start, struct timespec end)
{
  struct timespec temp;
  temp.tv_sec = end.tv_sec - start.tv_sec;
  temp.tv_nsec = end.tv_nsec - start.tv_nsec;
  if (temp.tv_nsec < 0) {
    temp.tv_sec = temp.tv_sec - 1;
    temp.tv_nsec = temp.tv_nsec + 1000000000;
  }
  return (((double)temp.tv_sec) + ((double)temp.tv_nsec)*1.0e-9);
}
double wakeup_delay()
{
  double meas = 0; int j, i;
  struct timespec time_start, time_stop;
  double quasi_random = 0;
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time_start);
  j = 10000000;
  while (meas < 1.0) {
    for (i=1; i<j; i++) {
      /* This iterative calculation uses a chaotic map function, specifically
         the complex quadratic map (as in Julia and Mandelbrot sets), which is
         unpredictable enough to prevent compiler optimisation. */
      quasi_random = quasi_random*quasi_random - 1.923432;
    }
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time_stop);
    meas = interval(time_start, time_stop);
    j *= 2; /* Twice as much delay next time, until we've taken 1 second */
  }
  return quasi_random;
}
/*--------------------------------------------------------------------
 * Global Variables
 */
//Defines size of strings to be compared
/*--------------------------------------------------------------------
 * Global Variables
 */
//Defines size of strings to be compared
long long int m  = 3; //Columns - Size of string a
long long int mm = 11;  //Lines - Size of string b
long long int n  = 7;  //Lind for each chunck

//Defines scores
int matchScore = 5;
int missmatchScore = -3;
int gapScore = -4; 

//Strings over the Alphabet Sigma
char *a, *b;

/* End of global variables */

double overT = 0;
pthread_mutex_t lock;

typedef struct
{
    long long int start_ind;
    int* H;
}similarity_data;

/*--------------------------------------------------------------------
 * Function:    main
 */
int main(int argc, char* argv[]) {
    
    int NumOfThreads = NumOfThs;
    if(argc>1){
        mm = atoi(argv[1]);
        n = atoi(argv[2]); 
        NumOfThreads = atoi(argv[3]); 
        int temp;
        if( mm<n){
            temp = mm;
            mm = n;
            n = temp;
        }
    }
    else{
        mm = 1000;
        n = 100;
        NumOfThreads =1;
    }

    struct timespec time_start, time_stop;
    #ifdef DEBUG
    printf("\nMatrix[%lld][%lld]\n", n, m);
    #endif

    //Allocates a and b
    a = malloc(mm * sizeof(char));
    b = malloc(n * sizeof(char));
    //Because now we have zeros
    m = mm/NumOfThreads;
    m++;
    n++;
    
    //Allocates similarity matrix H
    int *H;
    H = calloc(NumOfThreads * m * n, sizeof(int));

    //Allocates predecessor matrix P
    int *P;
    P = calloc(NumOfThreads * m * n, sizeof(int));


    //Gen rand arrays a and b
    generate();
    // a[0] =   'C';
    // a[1] =   'G';
    // a[2] =   'T';
    // a[3] =   'G';
    // a[4] =   'A';
    // a[5] =   'A';
    // a[6] =   'T';
    // a[7] =   'T';
    // a[8] =   'C';
    // a[9] =   'A';
    // a[10] =  'T';

    // b[0] =   'G';
    // b[1] =   'A';
    // b[2] =   'C';
    // b[3] =   'T';
    // b[4] =   'T';
    // b[5] =   'A';
    // b[6] =   'C';



    //Start position for backtrack
    long long int maxPos         = 0;
    long long int maxPos_max_len = 0;

    //Calculates the similarity matrix
    int i, j;

    double ww = wakeup_delay();
    pthread_mutex_init(&lock, NULL);
    clock_gettime(CLOCK_REALTIME, &time_start);
    double initialTime;// = omp_get_wtime();
    #ifdef pragmas
    #pragma GCC ivdep
    #endif


    #ifdef DEBUG
    printf("\n a string:\n");
    for(i=0; i<m-1; i++)
        printf("%c ",a[i]);
    printf("\n b string:\n");
    for(i=0; i<n-1; i++)
        printf("%c ",b[i]);
    printf("\n");
    #endif
     double t1, t2, overall=100000000000;
    int it;
    double finalTime;
    
    for(it=0; it<NumOfTest; it++){
        initialTime = omp_get_wtime();
        pthread_t threads[NumOfThreads];
        similarity_data thread_data_array[NumOfThreads];
        int t;
        int rc;          
        for(t=0; t<NumOfThreads; t++){ 
            thread_data_array[t].start_ind = t*m;  
            thread_data_array[t].H         = (H+(t*m*n));        
            rc = pthread_create(&threads[t], NULL, enitre_simiraity_pth_worker, (void*) &thread_data_array[t]);
            if (rc) {
                printf("ERROR; return code from pthread_create() is %d\n", rc);
                exit(-1);
            }
        }
        
        for (t = 0; t<NumOfThreads; t++) {
            if (pthread_join(threads[t], NULL)){ 
                printf("ERROR; code on return from join is %d\n", rc);
                exit(-1);
            }
        }
        finalTime = omp_get_wtime();
        double temp = finalTime-initialTime;
        if (temp<overall)
            overall = temp;
    }
    
    //Gets final time
    
    clock_gettime(CLOCK_REALTIME, &time_stop);
    printf("\nww:%f Best Function time V8 for n(%lld), m(%lld) and threads(%d): %f", ww, n-1, m-1, overall, NumOfThreads);
    printf("\nGCUPS max: %f", 1e-9*(m-1)*(n-1)*NumOfThreads/overall);
    double mean_time = interval(time_start, time_stop)/NumOfTest;
    printf("\nClock wall: %f\n", mean_time);
    printf("\nGCUPS mean: %f\n", 1e-9*(m-1)*(n-1)*NumOfThreads/mean_time);
    
    FILE *fp;
    fp = fopen("Results.txt", "a");
    fprintf(fp, "\nElapsed time V8 for n(%lld), m(%lld) and threads(%d): %f\n", n-1, m-1, ((finalTime - initialTime)-overall)/NumOfTest, NumOfThreads);
    fclose(fp);


    backtrack(P, maxPos, maxPos_max_len);
    #ifdef DEBUG
    printf("\nSimilarity Matrix:\n");
    printMatrix(H);

    printf("\nPredecessor Matrix:\n");
    printPredecessorMatrix(P);
    #endif

    pthread_mutex_destroy(&lock);

    //Frees similarity matrixes
    free(H);
    free(P);

    //Frees input arrays
    free(a);
    free(b);

    return 0;
}  /* End of main */


/*--------------------------------------------------------------------
 * Function:    SimilarityScore
 * Purpose:     Calculate  the maximum Similarity-Score H(i,j)
 */
void similarityScore(long long int ind, long long int ind_u, long long int ind_d, long long int ind_l, long long int ii, long long int jj, int* H, int* P, long long int max_len, long long int* maxPos, long long int *maxPos_max_len) {

    int up, left, diag;

    //Get element above
    up = H[ind_u] + gapScore;

    //Get element on the left
    left = H[ind_l] + gapScore;

    //Get element on the diagonal
    diag = H[ind_d] + matchMissmatchScore(ii, jj);

    //Calculates the maximum
    int max = NONE;
    //int pred = NONE;
    /* === Matrix ===
     *      a[0] ... a[n] 
     * b[0]
     * ...
     * b[n]
     *
     * generate 'a' from 'b', if '←' insert e '↑' remove
     * a=GAATTCA
     * b=GACTT-A
     * 
     * generate 'b' from 'a', if '←' insert e '↑' remove
     * b=GACTT-A
     * a=GAATTCA
    */
    
    if (diag > max) { //same letter ↖
        max = diag;
       // pred = DIAGONAL;
    }

    if (up > max) { //remove letter ↑ 
        max = up;
      //  pred = UP;
    }
    
    if (left > max) { //insert letter ←
        max = left;
    //    pred = LEFT;
    }
    //Inserts the value in the similarity and predecessor matrixes
    H[ind] = max;
   /* P[ind] = pred;

    //Updates maximum score to be used as seed on backtrack 
    if (max > H[*maxPos]) {
        *maxPos = ind;
        *maxPos_max_len = max_len;
    }*/

}  /* End of similarityScore */

void* enitre_simiraity_pth_worker(void* in){

    similarity_data *inss = (similarity_data *) in;
    long long int start_ind = inss -> start_ind;
    int * H                 = inss -> H;
    int * P;
    //Allocates similarity matrix H

    //Start position for backtrack
    long long int maxPos         = 0;
    long long int maxPos_max_len = 0;

    //Calculates the similarity matrix
    long long int i, j;

    double tt1=0, tt2=0, overall=0;
    double t1=0, t2=0, overal=0;

    long long int ind   = 3;
    long long int indd  = 0;
    long long int indul = 1;
    long long int ind_u, ind_d, ind_l; 
   __m256i offset =_mm256_setr_epi32(0, 1, 2, 3, 4, 5, 6, 7);
    for (i = 2; i < m+n-1; i++) { //Lines
        long long int max_len;
        long long int ii,jj;
        long long int j_start, j_end;
        if (i<n){
		    max_len = i+1;
            j_start = 1;
            j_end   = max_len-1;
            ind_u   = ind - max_len;
            ind_l   = ind - max_len + 1;
            ind_d   = ind - (max_len<<1) + 2;
        }
	    else if (i>=m){
		    max_len = m+n-1-i;
            j_start = 0;   
            j_end   = max_len;     
            ind_u   = ind - max_len - 1;
            ind_l   = ind - max_len; 
            ind_d   = ind - (max_len<<1) - 2;
        }
	    else{
		    max_len   = n;
            j_start = 1;
            j_end   = max_len;
            ind_u     = ind - max_len - 1;
            ind_l     = ind - max_len; 
            if(i>n)
                ind_d = ind - (max_len<<1) - 1;
            else
                ind_d = ind - (max_len<<1);
        }  
       
        __m256i* Hu = (__m256i*) (H+ind_u+j_start);
        __m256i* Hl = (__m256i*) (H+ind_l+j_start);
        __m256i* Hd = (__m256i*) (H+ind_d+j_start);
        __m256i* HH = (__m256i*) (H+ind+j_start);
        __m256i* PP = (__m256i*) (P+ind+j_start);
        /*uintptr_t addr = (uintptr_t)Hu;
        int offs = addr & 0x1f;
        if (offs != 0) {
            Hu = (__m256i*)((uintptr_t)Hu - offs + 32);
            Hl = (__m256i*)((uintptr_t)Hl - offs + 32);
            Hd = (__m256i*)((uintptr_t)Hd - offs + 32);
            HH = (__m256i*)((uintptr_t)HH - offs + 32);
            PP = (__m256i*)((uintptr_t)PP - offs + 32);
        }*/
        //int Vsize = 256/sizeof(typeof(H));
      //  #pragma gcc ivdep 
        for (j = j_start; j <j_end-Vsize+1; j+=Vsize) { //Columns  

           __m256i Joffset = _mm256_add_epi32(offset,_mm256_set1_epi32(j));
           __m256i Ioffset = _mm256_set1_epi32(i);
           __m256i mask    = _mm256_set1_epi32(-(int)(m>i));
           __m256i I_J     = _mm256_sub_epi32(Ioffset, Joffset);
           __m256i M_1     = _mm256_set1_epi32(m-1);
           __m256i M_1_J   = _mm256_sub_epi32(M_1,Joffset);
           __m256i II      = _mm256_blendv_epi8(M_1_J, I_J, mask);
           __m256i IJ      = _mm256_add_epi32(Joffset, Ioffset);
           __m256i I_MJ1   = _mm256_sub_epi32(IJ,M_1);
           __m256i JJ      = _mm256_blendv_epi8(I_MJ1, Joffset, mask);
                II         = _mm256_add_epi32(II, _mm256_set1_epi32(start_ind));
           similarityScoreIntrinsic(HH, Hu, Hd, Hl, PP, II, JJ, H, ind+j, max_len, &maxPos, &maxPos_max_len);
           Hu++;
           Hl++;
           Hd++;
           HH++;
           PP++;
        }
       
        for(;j<j_end; j++){
            if (i<m){
                ii = i-j+start_ind;
                jj = j;
            }
            else{
                ii = m-1-j+start_ind;
                jj = i-m+j+1;
            }      
            similarityScore(ind+j, ind_u+j, ind_d+j, ind_l+j, ii, jj, H, P, max_len, &maxPos, &maxPos_max_len);
        }
        ind += max_len;
    }
}


void similarityScoreIntrinsic(__m256i* HH,__m256i* Hu,__m256i* Hd,__m256i* Hl,__m256i* PP,__m256i ii,__m256i jj, int* H, long long int ind, long long int max_len, long long int* maxPos, long long int *maxPos_max_len) {

   
   __m256i up, left, diag;

    __m256i HHu = _mm256_loadu_si256(Hu);
    __m256i HHd = _mm256_loadu_si256(Hd);
    __m256i HHl = _mm256_loadu_si256(Hl);

    //Get element above
    up                     =_mm256_add_epi32(HHu,_mm256_set1_epi32(gapScore));

    //Get element on the left
    left                   =_mm256_add_epi32(HHl,_mm256_set1_epi32(gapScore));

    //Get element on the diagonal
    __m256i A              =_mm256_i32gather_epi32((void*) a, _mm256_sub_epi32(ii,_mm256_set1_epi32(1)), sizeof(char));
    __m256i B              =_mm256_i32gather_epi32((void*) b, _mm256_sub_epi32(jj,_mm256_set1_epi32(1)), sizeof(char));
   A                      = _mm256_slli_epi32(A,24);
   B                      = _mm256_slli_epi32(B,24);
   __m256i mask           = _mm256_cmpeq_epi32(A, B);

   __m256i MATCHSCORE     =_mm256_set1_epi32(matchScore);
   __m256i MISSMATCHSCORE =_mm256_set1_epi32(missmatchScore);
  __m256i MATCHMISS       = _mm256_blendv_epi8(MISSMATCHSCORE, MATCHSCORE, mask);
    diag                  =_mm256_add_epi32(HHd, MATCHMISS);


    //Calculates the maximum
   __m256i max  =_mm256_set1_epi32(NONE);

    /* === Matrix ===
     *      a[0] ... a[n] 
     * b[0]
     * ...
     * b[n]
     *
     * generate 'a' from 'b', if '←' insert e '↑' remove
     * a=GAATTCA
     * b=GACTT-A
     * 
     * generate 'b' from 'a', if '←' insert e '↑' remove
     * b=GACTT-A
     * a=GAATTCA
    */
   //same letter ↖
    mask    = _mm256_cmpgt_epi32(diag, max);
    max     = _mm256_blendv_epi8(max, diag, mask);

    //remove letter ↑ 
    mask    = _mm256_cmpgt_epi32(up, max);
    max     = _mm256_blendv_epi8(max, up, mask);

    //insert letter ←
    mask    = _mm256_cmpgt_epi32(left, max);
    max     = _mm256_blendv_epi8(max, left, mask);


    //Inserts the value in the similarity and predecessor matrixes
    _mm256_storeu_si256(HH, max);

}  /* End of similarityScore */


/*--------------------------------------------------------------------
 * Function:    matchMissmatchScore
 * Purpose:     Similarity function on the alphabet for match/missmatch
 */
int matchMissmatchScore(long long int i, long long int j) {
    if (a[i-1] == b[j-1])
        return matchScore;
    else
        return missmatchScore;
}  /* End of matchMissmatchScore */

/*--------------------------------------------------------------------
 * Function:    backtrack
 * Purpose:     Modify matrix to print, path change from value to PATH
 */
void backtrack(int* P, long long int maxPos, long long int maxPos_max_len) {
    //hold maxPos value
    long long int predPos;
    long long int predMaxLen;
    #ifdef pragmas
    #pragma GCC ivdep
    #endif
    //backtrack from maxPos to startPos = 0 
    long long int first_sec = (n*(n+1))/2;
    long long int last_sec  = n*m - (n*(n-1))/2;
    long long int ind_u, ind_d, ind_l; 
    bool diagCompensate     = 0;
    do {
        if (maxPos<first_sec){
            if(diagCompensate){
                if(maxPos<first_sec-n)
                    maxPos_max_len --;
                diagCompensate = 0;
            }
            ind_u      = maxPos - maxPos_max_len;
            ind_l      = maxPos - maxPos_max_len + 1;
            ind_d      = maxPos - (maxPos_max_len<<1) + 2;
            predMaxLen = maxPos_max_len-1;
        }
	    else if (maxPos>=last_sec){
            if(diagCompensate){
                maxPos_max_len ++;
                diagCompensate = 0;
            }
            ind_u      = maxPos - maxPos_max_len - 1;
            ind_l      = maxPos - maxPos_max_len; 
            ind_d      = maxPos - (maxPos_max_len<<1) - 2;
            predMaxLen = maxPos_max_len+1;
        }
	    else{
            if(diagCompensate){
                if(maxPos>=last_sec-n)
                    maxPos_max_len ++;
                diagCompensate = 0;
            }
            ind_u      = maxPos - n - 1;
            ind_l      = maxPos - n;
            predMaxLen = maxPos_max_len; 
            if(maxPos>=first_sec+n)
                ind_d  = maxPos - (n<<1) - 1;
            else
                ind_d  = maxPos - (n<<1);
        }

        if(P[maxPos] == DIAGONAL){
            predPos        = ind_d;
            diagCompensate = 1;
        }
        else if(P[maxPos] == UP)
            predPos        = ind_u;
        else if(P[maxPos] == LEFT)
            predPos        = ind_l;
        
        P[maxPos]*=PATH;
        maxPos         = predPos;
        maxPos_max_len = predMaxLen;
    } while(P[maxPos] != NONE);
}  /* End of backtrack */


/*--------------------------------------------------------------------
 * Function:    printMatrix
 * Purpose:     Print Matrix
 */
void printMatrix(int* matrix) {
    long long int i, j, ind;
    printf(" \t \t");
    for(i=0; i<m-1; i++)
        printf("%c\t",a[i]);
    printf("\n");
    for (i = 0; i < n; i++) { //Lines
        for (j = -1; j < m; j++) {
            if(i+j<n)
                ind = (i+j)*(i+j+1)/2 + i;
            else if(i+j<m)
                ind = (n+1)*(n)/2 + (i+j-n)*n + i;
            else
                ind = (i*j) + ((m-j)*(i+(i-(m-j-1))))/2 + ((n-i)*(j+(j-(n-i-1))))/2 + (m-j-1);
            if(i+j<0)
                printf(" \t");
            else if(j==-1 && i>0)
                printf("%c\t",b[i-1]); 
            else
                printf("%d\t", matrix[ind]);
        }
        printf("\n");
    }
}  /* End of printMatrix */

/*--------------------------------------------------------------------
 * Function:    printPredecessorMatrix
 * Purpose:     Print predecessor matrix
 */
void printPredecessorMatrix(int* matrix) {
    long long int i, j, ind;
    printf("    ");
    for(i=0; i<m-1; i++)
        printf("%c ",a[i]);
    printf("\n");
    for (i = 0; i < n; i++) { //Lines
        for (j = -1; j < m; j++) {
            if(i+j<n)
                ind = (i+j)*(i+j+1)/2 + i;
            else if(i+j<m)
                ind = (n+1)*(n)/2 + (i+j-n)*n + i;
            else
                ind = (i*j) + ((m-j)*(i+(i-(m-j-1))))/2 + ((n-i)*(j+(j-(n-i-1))))/2 + (m-j-1);
            if(i+j<0)
                printf("  ");
            else if(j==-1 && i>0)
                printf("%c ",b[i-1]); 
            else{
                if(matrix[ind] < 0) {
                    printf(BOLDRED);
                    if (matrix[ind] == -UP)
                        printf("↑ ");
                    else if (matrix[ind] == -LEFT)
                        printf("← ");
                    else if (matrix[ind] == -DIAGONAL)
                        printf("↖ ");
                    else
                        printf("- ");
                    printf(RESET);
                } else {
                    if (matrix[ind] == UP)
                        printf("↑ ");
                    else if (matrix[ind] == LEFT)
                        printf("← ");
                    else if (matrix[ind] == DIAGONAL)
                        printf("↖ ");
                    else
                        printf("- ");
            }
            }
        }
        printf("\n");
    }

}  /* End of printPredecessorMatrix */


/*--------------------------------------------------------------------
 * Function:    generate
 * Purpose:     Generate arrays a and b
 */
 /*--------------------------------------------------------------------
 * Function:    generate
 * Purpose:     Generate arrays a and b
 */
 void generate(){
    //Generates the values of a
    int i;
    for(i=0;i<m;i++){
        int aux=rand()%24;
        if(aux==0)
            a[i]='A';
        else if(aux==1)
            a[i]='R';
        else if(aux==2)
            a[i]='N';
        else if(aux==3)
            a[i]='D';
        else if(aux==4)
            a[i]='C';
        else if(aux==5)
            a[i]='Q';
        else if(aux==6)
            a[i]='E';
        else if(aux==7)
            a[i]='G';
        else if(aux==8)
            a[i]='H';
        else if(aux==9)
            a[i]='I';
        else if(aux==10)
            a[i]='L';
        else if(aux==11)
            a[i]='K';
        else if(aux==12)
            a[i]='M';
        else if(aux==13)
            a[i]='F';
        else if(aux==14)
            a[i]='P';
        else if(aux==15)
            a[i]='S';
        else if(aux==16)
            a[i]='T';
        else if(aux==17)
            a[i]='W';
        else if(aux==18)
            a[i]='Y';
        else if(aux==19)
            a[i]='V';
        else if(aux==20)
            a[i]='B';
        else if(aux==21)
            a[i]='J';
        else if(aux==22)
            a[i]='Z';
        else 
            a[i]='X';
    }

    //Generates the values of b
    for(i=0;i<n;i++){
        int aux=rand()%24;
        if(aux==0)
            b[i]='A';
        else if(aux==1)
            b[i]='R';
        else if(aux==2)
            b[i]='N';
        else if(aux==3)
            b[i]='D';
        else if(aux==4)
            b[i]='C';
        else if(aux==5)
            b[i]='Q';
        else if(aux==6)
            b[i]='E';
        else if(aux==7)
            b[i]='G';
        else if(aux==8)
            b[i]='H';
        else if(aux==9)
            b[i]='I';
        else if(aux==10)
            b[i]='L';
        else if(aux==11)
            b[i]='K';
        else if(aux==12)
            b[i]='M';
        else if(aux==13)
            b[i]='F';
        else if(aux==14)
            b[i]='P';
        else if(aux==15)
            b[i]='S';
        else if(aux==16)
            b[i]='T';
        else if(aux==17)
            b[i]='W';
        else if(aux==18)
            b[i]='Y';
        else if(aux==19)
            b[i]='V';
        else if(aux==20)
            b[i]='B';
        else if(aux==21)
            b[i]='J';
        else if(aux==22)
            b[i]='Z';
        else 
            b[i]='X';
    }
} /* End of generate */



/*--------------------------------------------------------------------
 * External References:
 * http://vlab.amrita.edu/?sub=3&brch=274&sim=1433&cnt=1
 * http://pt.slideshare.net/avrilcoghlan/the-smith-waterman-algorithm
 * http://baba.sourceforge.net/
 */