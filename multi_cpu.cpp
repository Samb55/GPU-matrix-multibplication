#include <fstream>
#include <iostream>
#include <stdio.h>
#include <string>
#include <sstream>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <ctime>
#include <vector>
#include <ctime>
#define BILLION 1000000000L

long runTime_Nsec;
double runTime_Sec;
struct timespec startTime;
struct timespec endTime;

float GFLOPS;
double numberOfOperations;
using namespace std;

struct matrix{
    unsigned int rows;
    unsigned int cols;
};


//ijk algorithm no blocking
// This part is only for checking the correctness of the blocking version algorithms
void dgemm_1 (float *a, float *b, float *c,int m,int n, int z) {
    int i,j,k;
    for (i=0; i<m; i++)
        for (j=0; j<n; j++) {
            register float r = c[i*n+j] ;
            for (k=0; k<z; k++)
                r += a[i*z+k] * b[k*z+j];
            c[i*n+j] = r;
        }
}
//ijk algorithm blocked
void dgemm_1B (float *a, float *b, float *c, int m,int n, int z , int B) {
    int i,j,k,i1,j1,k1;
    for(i=0; i<m; i+=B)
        for(j=0; j<n; j+=B)
            for(k=0; k<z; k+=B)
                for(i1=i; i1<i+B; i1++)
                    for(j1=j; j1<j+B; j1++){
                        register float r = c[i1*n+j1];
                        for(k1=k; k1<k+B; k1++)
                            r+=a[i1*z+k1]*b[k1*z+j1];
                        c[i1*n+j1]=r;
                    }
}
//jik algorithm blocked
void dgemm_2B (float *a, float *b, float *c, int m,int n, int z,int B) {
    int i,j,k,i1,j1,k1;
    for(j=0; j<m; j+=B)
        for(i=0; i<n; i+=B)
            for(k=0; k<z; k+=B)
                for(j1=j; j1<j+B; j1++)
                    for(i1=i; i1<i+B; i1++){
                        register float r = c[i1*n+j1];
                        for(k1=k; k1<k+B; k1++)
                            r+=a[i1*z+k1]*b[k1*z+j1];
                        c[i1*n+j1]=r;
                    }
}

//kij algorithm blocked
void dgemm_3B (float *a, float *b, float *c, int m,int n, int z,int B) {
    int i,j,k,i1,j1,k1;
    for(k=0; k<m; k+=B)
        for(i=0; i<n; i+=B)
            for(j=0; j<z; j+=B)
                for(k1=k; k1<k+B; k1++)
                    for(i1=i; i1<i+B; i1++){
                        register float r = a[i1*z+k1];
                        for(j1=j; j1<j+B; j1++)
                            c[i1*n+j1]+=r*b[k1*z+j1];
                    }
}
//ikj algorithm blocked
void dgemm_4B (float *a, float *b, float *c,int m,int n, int z,int B) {

    int i,j,k,i1,j1,k1;
    for(i=0; i<m; i+=B)
        for(k=0; k<n; k+=B)
            for(j=0; j<z; j+=B)
                for(i1=i; i1<i+B; i1++)
                    for(k1=k; k1<k+B; k1++){
                        register float r = a[i1*n+k1];
                        for(j1=j; j1<j+B; j1++)
                            c[i1*z+j1]+=r*b[k1*z+j1];
                    }

}
void dgemm44(float *a, float *b,float *c,int m ,int n,int z, int B)
{

    int i, j, k, i1, j1, k1;
    for (j = 0; j < m; j += B)
        for (k = 0; k < n; k += B)
            for (i = 0; i < z; i += B)
                for (j1 = j; j1 < j+B; j1++)
                    for (k1 = k; k1 < k+B; k1++)
                    {
                        register float r = b[k1*z+j1];
                        for (i1 = i; i1 < i+B; i1++)
                            c[i1*n+j1] += a[i1*z+k1] * r;
                    }

}

//jki algorithm blocked
void dgemm_5B (float *a, float *b, float *c, int m ,int n,int z,int B) {
    int i,j,k,i1,j1,k1;
    for(j=0; j<m; j+=B)
        for(k=0; k<n; k+=B)
            for(i=0; i<z; i+=B)
                for(j1=j; j1<j+B; j1++)
                    for(k1=k; k1<k+B; k1++){
                        register float r = b[k1*z+j1];
                        for(i1=i; i1<i+B; i1++)
                            c[i1*n+j1]+=r*a[i1*z+k1];
                    }
}
//kji algorithm blocked
void dgemm_6B (float *a, float *b, float *c,int m ,int n,int z,int B) {
    int i,j,k,i1,j1,k1;
    for(k=0; k<m; k+=B)
        for(j=0; j<n; j+=B)
            for(i=0; i<z; i+=B)
                for(k1=k; k1<k+B; k1++)
                    for(j1=j; j1<j+B; j1++){
                        register float r = b[k1*z+j1];
                        for(i1=i; i1<i+B; i1++)
                            c[i1*n+j1]+=r*a[i1*z+k1];
                    }
}



 
void calculateElapsedTimeAndGFLOPs(int m, int n, int k){
    runTime_Nsec = BILLION * (endTime.tv_sec-startTime.tv_sec) + endTime.tv_nsec - startTime.tv_nsec;
    runTime_Sec = runTime_Nsec * pow(10,(-9));
    numberOfOperations = 2.0 * (m*n);
    GFLOPS = numberOfOperations / runTime_Nsec;
}  
   
int main(int argc, char* argv[])
{
    if(argc != 4) {//there should be four arguments
        printf("missing argumnet");
        return 1; //exit and return an error
    }


    ifstream infile_A, infile_B;	//reading the input matrices

    // *****************************************************************************
    //                                   Matrix A                                 //
    //******************************************************************************

    infile_A.open(argv[1],ios::binary|ios::in|ios::ate);

    //getting end and beginning of the file
    infile_A.seekg(0,ios::end);
    infile_A.seekg(0,ios::beg);

    //memory allocation
    matrix M_A;
    infile_A.read(reinterpret_cast<char*>(&M_A),2*sizeof(unsigned int));


    float* array_A=(float*)malloc(M_A.rows*M_A.cols*sizeof(float));	//column major
    infile_A.read(reinterpret_cast<char*>(array_A),M_A.rows*M_A.cols);

    infile_A.close();
    //print(array_A,M_A.rows,M_A.cols);

    // *****************************************************************************
    //                                   Matrix B                                 //
    //******************************************************************************

    infile_B.open(argv[2],ios::binary|ios::in|ios::ate);

    //getting end and beginning of the file
    infile_B.seekg(0,ios::end);
    infile_B.seekg(0,ios::beg);

    //memory allocation
    matrix M_B;
    infile_B.read(reinterpret_cast<char*>(&M_B),2*sizeof(unsigned int));


    float* array_B=(float*)malloc(M_B.rows*M_B.cols*sizeof(float));	//column major
    infile_B.read(reinterpret_cast<char*>(array_B),M_B.rows*M_B.cols);
    //print(array_B,M_B.rows,M_B.cols);
    infile_B.close();
    // *****************************************************************************
    //                                  Multblication                             //
    //******************************************************************************
    float* D1=(float*)malloc(M_A.rows*M_B.cols*sizeof(float));
    float* D2=(float*)malloc(M_A.rows*M_B.cols*sizeof(float));
    float* D3=(float*)malloc(M_A.rows*M_B.cols*sizeof(float));
    float* D4=(float*)malloc(M_A.rows*M_B.cols*sizeof(float));
    float* D5=(float*)malloc(M_A.rows*M_B.cols*sizeof(float));
    float* D6=(float*)malloc(M_A.rows*M_B.cols*sizeof(float));
    float* D7=(float*)malloc(M_A.rows*M_B.cols*sizeof(float));

    clock_gettime(CLOCK_REALTIME, &startTime);

    dgemm_1(array_A, array_B, D1,M_A.rows,M_B.cols,M_A.cols);
    clock_gettime(CLOCK_REALTIME, &endTime);
    calculateElapsedTimeAndGFLOPs(M_A.rows, M_B.cols, M_A.cols);
    


 printf(" runTime = %f s \ns\tnumber of FloatingPoint Operations = %.0f\tGFLOPS = %f\n" ,runTime_Sec ,(long long unsigned int)runTime_Nsec ,numberOfOperations ,GFLOPS); 

    printf("------------------------------|ijk algorithm done no block|---------------------------------------\n");

 for(int B=8;B<=64;B=B*2){

        printf(" B= %d\n", B);
        clock_gettime(CLOCK_REALTIME, &startTime);
        dgemm_1B(array_A, array_B, D1, M_A.rows, M_B.cols, M_A.cols, B);
        clock_gettime(CLOCK_REALTIME, &endTime);
        calculateElapsedTimeAndGFLOPs(M_A.rows, M_B.cols, M_A.cols);  

 printf(" runTime = %f s \t %llu ns\tnumber of FloatingPoint Operations = %.0f\tGFLOPS = %f\n" ,runTime_Sec ,(long long unsigned int)runTime_Nsec ,numberOfOperations ,GFLOPS);                 


        printf("------------------------------|ijk algorithm done|---------------------------------------\n");

        clock_gettime(CLOCK_REALTIME, &startTime);
        dgemm_2B(array_A, array_B, D2, M_A.rows, M_B.cols, M_B.cols, B);
        clock_gettime(CLOCK_REALTIME, &endTime);
        calculateElapsedTimeAndGFLOPs(M_A.rows, M_B.cols, M_A.cols);

 printf(" runTime = %f s \t %llu ns\tnumber of FloatingPoint Operations = %.0f\tGFLOPS = %f\n" ,runTime_Sec ,(long long unsigned int)runTime_Nsec ,numberOfOperations ,GFLOPS);                 



        printf("------------------------------|jik algorithm done|---------------------------------------\n");

        clock_gettime(CLOCK_REALTIME, &startTime);
        dgemm_3B(array_A, array_B, D3, M_A.rows, M_B.cols, M_A.cols, B);
        clock_gettime(CLOCK_REALTIME, &endTime);
        calculateElapsedTimeAndGFLOPs(M_A.rows, M_B.cols, M_A.cols);



 printf(" runTime = %f s \t %llu ns\tnumber of FloatingPoint Operations = %.0f\tGFLOPS = %f\n" ,runTime_Sec ,(long long unsigned int)runTime_Nsec ,numberOfOperations ,GFLOPS);                 

        printf("------------------------------|kij algorithm done|---------------------------------------\n");


        clock_gettime(CLOCK_REALTIME, &startTime);
        dgemm_4B(array_A, array_B, D4, M_A.rows, M_B.cols, M_A.cols, B);
        clock_gettime(CLOCK_REALTIME, &endTime);
        calculateElapsedTimeAndGFLOPs(M_A.rows, M_B.cols, M_A.cols);

 printf(" runTime = %f s \t %llu ns\tnumber of FloatingPoint Operations = %.0f\tGFLOPS = %f\n" ,runTime_Sec ,(long long unsigned int)runTime_Nsec ,numberOfOperations ,GFLOPS);                 



        printf("------------------------------|ikj algorithm done|---------------------------------------\n");

        clock_gettime(CLOCK_REALTIME, &startTime);
        dgemm_5B(array_A, array_B, D5, M_A.rows, M_B.cols, M_A.cols, B);
        clock_gettime(CLOCK_REALTIME, &endTime);
        calculateElapsedTimeAndGFLOPs(M_A.rows, M_B.cols, M_A.cols);




 printf(" runTime = %f s \t %llu ns\tnumber of FloatingPoint Operations = %.0f\tGFLOPS = %f\n" ,runTime_Sec ,(long long unsigned int)runTime_Nsec ,numberOfOperations ,GFLOPS);                 


        printf("------------------------------|jki algorithm done|---------------------------------------\n");

        clock_gettime(CLOCK_REALTIME, &startTime);
        dgemm_6B(array_A, array_B, D6, M_A.rows, M_B.cols, M_B.rows, B);
        clock_gettime(CLOCK_REALTIME, &endTime);
        calculateElapsedTimeAndGFLOPs(M_A.rows, M_B.cols, M_A.cols);



       printf(" runTime = %f s \t %llu ns\tnumber of FloatingPoint Operations = %.0f\tGFLOPS = %f\n" ,runTime_Sec ,(long long unsigned int)runTime_Nsec ,numberOfOperations ,GFLOPS);                 
 


        printf("------------------------------|kji algorithm done|---------------------------------------\n");
}


/*  
     // *****************************************************************************
    //                                   Saving the result                        //
    //******************************************************************************

    ofstream ofile(argv[3], ios::binary);

    ofile.write((char*) &M_A.rows, sizeof(unsigned int));
    ofile.write((char*) &M_B.cols, sizeof(unsigned int));
    ofile.write((char*) D1, M_A.rows*M_B.cols*sizeof(float))	;
*/

    return 0;
}

