//header files included
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
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cublas_v2.h>

//declaring the tile width and height
//for tile based matrix multiplication
#define TILE_WIDTH 64
#define TILE_HEIGHT 64

//Namespace for std
using namespace std;

//structure declaration for storing rows and columns for a matrix
struct matrix{
    unsigned int rows;	//storing rows of a matrix
    unsigned int cols;	//storing columns of a matrix
};








//handle error alias name declaration
#define HANDLE_ERROR( err ) (HandleError( err, __FILE__, __LINE__ ))

//global kernal for matrix multiplication, takes in input matrices and sizes, and multiplies them
//matrix multiplication is being done tile by tile
__global__ void matrix_mult(float* array1, unsigned int rows1, unsigned int cols1, float* array2, unsigned int rows2, unsigned int cols2, float* array3)
{
    //shared memory takes one tile at a time
    __shared__ float S1[TILE_WIDTH][TILE_HEIGHT];	//to store tiles for array 1
    __shared__ float S2[TILE_HEIGHT][TILE_WIDTH];	//to store tiles for array 2

    //threads x and y index for the current block
    unsigned int tx=threadIdx.x;
    unsigned int ty=threadIdx.y;

    unsigned int c=blockIdx.x*blockDim.x + threadIdx.x;	//row value using x-index of current thread
    unsigned int r=blockIdx.y*blockDim.y + threadIdx.y;	//column value using y-index of current thread

    unsigned int idx=c*rows1+r;				//column major index, using row and column value

    float val=0;		//register to store multiplication result initialized to zero

    for(int m=0; m<1+((rows2-1)/TILE_WIDTH);m++)	//going over all tiles one by one, with each m
    {

        int var1=m*TILE_WIDTH+tx ;		//x thread value for current tile
        int var2=m*TILE_WIDTH+ty ;		//y thread value for current tile

        //copying a tile from array1
        if (r < rows1 && var1 < rows2)		//if the value is associated to a valid matrix coordinate in array1 then store it to shared memory S1
            S1[ty][tx]=array1[r + var1*rows1];//storing a "valid" value from array to shared memory
        else
            S1[ty][tx]=0;					//storing zero, since there is no valid value
        __syncthreads();						//syncing all threads once shared memory S1 is stored

        //copying a tile from array2
        if(c < cols2 && var2 < rows2)	//if value is associates to a valid matrix coordinate in array2 then store it to shared memory S2
            S2[ty][tx]=array2[var2+rows2*c];	//storing the valid value
        else
            S2[ty][tx]=0;		//storing zero, since no valid value
        __syncthreads();		//synchronizing threads


        for(int i=0; i<TILE_WIDTH;i++)	//going over entire tile, ty row in S1 and tx column in S2
            val+=S1[ty][i]*S2[i][tx];	//and multiplying elements
        __syncthreads();		//synchronizing threads

    }

    if(r < rows1 && c< cols2)	//removing degenerate cases
        array3[idx]=val;	//saving multiplication result to global memory

}

int main(int argc, char* argv[])
{
    if(argc != 4) //there should be four arguments, Usage: prog matrix1.mtx matrix2.mtx matrix3.mtx
        return 1; //exit and return an error

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

    infile_B.close();

    if(M_A.cols!=M_B.rows) //checking if the two matrices can be multiplied
    {
        cout<<"Illegal matrix sizes: "<<M_A.cols<<" != "<<M_B.rows<<endl;
        return 1;
    }
    // *****************************************************************************
    //                                   allocate to the host                      //
    //******************************************************************************
    float* array_C=(float*)malloc(M_A.rows*M_B.cols*sizeof(float));//array to store gpu result in column major format
    //GPU DEVICE PROPERTIES and selecting a GPU for calculation
    int nDevices;
    cudaGetDeviceCount(&nDevices);

    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, 0);	//using GPU0

    //BLOCK AND GRID SIZE DECLARATION
    float thread_block=sqrt(prop.maxThreadsPerBlock);	//2D blocks used
    dim3 DimGrid(ceil(M_B.cols/thread_block),ceil(M_A.rows/thread_block),1); //image saved as a 2D grid
    dim3 DimBlock(thread_block,thread_block,1);

    size_t Sbytes = 2* DimBlock.x * DimBlock.y ;	//2 arrays used in the calculation, hence 2 * DimBlock.x * DimBlock.y

    //Checking if sufficient shared memory available or not

    if(prop.sharedMemPerBlock < Sbytes){
        std::cout<<"ERROR: insufficient shared memory"<<std::endl;
        exit(1);
    }

    // *****************************************************************************
    //                                   allocate to the GPU                      //
    //******************************************************************************

    float *array_A_gpu, *array_B_gpu, *array_C_gpu;//gpu arrays declared

    cudaMalloc(&array_A_gpu,M_A.rows*M_A.cols*sizeof(float)); //allocate space to store arrayA

    cudaMalloc(&array_B_gpu,M_B.rows*M_B.cols*sizeof(float)); //allocate space to store arrayB

    cudaMalloc(&array_C_gpu,M_A.rows*M_B.cols*sizeof(float)); //allocate space to store gpu result

    //COPY TO GPU MEMORY
    cudaMemcpy(array_A_gpu, array_A, M_A.rows*M_A.cols*sizeof(float), cudaMemcpyHostToDevice);//copy arrayA to gpu

    cudaMemcpy(array_B_gpu, array_B, M_B.rows*M_B.cols*sizeof(float), cudaMemcpyHostToDevice);//copy arrayB to gpu

    cudaMemcpy(array_C_gpu, array_C, M_A.rows*M_B.cols*sizeof(float), cudaMemcpyHostToDevice);//copy arrayC to gpu

    // *****************************************************************************
    //                                   allocate to the GPU                      //
    //******************************************************************************

    //time measurement for matrix multiplication
    cudaEvent_t start1, stop1;

    cudaEventCreate(&start1);
    cudaEventCreate(&stop1);

    //MATRIX MULTIPLICATION USING KERNEL
    cudaEventRecord(start1);
    matrix_mult<<<DimGrid, DimBlock, Sbytes>>>(array_A_gpu,M_A.rows,M_A.cols,array_B_gpu,M_B.rows,M_B.cols,array_C_gpu);//calling the kernel
    cudaEventRecord(stop1);

    cudaEventSynchronize(stop1);

    float milliseconds1 = 0;//storing the execution time in milliseconds

    cudaEventElapsedTime(&milliseconds1, start1, stop1);//get the time in milliseconds
    float msecPerMatrixMul = milliseconds1;
    double flopsPerMatrixMul = 2.0 * (double) M_A.rows *(double) M_B.cols *(double) M_A.cols;
    double gigaFlops = (flopsPerMatrixMul * 1.0e-9f) / (msecPerMatrixMul / 1000.0f);
    printf("Performance= %.2f GFlop/s, Time= %.3f msec",gigaFlops,msecPerMatrixMul);

    //copy to CPU MEMORY
    cudaMemcpy(array_C, array_C_gpu, M_A.rows*M_B.cols*sizeof(float), cudaMemcpyDeviceToHost);//copying result of multiplication from gpu to cpu


    // *****************************************************************************
    //                                   Saving the result                        //
    //******************************************************************************
    //SAVING THE OUTPUT MATRIX
    ofstream ofile(argv[3], ios::binary);

    ofile.write((char*) &M_A.rows, sizeof(unsigned int));//writing the rows
    ofile.write((char*) &M_B.cols, sizeof(unsigned int));//writing the cols
    ofile.write((char*) array_C , M_A.rows*M_B.cols*sizeof(float));//writing all elements


    return 0;
}
