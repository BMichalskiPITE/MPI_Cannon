#include "mpi.h"
#include <stdio.h>



//na chwile obecna proponuje zhardcodowac maciez 3x3,
//jak algorytm zadziala to zrobic uniklanie nxm x mxl
int matrixA[3][3] = {{1,2,1},{2,5,3},{2,4,1}};
int matrixB[3][3] = {{2,1,0},{1,1,1},{3,2,2}};
//wynik AxB = {{7,5,4},{18,13,11},{11,8,6}}
int resultMatrix[3][3];
struct MultiplicationHelper {
	int A;
	int B;
};

int main(int argc, char *argv[])
{
    int    myid, numprocs;
    int    namelen;
    char   processor_name[MPI_MAX_PROCESSOR_NAME];
	int i,j,l;
	struct MultiplicationHelper AxBElements[3][3];
	for(l = 0; l < 3; ++l){
		printf("\n%d\n",l);
		for(i = 0; i < 3; ++i){
			for(j = 0; j < 3; ++j){
				AxBElements[i][j].A = matrixA[i][(i+j+l) % (3)];
				AxBElements[i][j].B = matrixB[(i+j+l) %3][j];
				printf("%d|%d\t", AxBElements[i][j].A,AxBElements[i][j].B);
				resultMatrix[i][j] +=AxBElements[i][j].A * AxBElements[i][j].B;
			}
			printf("\n");
		}
	}
	
	printf("result\n");	
	for(i = 0; i < 3; ++i){
		for(j = 0; j < 3; ++j){
			printf("%d\t", resultMatrix[i][j]);
		}
		printf("\n");
	}
		
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD,&myid);
    MPI_Get_processor_name(processor_name,&namelen);

    

    MPI_Finalize();

    return 0;
}

