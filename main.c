#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>



//FROM BOOK START

//define zdefiniowane na stronie 120 w ksiazce
#define BLOCK_LOW(id,p,n) ((id)*(n)/(p))
#define BLOCK_HIGH(id,p,n) (BLOCK_LOW((id)+1,p,n)-1)
#define BLOCK_SIZE(id,p,n) (BLOCK_HIGH(id,p,n)-BLOCK_LOW(id,p,n)+1)

#define PTR_SIZE (sizeof(void*))

#define mpitype MPI_DOUBLE

//funkcja na stronie 487 w ksiazce
void *my_malloc(
    int id, //process rank
    int bytes) //bytes to allocate
{
    void *buffer;
    //printf("bytes: %d", bytes);
    if((buffer=malloc((size_t)bytes))==NULL){
        printf("Error: Malloc failed for process %d\n", id);
        fflush(stdout);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    return buffer;
}
      
//funkcja na stronie 487 w ksiazce      
int get_size(MPI_Datatype t){
    if(t==MPI_BYTE) return sizeof(char);
    if(t == MPI_DOUBLE) return sizeof(double);
    if(t == MPI_FLOAT) return sizeof(float);
    if(t == MPI_INT) return sizeof(int);
    printf ("Error: Unrecognized argwrlent to ' get_size' \n") ;
    fflush (stdout);
    MPI_Abort (MPI_COMM_WORLD, 1);
}

//funkcja na stronie 495 w ksiazce
void read_row_striped_matrix(
    char *s, // IN file name
    void ***subs, // OUT 2d submatrix indicex
    void **storage, // OUT submatrix stored 
    MPI_Datatype dtype, //matrix element type
    int *m, //matrix rows
    int *n, //matrix cols
    MPI_Comm comm)
{
    int datum_size; // size of matrix element
    int i;
    int id; //process rank
    FILE *infileptr; //input file pointer
    int local_rows; // rows on this proc
    void ** lptr; // pointer into 'subs'
    int p; // number of processes
    void *rptr; //pointer into 'storage'
    MPI_Status status; //result of receive
    double x; //result of read
    
    MPI_Comm_size(comm, &p);
    MPI_Comm_rank(comm, &id);
    datum_size = get_size(dtype);
    
    /* Process p~'l opens file, reads size of matrix,
    and broadcasts matrix dimensions to other procs ~!*/
    if (id == (p-1)){
        infileptr = fopen(s, "r");
        if(infileptr == NULL) *m = 0;
        else{
            //fread (m, sizeof(int), 1, infileptr);
            fscanf(infileptr, "%d", m);
            //fread (n, sizeof(int), 1, infileptr);
            fscanf(infileptr, "%d", n);
        }
    }
    MPI_Bcast(m, 1, MPI_INT, p-1, comm);
    
    if(!(*m)) MPI_Abort (MPI_COMM_WORLD, 2);
    
    MPI_Bcast (n, 1, MPI_INT, p-1 ,comm);
    
    local_rows = BLOCK_SIZE(id, p, *m);
    
    /*(1* r:ynamically allocate matrix. Allo"w double subscripting
    through' a'. */
    //printf("l rows: %d", local_rows);
    
    *storage = (void*)my_malloc(id,
        local_rows * *n * datum_size);
    *subs = (void**) my_malloc(id, local_rows * PTR_SIZE);
    
    lptr = (void*) &(*subs[0]);
    rptr = (void *) *storage;
    for(i = 0; i< local_rows; i++){
        *(lptr++) = (void*) rptr;
        rptr += *n*datum_size;
    }
    
    /* Process p-l reads blocks ot rows from file and
    sends each block to the correct destination proces:>.
    The last block it keeps. */

    if(id == (p-1)){
        for (i = 0; i <= p-1; i++) {				
            double *temp = *storage;
            int j;
            for(j = 0; j < BLOCK_SIZE(i,p,*m) * *n; j++)		
            {
                fscanf(infileptr, "%lf", &x);
                temp[j] = x;
            }
            if(i != p-1){
                MPI_Send (*storage, BLOCK_SIZE(i,p,*m) * *n, dtype, i, 0, comm);      
            }            
        }
        fclose (infileptr);
    }else{
        MPI_Recv (*storage, local_rows* *n, dtype, p-1,
            0, comm, &status);
    }
}

// FROM BOOK END

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

int size = 3;

int getRowById(int id){
    return id/size;
}

int getColById(int id){
    return id%size;
}

int getIdToSendLeft(int id){
    if(id%size == 0)
        return id+size-1;
    else
        return id-1;
}

int getIdToRecvRight(int id){
    if((id+1)%size == 0)
        return id-size+1;
    else
        return id+1;
}

int getIdOfLastInRow(int id){
    int i = id;
    while((i+1)%size != 0)
        i++;
    return i;
}

int getIdToSendUp(int id){
    if(id-size < 0)
        return id+(size*size - size);
    else
        return id-size;
}

int getIdToRecvDown(int id){
    if(id+size > size*size - 1)
        return id%size;
    else
        return id+size;
}

int getIdLeftByDist(int id, int dist){
    int leftId = id;
    int i;
    for(i = 0; i  <dist; i++)
        leftId = getIdToSendLeft(leftId);
    return leftId;
}

int getIdRightByDist(int id, int dist){
    int rightId = id;
    int i;
    for(i = 0; i  <dist; i++)
        rightId = getIdToRecvRight(rightId);
    return rightId;
}

int getIdUpByDist(int id, int dist){
    int upId = id;
    int i;
    for(i = 0; i  <dist; i++)
        upId = getIdToSendUp(upId);
    return upId;
}

int getIdDownByDist(int id, int dist){
    int downId = id;
    int i;
    for(i = 0; i  <dist; i++)
        downId = getIdToRecvDown(downId);
    return downId;
}
    
int main(int argc, char *argv[])
{      
    int    myid, numprocs;
    int    namelen;
    char   processor_name[MPI_MAX_PROCESSOR_NAME];
	int i,j,l;
	struct MultiplicationHelper AxBElements[3][3];
    
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD,&myid);
    MPI_Get_processor_name(processor_name,&namelen);
    
    //READ FILE
    double **a; //a matrix
    double *storageA; //matrix elements stored here
    double **b; //b matrix
    double *storageB; //matrix elements stored here
    int m; //rows in matrix
    int n; //columns in matrix
    
    //pobieranie macierzy A z pliku
    //dane sa pobierane i zapisywane w ostatnim elemencie wiersza
    read_row_striped_matrix(argv[1], (void*) &a,
        (void*) &storageA, mpitype, &m, &n, MPI_COMM_WORLD);
        
    //pobieranie macierzy B z pliku
    //dane sa pobierane i zapisywane w ostatnim elemencie wiersza
    read_row_striped_matrix(argv[2], (void*) &b,
        (void*) &storageB, mpitype, &m, &n, MPI_COMM_WORLD);
        
    int myI = getRowById(myid);
    int myJ = getColById(myid);
    
    int myA;
    int myB;
        
        
    //TODO
    //KAZDY PROCES KTORY JEST 'OSTATNI W DANYM WIERSZU' ZAWIERA WSZYSTKIE ELEMENTY MACIERZY Z WIERSZA
    //NALEZY ZROBIC DLA DOWOLNEJ WIELKOSCI MACIERZY ZEBY WYSYLAL WARTOSCI
    if((myid+1)%size==0){
        
        //TODO na double
        myA = storageA[2];
        int tosendA1 = storageA[0];
        int tosendA2 = storageA[1];
        

        //TODO
        //W PETLI WYSYLANIE WARTOSCI DO WSZYSTKICH ELEMENTOW WIERSZA
        MPI_Send(&tosendA1, 1, MPI_INT, getIdToSendLeft(myid)-1, 0, MPI_COMM_WORLD);
        MPI_Send(&tosendA2, 1, MPI_INT, getIdToSendLeft(myid), 0, MPI_COMM_WORLD);
        myA = storageA[2];
    }else {
        MPI_Recv(&myA, 1, MPI_INT, getIdOfLastInRow(myid), 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    
    if((myid+1)%size==0){
        myB = storageB[2];
        int tosendB1 = storageB[0];
        int tosendB2 = storageB[1];
        
        //TODO
        //W PETLI WYSYLANIE WARTOSCI DO WSZYSTKICH ELEMENTOW WIERSZA
        MPI_Send(&tosendB1, 1, MPI_INT, getIdToSendLeft(myid)-1, 0, MPI_COMM_WORLD);
        MPI_Send(&tosendB2, 1, MPI_INT, getIdToSendLeft(myid), 0, MPI_COMM_WORLD);
        myB = storageB[2];
    }else {
        MPI_Recv(&myB, 1, MPI_INT, getIdOfLastInRow(myid), 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    
    
    //TUTAJ NASTEPUJE PRZESUWANIE W LEWO ELEMENTOW MACIERZY A O ILOSC ROWNA NUMEROWI WIERSZA (NUMEROWANIE OD 0)
    int distToSendHorizontally = getRowById(myid);
    if(distToSendHorizontally!=0){
        MPI_Send(&myA, 1, MPI_INT, getIdLeftByDist(myid, distToSendHorizontally), 0, MPI_COMM_WORLD);
        MPI_Recv(&myA, 1, MPI_INT, getIdRightByDist(myid, distToSendHorizontally), 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    
    //TUTAJ NASTEPUJE PRZESUWANIE W GORE ELEMENTOW MACIERZY B O ILOSC ROWNA NUMEROWI KOLUMNY (NUMEROWANIE OD 0)
    int distToSendVertically = getColById(myid);
    if(distToSendVertically!=0){
        MPI_Send(&myB, 1, MPI_INT, getIdUpByDist(myid, distToSendVertically), 0, MPI_COMM_WORLD);
        MPI_Recv(&myB, 1, MPI_INT, getIdDownByDist(myid, distToSendVertically), 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    
    
    //FRAGMENT STAREGO KODU. NALEZY ODKOMENTOWAC ABY UZYC ZAHARDKODOWANYCH MACIERZY 3X3 Z TEGO KODU
    /*for(i = 0; i < 3; ++i){
        for(j = 0; j < 3; ++j){
            AxBElements[i][j].A = matrixA[i][(j+i) % (3)];
            AxBElements[i][j].B = matrixB[(i+j) %3][j];
        }
    }
        
    for(i = 0; i < 3; ++i){
		for(j = 0; j < 3; ++j){    
            matrixA[i][j] = AxBElements[i][j].A;
            matrixB[i][j] = AxBElements[i][j].B;
        }
	}
    
     myA = matrixA[myI][myJ];
     myB = matrixB[myI][myJ];
     */
     
    int myRes = 0;
    
    //ROZPOCZYNA SIE ITERACJA OBLICZANIA
    //W KAZDYM KROKU LICZONY JEST ILOCZYN OBECNEJ WARTOSCI A ORAZ B
    //A NASTEPNIE WARTOSCI MACIERZY A SA PRZESYLANE W LEWO, A WARTOSCI MACIERZY B W GORE
    myRes += myA*myB;
    for(l = 1; l < 3; ++l){
        MPI_Send(&myA, 1, MPI_INT, getIdToSendLeft(myid), 0, MPI_COMM_WORLD);
        MPI_Send(&myB, 1, MPI_INT, getIdToSendUp(myid), 0, MPI_COMM_WORLD);
        MPI_Recv(&myA, 1, MPI_INT, getIdToRecvRight(myid), 0, MPI_COMM_WORLD,
             MPI_STATUS_IGNORE);
        MPI_Recv(&myB, 1, MPI_INT, getIdToRecvDown(myid), 0, MPI_COMM_WORLD,
             MPI_STATUS_IGNORE);
        myRes += myA*myB;
    }
    
    //NASTEPUJE ZEBRANIE WSZYSTKICH WYNIKOW DO PROCESU 0       
    if(myid==0){
        resultMatrix[0][0] = myRes;
        for(i = 1; i < numprocs; i++){
            int res;
            MPI_Recv(&res, 1, MPI_INT, i, 0, MPI_COMM_WORLD,
                MPI_STATUS_IGNORE);
            resultMatrix[getRowById(i)][getColById(i)] = res;
        }
    }else{
        MPI_Send(&myRes, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
    }
	
    //WYPISANIE WYNIKU
    if(myid == 0){
        printf("result\n");	
        for(i = 0; i < 3; ++i){
            for(j = 0; j < 3; ++j){
                printf("%d\t", resultMatrix[i][j]);
            }
            printf("\n");
        }
    }

    

    MPI_Finalize();

    return 0;
}

