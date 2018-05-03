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
#define MPI_TYPE MPI_FLOAT
#define TYPE float
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

//zeby przechowywac wymiary macierzy razem
struct MatrixDim{
	int row;
	int col;
};
	
int getRowById(int id, struct MatrixDim m){
    return id/m.col;
}

int getColById(int id, struct MatrixDim m){
    return id%m.row;
}

int getIdToSendLeft(int id, struct MatrixDim m){
    if(id%m.col == 0)
        return id+m.col-1;
    else
        return id-1;
}

int getIdToRecvRight(int id, struct MatrixDim m){
    if((id+1)%m.col == 0)
        return id-m.col+1;
    else
        return id+1;
}

int getIdOfLastInRow(int id, struct MatrixDim m){
    int i = id;
    while((i+1)%m.col != 0)
        i++;
    return i;
}

int getIdToSendUp(int id, struct MatrixDim m){
    if(id-m.row < 0)
        return id+(m.row*m.col - m.col);
    else
        return id-m.col;
}

int getIdToRecvDown(int id, struct MatrixDim m){
    if(id+m.row > m.row*m.col - 1)
        return id%m.col;
    else
        return id+m.col;
}

int getIdLeftByDist(int id, int dist, struct MatrixDim m){
    int leftId = id;
    int i;
    for(i = 0; i  <dist; i++)
        leftId = getIdToSendLeft(leftId,m);
    return leftId;
}

int getIdRightByDist(int id, int dist, struct MatrixDim m){
    int rightId = id;
    int i;
    for(i = 0; i  <dist; i++)
        rightId = getIdToRecvRight(rightId,m);
    return rightId;
}

int getIdUpByDist(int id, int dist, struct MatrixDim m){
    int upId = id;
    int i;
    for(i = 0; i  <dist; i++)
        upId = getIdToSendUp(upId,m);
    return upId;
}

int getIdDownByDist(int id, int dist, struct MatrixDim m){
    int downId = id;
    int i;
    for(i = 0; i  <dist; i++)
        downId = getIdToRecvDown(downId,m);
    return downId;
}

int main(int argc, char *argv[])
{      
	int    myid, numprocs;
    int    namelen;
    char   processor_name[MPI_MAX_PROCESSOR_NAME];
	int i,j,l;
    
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD,&myid);
    MPI_Get_processor_name(processor_name,&namelen);
    //READ FILE
    double **a; //a matrix
    double *storageA; //matrix elements stored here
    double **b; //b matrix
    double *storageB; //matrix elements stored here
    struct MatrixDim A, B, RES;
    
    //pobieranie macierzy A z pliku
    //dane sa pobierane i zapisywane w ostatnim elemencie wiersza
    read_row_striped_matrix(argv[1], (void*) &a,
        (void*) &storageA, mpitype, &A.row, &A.col, MPI_COMM_WORLD);

    //pobieranie macierzy B z pliku
    //dane sa pobierane i zapisywane w ostatnim elemencie wiersza
    read_row_striped_matrix(argv[2], (void*) &b,
        (void*) &storageB, mpitype, &B.row, &B.col, MPI_COMM_WORLD);
        
    RES.row = A.row;
    RES.col = B.col;
    //TODO
    // jeżeli n != n2 MPI ABORD, bo nie mozna mnozyc
    if(myid == 0 && A.col != B.row){
    	MPI_Abort (MPI_COMM_WORLD, 1);
    }
	if(myid == 0){
		printf("%d\t%d\t%d\t%d\t%d\t%d\n",A.col, A.row, B.col, B.row, RES.col, RES.row);
	}
    TYPE myA;
    TYPE myB;
    //KAZDY PROCES KTORY JEST 'OSTATNI W DANYM WIERSZU' ZAWIERA WSZYSTKIE ELEMENTY MACIERZY Z WIERSZA
    if((myid+1)%A.col==0 && myid < A.row*A.col){
        //W PETLI WYSYLANIE WARTOSCI DO WSZYSTKICH ELEMENTOW WIERSZA
        for(i = 0; i < A.col-1; ++i) {
        	TYPE tosend = storageA[A.col-2 -i];
        	MPI_Send(&tosend, 1, MPI_TYPE, getIdToSendLeft(myid,A)-i, 0, MPI_COMM_WORLD);
        }
        myA = storageA[A.col-1];
    }else 
    //Jesli beda macierze innych wymiarow to recv nie zablokuje programu bo bedzietyle oczekiwan ile elementow macierzy A, to samo tyczy sie nizej B
    if(myid < A.col*A.row && (myid+1)%A.col != 0){
        MPI_Recv(&myA, 1, MPI_TYPE, getIdOfLastInRow(myid,A), 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    
    if((myid+1)%B.col==0 && myid < B.col*B.row){
        //W PETLI WYSYLANIE WARTOSCI DO WSZYSTKICH ELEMENTOW WIERSZA
        for(i = 0; i < B.col-1; ++i) {
        	TYPE tosend = storageB[B.col-2 -i];
        	MPI_Send(&tosend, 1, MPI_TYPE, getIdToSendLeft(myid,B)-i, 0, MPI_COMM_WORLD);
        }
        myB = storageB[B.col-1];
    }else if(myid < B.row*B.col && (myid+1)%B.col !=0){
        MPI_Recv(&myB, 1, MPI_TYPE, getIdOfLastInRow(myid,B), 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    if(myid == 0) printf("thyju");
    //TUTAJ NASTEPUJE PRZESUWANIE W LEWO ELEMENTOW MACIERZY A O ILOSC ROWNA NUMEROWI WIERSZA (NUMEROWANIE OD 0)
    int distToSendHorizontally = getRowById(myid,A);
    if(distToSendHorizontally!=0){
        MPI_Send(&myA, 1, MPI_TYPE, getIdLeftByDist(myid, distToSendHorizontally,A), 0, MPI_COMM_WORLD);
        MPI_Recv(&myA, 1, MPI_TYPE, getIdRightByDist(myid, distToSendHorizontally,A), 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    
    //TUTAJ NASTEPUJE PRZESUWANIE W GORE ELEMENTOW MACIERZY B O ILOSC ROWNA NUMEROWI KOLUMNY (NUMEROWANIE OD 0)
    int distToSendVertically = getColById(myid,B);
    if(distToSendVertically!=0){
        MPI_Send(&myB, 1, MPI_TYPE, getIdUpByDist(myid, distToSendVertically,B), 0, MPI_COMM_WORLD);
        MPI_Recv(&myB, 1, MPI_TYPE, getIdDownByDist(myid, distToSendVertically,B), 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
     
    TYPE myRes = 0;
    
    //ROZPOCZYNA SIE ITERACJA OBLICZANIA
    //W KAZDYM KROKU LICZONY JEST ILOCZYN OBECNEJ WARTOSCI A ORAZ B
    //A NASTEPNIE WARTOSCI MACIERZY A SA PRZESYLANE W LEWO, A WARTOSCI MACIERZY B W GORE
    myRes += myA*myB;
    for(l = 1; l < A.col; ++l){
        MPI_Send(&myA, 1, MPI_TYPE, getIdToSendLeft(myid,A), 0, MPI_COMM_WORLD);
        MPI_Send(&myB, 1, MPI_TYPE, getIdToSendUp(myid,B), 0, MPI_COMM_WORLD);
        MPI_Recv(&myA, 1, MPI_TYPE, getIdToRecvRight(myid,A), 0, MPI_COMM_WORLD,
             MPI_STATUS_IGNORE);
        MPI_Recv(&myB, 1, MPI_TYPE, getIdToRecvDown(myid,B), 0, MPI_COMM_WORLD,
             MPI_STATUS_IGNORE);
        myRes += myA*myB;
    }
    
	//przeniesc do MPI_INIT??
	// i zwalniać w MPI_FINALIZE?
	//TODO nie zwalniane jest
    TYPE **resultMatrix = (TYPE **)malloc(RES.row * sizeof(TYPE *));
    for (i=0; i<RES.row; i++)
         resultMatrix[i] = (TYPE *)malloc(RES.col * sizeof(TYPE));
    
    //NASTEPUJE ZEBRANIE WSZYSTKICH WYNIKOW DO PROCESU 0       
    int size = RES.col*RES.row;
    if(size > numprocs) size = numprocs;
        
    if(myid==0){
        resultMatrix[0][0] = myRes;
        for(i = 1; i < size; i++){
            TYPE res;
            MPI_Recv(&res, 1, MPI_TYPE, i, 0, MPI_COMM_WORLD,
                MPI_STATUS_IGNORE);
            //printf("r: %d\tc: %d\tw: %.2f\n",getRowById(i,RES),getColById(i,RES), res);
            resultMatrix[getRowById(i,RES)][getColById(i,RES)] = res;
        }
    }else if(myid < size){
    	
        MPI_Send(&myRes, 1, MPI_TYPE, 0, 0, MPI_COMM_WORLD);
    }
	
    //WYPISANIE WYNIKU
    if(myid == 0){
        printf("result\n");	
        for(i = 0; i < RES.row; ++i){
            for(j = 0; j < RES.col; ++j){
                printf("%.2f\t", resultMatrix[i][j]);
            }
            printf("\n");
        }
    }
    

    MPI_Finalize();
    
    return 0;
}

