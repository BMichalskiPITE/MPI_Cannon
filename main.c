#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "mpe.h"

#define START_BCAST 0
#define END_BCAST 1
#define START_ALLRED 2
#define END_ALLRED 3
#define START_RECV 4
#define END_RECV 5
#define START_SEND 6
#define END_SEND 7


int MPI_Init( int *argc, char ***argv )
{

    int ret, myid;
    ret = PMPI_Init(argc, argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
	MPE_Init_log();
    if(myid == 0){
		MPE_Describe_state(START_BCAST, END_BCAST,"broadcast","red");
		MPE_Describe_state(START_ALLRED, END_ALLRED,"all reduce","yellow");
		MPE_Describe_state(START_RECV, END_RECV,"receive","blue");
		MPE_Describe_state(START_SEND, END_SEND,"send","white");
	}
	MPE_Start_log();
    return ret;

}

int MPI_Finalize( void )
{

    int ret;
    MPE_Finish_log("mpe_logs");
    ret = PMPI_Finalize();

    return ret;

}

int MPI_Bcast( void *buf, int count, MPI_Datatype datatype,  
	       int root, MPI_Comm comm ) 
{ 
    int ret, myid;
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
    MPE_Log_event(START_BCAST,myid, "Start\tbroadcast");
    ret = PMPI_Bcast( buf, count, datatype, root, comm ); 
    MPE_Log_event(END_BCAST,myid, "End\tbroadcast");
    return ret; 
} 

int MPI_Send(const void *buf, int count, MPI_Datatype datatype, int dest, int tag,
                    MPI_Comm comm)
{ 
    int ret, myid;
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
    MPE_Log_event(START_SEND,myid, "Start\tsend");
    ret = PMPI_Send(buf, count, datatype, dest, tag, comm );
    MPE_Log_event(END_SEND,myid, "End\tsend");
    return ret; 
} 

int MPI_Recv(void *buf, int count, MPI_Datatype datatype, int source, int tag,
                           MPI_Comm comm, MPI_Status *status)                           
{ 
    int ret, myid;
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
    MPE_Log_event(START_RECV,myid, "Start\treceive");
    ret = PMPI_Recv(buf, count, datatype, source, tag, comm, status);
    MPE_Log_event(END_RECV,myid, "End\treceive");
    return ret; 
} 


int MPI_Allreduce(const void *sendbuf, void *recvbuf, int count,
                         MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
{ 
    int ret, myid;
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
    MPE_Log_event(START_ALLRED,myid, "Start\tallreduce");
    ret = PMPI_Allreduce(sendbuf, recvbuf, count, datatype, op, comm);
    MPE_Log_event(END_ALLRED,myid, "End\tallreduce");
    return ret; 
} 

//FROM BOOK START

//define zdefiniowane na stronie 120 w ksiazce
#define BLOCK_LOW(id,p,n) ((id)*(n)/(p))
#define BLOCK_HIGH(id,p,n) (BLOCK_LOW((id)+1,p,n)-1)
#define BLOCK_SIZE(id,p,n) (BLOCK_HIGH(id,p,n)-BLOCK_LOW(id,p,n)+1)

#define PTR_SIZE (sizeof(void*))

#define mpitype MPI_DOUBLE
#define MPI_TYPE MPI_DOUBLE
#define TYPE double
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

double **alloc_2d(int rows, int cols) {
    double *data = (double *)malloc(rows*cols*sizeof(double));
    double **array= (double **)malloc(rows*sizeof(double*));
    for (int i=0; i<rows; i++)
        array[i] = &(data[cols*i]);

    return array;
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

int getProcLeft(int proc, int procsDim){
    if(proc%procsDim == 0)
        return proc+procsDim-1;
    else
        return proc-1;
}

int getProcLeftByDist(int proc, int pos, int dist, int dim, int procsDim){
    int colsPerProc = dim/procsDim;
    int toLeftBy = dist - pos;
    if(toLeftBy <=0)
        return proc;
    
    toLeftBy = 1+(toLeftBy-1)/colsPerProc;
    
    int i;
    int leftProc = proc;
    for(i = 0; i  < toLeftBy; i++)
        leftProc = getProcLeft(leftProc,procsDim);
    return leftProc;
}

int getProcTop(int proc, int procsDim){
    proc = proc - procsDim;
    if(proc >=0)
        return proc;
    else
        return procsDim*procsDim-fabs(proc);
   /* if(proc%procsDim == 0)
        return proc+procsDim-1;
    else
        return proc-1;*/
}

int getProcTopByDist(int proc, int pos, int dist, int dim, int procsDim){
    int rowsPerProc = dim/procsDim;
    int toTopBy = dist - pos;
    if(toTopBy <=0)
        return proc;
    
    toTopBy = 1+(toTopBy-1)/rowsPerProc;
    
    int i;
    int topProc = proc;
    for(i = 0; i  < toTopBy; i++)
        topProc = getProcTop(topProc,procsDim);
    return topProc;
}

int getProcDown(int proc, int procsDim){
    proc = (proc+procsDim);
    if(proc <= procsDim*procsDim - 1)
        return proc;
    else
        return (proc)%(procsDim*procsDim);
}

int getProcDownByDist(int proc, int pos, int dist, int dim, int procsDim){

    int rowsPerProc = dim/procsDim;
    int toDownBy = dist -(rowsPerProc - pos - 1); 
    if(toDownBy <=0)
        return proc;
    
    toDownBy = 1+(toDownBy-1)/rowsPerProc;
    
    int i;
    int downProc = proc;
    for(i = 0; i  < toDownBy; i++)
        downProc = getProcDown(downProc,procsDim);
    return downProc;
}

int getNewLeftPos(int pos, int dist, int colsPerProc){
    pos = (pos-dist);
    if(pos>=0)
        return pos;
    int i = 0;
    while( pos < 0)
        pos += colsPerProc;
    
    return ((pos)%colsPerProc);
}

int getProcRight(int proc, int procsDim){
    if((proc+1)%procsDim == 0)
        return proc-procsDim+1;
    else
        return proc+1;
}

int getProcRightByDist(int proc, int pos, int dist, int dim, int procsDim){
    int colsPerProc = dim/procsDim;
    int toRightBy = dist - (colsPerProc - pos -1);
    if(toRightBy <=0)
        return proc;
    
    toRightBy = 1+(toRightBy-1)/colsPerProc;
    
    int i;
    int rightProc = proc;
    for(i = 0; i  < toRightBy; i++)
        rightProc = getProcRight(rightProc,procsDim);
    return rightProc;
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
    
    //PRZESYLANIE PO WCZYTANIU START!
    int procSqrt = sqrt(numprocs); 
    int myRows = A.row / procSqrt; //TAKI JEST WYMIAR POJEDYNCZEGO BLOKU!!
    //MUSIMY TERAZ ZAALOKOWAC!!

    double **mA = alloc_2d(myRows, myRows);
     
    double **nA = alloc_2d(myRows, myRows);
     
    double **mB = alloc_2d(myRows, myRows);
     
    double **nB = alloc_2d(myRows, myRows);
     
    double **res = alloc_2d(myRows, myRows);
	
	for(i = 0; i<myRows; i++)
    {
        for(j = 0; j<myRows; j++)
        {
            res[i][j] = 0;
        }
    }
    
    int matrixDim = A.col;

    
    int proc = 0;
    int meProc = (myid%A.col);
    //ITEROWANIE PO PROCESACH W DANYM WIERSZU BLOKU
    int rowOfProcess = myid/procSqrt; // jestem w danym wierszu procesow
    int rowsPerProc = matrixDim/procSqrt;
    int colsPerProc = rowsPerProc;
    int currentRow = 0;
    int currentRowForSend = 0;
    int colOffset = (myid%procSqrt)*myRows; //offset kolumny, czyli pierwszy proces ma 0, pozniej ile trzeba dodac do kolumny zeby miec pierwsza kolumne przynalezna do procesu
    
    
    //WYSYLANE MACIERZY A DO PROCESOW
        currentRowForSend = 0;
        //TERAZ TRZEBA WYSLAC WSZYSTKIE POSIADANE WIERSZE DO WSZYSTKICH PROCESOW!
        int rows = BLOCK_SIZE(myid, numprocs, A.col); // ILE PRZECZYTALEM WIERSZY MACIERZY
        i = 0;
        int row = 0; 
        for(row = 0; row<rows; row++)
        {
            for(i = 0; i<A.col; i+=myRows)
            {
                int procInRow = i/myRows; //NUMER PROCESU W WIERSZU DO WYSLANIA
                if(procInRow == meProc){
                    continue;
                }else{
                    int procNumber = rowOfProcess*procSqrt+procInRow;
                    MPI_Send(&a[currentRowForSend][i], myRows, MPI_TYPE, procNumber, row, MPI_COMM_WORLD);
                }
                
            }  
            currentRowForSend++;
        }
    
    
	//ODBIERANIE MACIERZY A OD PROCESOW
    for(proc = 0; proc < procSqrt; proc ++)
    {
        if(meProc == proc)
        { //JESLI MAM JUZ TE ELEMENTY CO POTRZEBUJE, TO NIE TRZEBA WYSYLAC, BIORE OD SIEBIE OD PRZECZYTANEJ MACIERZY DO mA
            int rows = BLOCK_SIZE(myid, numprocs, A.col); // ILE PRZECZYTALEM WIERSZY MACIERZY
            int i = 0;
            int row = 0;
            for(row = 0; row<rows; row++)
            {
                for(i = 0; i<myRows; i++)
                {
                    mA[currentRow][i] = a[row][colOffset+i];
                }   
                currentRow++;
            }
        }
        else
        {
            //TRZEBA POBRAC WARTOSCI MACIERZY Z PROCESU O ID rowOfProcess*procSqrt+proc!!
            int procNumber = rowOfProcess*procSqrt+proc;
            int rows = BLOCK_SIZE(procNumber, numprocs, A.col); // ILE PROCES procNumber PRZECZYTAL WIERSZY MACIUERZTY
            int i = 0;
            int row = 0;
            for(row = 0; row<rows; row++)
            {
                    currentRow++;

                    MPI_Recv(mA[currentRow-1], myRows, MPI_TYPE, procNumber, row, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            }
        } 
    }
    
    //WYSYLANE MACIRZY B, ANALOGICZNIE JAK A
        currentRowForSend = 0;
        currentRow = 0;
        //TERAZ TRZEBA WYSLAC WSZYSTKIE POSIADANE WIERSZE DO WSZYSTKICH PROCESOWW!
        rows = BLOCK_SIZE(myid, numprocs, B.col); // ILE PRZECZYTALEM WIERSZY MACIERZY
        i = 0;
        row = 0;
        for(row = 0; row<rows; row++)
        {
            for(i = 0; i<B.col; i+=myRows)
            {
                int procInRow = i/myRows; //NUMER PROCESU W WIERSZU DO WYSLANIA
                if(procInRow == meProc){
                    continue;
                }else{
                    int procNumber = rowOfProcess*procSqrt+procInRow;
                    MPI_Send(&b[currentRowForSend][i], myRows, MPI_TYPE, procNumber, row, MPI_COMM_WORLD);
                }
                
            }  
            currentRowForSend++;
        }
    
    
    for(proc = 0; proc < procSqrt; proc ++)
    {
        if(meProc == proc)
        { 
            int rows = BLOCK_SIZE(myid, numprocs, A.col); // ILE PRZECZYTALEM WIERSZY MACIERZY
            int i = 0;
            int row = 0;
            for(row = 0; row<rows; row++)
            {
                for(i = 0; i<myRows; i++)
                {
                    mB[currentRow][i] = b[row][colOffset+i];
                }   
                currentRow++;
            }
        }
        else
        {
            int procNumber = rowOfProcess*procSqrt+proc;
            int rows = BLOCK_SIZE(procNumber, numprocs, B.col); // ILE PROCES procNumber PRZECZYTAL WIERSZY MACIUERZTY
            int i = 0;
            int row = 0;
            for(row = 0; row<rows; row++)
            {

                    currentRow++;

                    MPI_Recv(mB[currentRow-1], myRows, MPI_TYPE, procNumber, row, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            }
        } 
    }
    
        


    
    
    //SKEWING A
	//ODPOWIEDNIE PRZESUWANIE WIERSZY MACIERZY A W LEWO , PRZESUWANIE O WARTOSC WIERSZA
    int colOfProcess = myid%procSqrt;
        row = 0;
        for(row = 0; row < myRows; row++)
        {
            int places = rowOfProcess*rowsPerProc+ row;
            if(places > 0){
                int place = 0;
                for(place = 0; place<myRows; place++){
                    int leftProc = getProcLeftByDist(colOfProcess, place, places, matrixDim, procSqrt);

                    int procNumber = rowOfProcess*procSqrt+leftProc;

                    int newLeftPlace = getNewLeftPos(place, places, colsPerProc);
                    MPI_Send(&mA[row][place], 1, MPI_TYPE, procNumber, newLeftPlace, MPI_COMM_WORLD);
                }
            
                for(place = 0; place<myRows; place++){
                    int rightProc = getProcRightByDist(colOfProcess, place, places, matrixDim, procSqrt);
                    int procNumber = rowOfProcess*procSqrt+rightProc;


                    MPI_Recv(&mA[row][place]/*&test*/, 1, MPI_TYPE, procNumber, place, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                }
            }
        }

    
    
    
    
    //SKEWING B
	//ODPOWIEDNIE PRZESUWANIE KOLUMN MACIERZY B W GORE , PRZESUWANIE O WARTOSC KOLUMNY
    colOfProcess = myid%procSqrt;
        int col = 0;
        //UWAGA METODA PODOBNA JAK SKEWING A, Z DROBNYMI ZMIANAMI
        for(col = 0; col < myRows; col++)
        {
            int places = colOfProcess*rowsPerProc+ col;
            if(places > 0){
                int place = 0;
                for(place = 0; place<myRows; place++){
                    int topProc = getProcTopByDist(myid, place, places, matrixDim, procSqrt);

                    int procNumber = topProc;

                    int newLeftPlace = getNewLeftPos(place, places, colsPerProc);
                    MPI_Send(&mB[place][col], 1, MPI_TYPE, procNumber, newLeftPlace, MPI_COMM_WORLD);
                }
            
                for(place = 0; place<myRows; place++){
                    int downProc = getProcDownByDist(myid, place, places, matrixDim, procSqrt);
                    int procNumber = downProc;
                    MPI_Recv(&mB[place][col], 1, MPI_TYPE, procNumber, place, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
            }
        }
    
	
//UWAGA macierz nA i nB JEST TRAKTOWANA JAKO NOWA MACIERZ TYMCZASOWA
//WPROWADZONA W CELU PRZESUWANIA PO JEDNYM MACIERZY A I B, PONIEWAZ MACIERZE A I B NIE SA PRZESYLANE PO JEDNYM, LECZ PELNYMI PARTIAMI	
	
    for(i=0; i<myRows ; i++)
            {
                for(j=0; j<myRows ; j++)
                {
                    nA[i][j] = mA[i][j];
                    nB[i][j] = mB[i][j];
                }
            }
    




	if(numprocs != 1)
	{
		//GDY LICZBA PROCESOW WIEKSZA NIZ 1
		for(l = 0; l < procSqrt; l++)
		{
			
			for(i=0; i<myRows ; i++)
            {//PRZEPISANIE WARTOSCI ODEBRANYCH DO TYMCZASOWEJ
				//POTRZEBNA SA TE STARE WARTOSCI, BO DO mA ZOSTANA ODEBRANE NOWE, WIEC TRZEBA ZNAC STARE ZEBY TE NOWE PRZESUWAC PO JEDNYM
                for(j=0; j<myRows ; j++)
                {
                    nA[i][j] = mA[i][j];
                    nB[i][j] = mB[i][j];
                }
            }   
            
			if(l==0)
			{//ZWYKLE LICZENIE ALE TYLKO W PIERWSZYM KROKU!!
				int k = 0;
				for(i=0; i<myRows ; i++)
				{
					for(j=0; j<myRows ; j++)
					{
						res[i][j] += nA[i][j] * nB[i][j];
					}
				}
			}	
        

			//PRZESLANIE CALEGO BLOKU DANYCH DO PROCESOW W LEWO I W GORE
            int procTop = getProcTop(myid, procSqrt);
            int procDown = getProcDown(myid, procSqrt);
            int procLeft = getProcLeft(myid, procSqrt);
            int procRight = getProcRight(myid, procSqrt);
            MPI_Send(&(mA[0][0]), myRows*myRows, MPI_TYPE, procLeft, 0, MPI_COMM_WORLD);
            MPI_Recv(&(mA[0][0]), myRows*myRows, MPI_TYPE, procRight, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Send(&(mB[0][0]), myRows*myRows, MPI_TYPE, procTop, 0, MPI_COMM_WORLD);
            MPI_Recv(&(mB[0][0]), myRows*myRows, MPI_TYPE, procDown, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		
		
		
		
			int k = 0;
			int maxIter = l==procSqrt-1 ? rowsPerProc - 1 : rowsPerProc; //W OSTATNIM KROKU NIE ROBIMY LICZENIA, BO ZOSTALO ZROBIONE NA SAMYM POCZATKU (GDY l==0)
			for(k = 0; k<maxIter; k++)
			{
			
				for(i=0; i<myRows ; i++)
				{
					for(j=0; j<myRows -1 ; j++)
					{
						nA[i][j] = nA[i][j+1]; //DOSUNIECIE JEDNEJ KOLUMNY A W LEWO
					}
					
				}
            
				for(i=0; i<myRows -1 ; i++)
				{
					for(j=0; j<myRows ; j++)
					{
						nB[i][j] = nB[i+1][j]; //DOSUNIECIE JEDNEJ KOLUMNY B W GORE
					}
					
				}
            
				for(i=0; i< myRows ; i++)
				{
					nA[i][myRows-1] = mA[i][k]; //OSTATNIEJ KOLUMNY NIE MAMY, MUSIMY POBRAC Z OTRZYMANYCH mA
					nB[myRows-1][i] = mB[k][i]; //OSTATNIEGO WIERSZA NIE MAMY, MUSIMY DOSUNAC Z OTRZYMANYCH mA
				}
			
				for(i=0; i<myRows ; i++)
				{
					for(j=0; j<myRows ; j++)
					{
						res[i][j] += nA[i][j] * nB[i][j]; //LICZENIE
					}
				}
				
			}

			for(i=0; i<myRows ; i++)
            {//GDY JUZ PRZESUNIETO WSZYSTKIE WARTOSCI, NALEZY JE WPISAC DO mA i mB ABY MOZNA BYLO JE WYSLAC W KOLEJNEJ ITERACJI DO POZOSTALYCH PROCESOR
                for(j=0; j<myRows ; j++)
                {
                    mA[i][j] = nA[i][j];
                    mB[i][j] = nB[i][j];
                }
            }
		}
 
		
	}else{
		//OSOBNA WERSJA ALGORYTMU DLA JEDNEGO PROCESU
		//ZWYKLE PRZESUWANIE I MNOZENIE
		int l = 0;
		for(l = 0; l < myRows ; l++)
		{
			for(i=0; i<myRows ; i++)
				{
					for(j=0; j<myRows ; j++)
					{
						nA[i][j] = mA[i][j];
						nB[i][j] = mB[i][j];
					}
				}

            for(i=0; i<myRows ; i++)
            {
                for(j=0; j<myRows ; j++)
                {
                    res[i][j] += nA[i][j] * nB[i][j];
                }
            }

            //trzeba w lewo i w gore o jeden
            for(i=0; i<myRows ; i++)
            {
                for(j=0; j<myRows ; j++)
                {
                    nA[i][j] = mA[i][j];
                    nB[i][j] = mB[i][j];
                }
            }
            
            
            for(i=0; i<myRows ; i++)
            {
                for(j=0; j<myRows -1 ; j++)
                {
                    nA[i][j] = mA[i][j+1];
                }
                
            }
            
            for(i=0; i<myRows -1 ; i++)
            {
                for(j=0; j<myRows ; j++)
                {
                    nB[i][j] = mB[i+1][j];
                }
                
            }
            
            for(i=0; i< myRows ; i++)
            {
                nA[i][myRows-1] = mA[i][0];
                nB[myRows-1][i] = mB[0][i];
            }
            
            for(i=0; i<myRows ; i++)
            {
                for(j=0; j<myRows ; j++)
                {
                    mA[i][j] = nA[i][j];
                    mB[i][j] = nB[i][j];
                }
            }
		}
	}		
    

    
	
	
	//NASTEPUJE ZEBRANIE WSZYSTKICH WYNIKOW DO PROCESU 0       
    int size = RES.col*RES.row;
    if(size > numprocs) size = numprocs;
	
    double** result;
    if(myid==0){
		result = alloc_2d(matrixDim, matrixDim);
		double** resrec = alloc_2d(myRows, myRows);
        for(i = 1; i < numprocs; i++){
            TYPE res;
            MPI_Recv(&(resrec[0][0]), myRows*myRows, MPI_TYPE, i, 0, MPI_COMM_WORLD,
                MPI_STATUS_IGNORE);
			int k, l;
			int procRow = i/procSqrt;
			int procCol = i%procSqrt;
			int jOffset = procCol*rowsPerProc;
			int iOffset = procRow*rowsPerProc;
			for(k=0; k<myRows; k++)
			{
				for(l=0; l<myRows; l++)
				{
					result[iOffset+k][jOffset+l] = resrec[k][l];
				}
			}
        }
		int k, l;
		for(k=0; k<myRows; k++)
			{
				for(l=0; l<myRows; l++)
				{
					result[k][l] = res[k][l];
				}
			}
    }else if(myid < size){
    	
        MPI_Send(&(res[0][0]), myRows*myRows, MPI_TYPE, 0, 0, MPI_COMM_WORLD);
    }
	
	if(myid==0){
		printf("RESULT IS: \n");
		for(i = 0; i<myRows*procSqrt; i++)
		{
			for(j = 0; j<myRows*procSqrt; j++)
			{
				printf("%f ", result[i][j]);
			}
			printf("\n");
		}
	}
    if(myid == 0){
    	FILE *f = fopen("result_matrix.txt", "w+");
    	if(f == NULL){
    		printf("Error opening file!\n");
    	} else {
    		fprintf(f, "%d %d ", myRows*procSqrt, myRows*procSqrt);
    		for(i = 0; i < myRows*procSqrt; ++i){
		        for(j = 0; j < myRows*procSqrt; ++j){
		            fprintf(f,"%.2f ", result[i][j]);
		        }
		    }
    	}
    	fclose(f);
    }	

    MPI_Finalize();
    return 0;
    /*
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
    
    return 0;*/
}

