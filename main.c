#include <upc.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define mpitype MPI_DOUBLE
#define MPI_TYPE MPI_DOUBLE
#define TYPE double

//kompilacja i uruchomienie:
// source /opt/nfs/config/source_bupc.sh
// upcc -gupc -network=mpi main.c -o main
// upcrun -n 4 ./main matrix.txt matrix.txt 
//w miejsce -n 4 WPISUJEMY ILOSC WATKOW

//ZMIENNE GLOBALNE, RACZEJ POWINNO SIE CHYBA Z NIMI COS ZROBIC
int rows_global;
int local_rows_global;


//lokalne przesuniecie w bloku.
//funkcja zwraca dla danego k, ktory jest to element w danym bloku
int localOffset(int k, int local_rows)
{
	int offsetByJHorizontal = (k%(local_rows*local_rows))%local_rows;
	
	int offsetByIVertical = (k%(local_rows*local_rows))/local_rows;
	
	return offsetByJHorizontal + offsetByIVertical*local_rows;
}

//map from real to upc
//metoda mapuje 'realne' polozenie indeksu do polozenia w 'upc' indeksu
/*przyklad:
	dla macierzy 4x4 i 4 procesow:
	realne:0 upc:0 (blok 0)
	realne:1 upc:4 (blok 0)
	realne:2 upc:8 (blok 0)
	realne:3 upc:12 (blok 0)
	realne:4 upc:1 (blok 1)
	realne:5 upc:5 (blok 1)
	realne:6 upc:9 (blok 1)
	realne:7 upc:13 (blok 1)
	realne:8 upc:2 (blok 2)
	realne:9 upc:6 (blok 2)
*/
int mp(int k)
{
	int blockNumber = (int)(k/(local_rows_global*local_rows_global));
	return (localOffset(k, local_rows_global))*THREADS + blockNumber;
}

//funkcja zwraca pozycje w wektorze na podstawie pozycji w macierzy
//rows : rozmiar macierzy
//local_rows : rozmiar macierzy danego procesu(bloku pojedynczego)
//i - pozycja w macierzy w pionie; j - pozycja w poziomie
int indMat2Vec(int rows, int local_rows, int i, int j)
{
	//offsetByBlockVertical : jaki offset w wektorze wynikowym powoduje polozenie bloku wzgledem pionu
	// i/local_rows : numer bloku liczac od gory
	// rows*local_rows : ilosc wszystkich elementow macierzy w jednym wierzu blokow
	int offsetByBlockVertical = (i/local_rows)*rows*local_rows;
	
	//offsetByBlockHorizontal : jaki offset w wektorze wynikowym powoduje polozenie bloku wzgledem poziomu
	// i/local_rows : numer bloku liczac od gory
	// local_rows*local_rows : ilosc wszystkich elementow macierzy w jednym bloku
	int offsetByBlockHorizontal = (j/local_rows)*local_rows*local_rows;
	
	//offset przez i w danym bloku
	int offsetByIVertical = (i%local_rows)*local_rows;
	
	//offset przez j w danym bloku
	int offsetByJHorizontal = (j%local_rows);
	
	//trzeba jeszcze zmapowac do indeksu upc
	return mp(offsetByBlockVertical + offsetByBlockHorizontal + offsetByIVertical + offsetByJHorizontal);
}

//PASOWALOBY UZYWANE WARTOSCI Z TEJ METODY ZAPISAC DO TABLICY PRZY URUCHAMIANIU PROGRAMU
//KAZDY PROCES MA SWOJA LOKALNA TABLICE
//I KAZDY PROCES UZYWALBY TEJ TABLICY ZAMIAST ZA KAZDYM RAZEM WYWOLYWAC TE FUNKCJE
//BO JEST DOSC SPORO OBLICZEN, ZWAZYWSZY ZE W KAZDEJ ITERACJI PO KILKA RAZY DLA KAZDEGO ELEMENTU MACIERZY JEST WYWOLYWANA TA FUNKCJA
int indM2V(int i, int j)
{
	return indMat2Vec(rows_global, local_rows_global, i, j);
}

int getBlockLeft(int block, int procsDim){
    if(block%procsDim == 0)
        return block+procsDim-1;
    else
        return block-1;
}

int getBlockTop(int block, int procsDim){
    block = block - procsDim;
    if(block >=0)
        return block;
    else
        return procsDim*procsDim-fabs(block);
}

//zwraca indeks 'realny' pierwszego elementu z bloku
int getFirstElementFromBlock(int block)
{
	return block*local_rows_global*local_rows_global;
}


void read_row_striped_matrix(
    char *s, // IN file name
	shared  TYPE ** storage,
    shared int *m, //matrix rows
    shared int *n //matrix cols
	)
{
    int i;
	int j;
    FILE *infileptr; //input file pointer
    shared int* local_rows = (shared int *) upc_all_alloc(THREADS, sizeof(int)); // rows on this proc
    int p; // number of processes
	p = THREADS;
	int procsqrt = sqrt(THREADS);


    if (MYTHREAD == 0){
        infileptr = fopen(s, "r");
        if(infileptr == NULL) *m = 0;
        else{
			int mm;
            fscanf(infileptr, "%d", &mm);
			*m = mm;
			int nn;
            fscanf(infileptr, "%d", &nn);
			*n = nn;
        }
		
		if(!(*m)) printf("Blad pliku! \n\n");
    }

    
	if(MYTHREAD == 0){
		*local_rows = (*m)/procsqrt;
	}


	upc_barrier;
	local_rows_global = *local_rows;
	rows_global = *m;

	*storage = (shared  TYPE*) upc_all_alloc( THREADS, local_rows_global*local_rows_global*sizeof(TYPE));
	
	

    if(MYTHREAD == 0){
		for (i = 0 ; i< (*n); i++) {
			for (j = 0; j<(*m); j++){
				double value;
				fscanf(infileptr, "%lf", &value);
				(*storage)[indMat2Vec((*n), *local_rows, i, j)] = value;
			}
		}
        
        fclose (infileptr);
    }
	
	upc_barrier;
}

void printMatrix(int size, shared  TYPE* matrix)
{
	int i, j;
	if (MYTHREAD == 0) {
		for (i=0; i<size; i++){
			for(j=0; j<size; j++){
				printf("\n\n");
				printf(" , %f ,  ", matrix[indM2V(i,j)]);
			}
		}
		printf("\n\n");
	}
}

void saveMatrix(int size, shared TYPE* matrix) {
	FILE *f = fopen("result_matrix.txt", "w+");
	if(f == NULL){
		printf("Error opening file!\n");
	} else {
		fprintf(f, "%d %d ", size, size);
		int i,j;
		for(i = 0; i < size; ++i){
	        for(j = 0; j < size; ++j){
	            fprintf(f,"%.0f ", matrix[indM2V(i,j)]);
	        }
	    }
	}
	fclose(f);
}

void copyMatrix(int size, shared  TYPE** toCopy, shared  TYPE** toPaste)
{
	int i, j;
	upc_forall (i=0; i<size; i++;&(*toCopy)[i]){
		for(j=0; j<size; j++){
			(*toPaste)[indM2V(i, j)] = (*toCopy)[indM2V(i, j)];
		}
	}
	upc_barrier;
}

int main(int argc, char *argv[]) {
  int i;
  int j;
  int k;
  int procSqrt = sqrt(THREADS);
  shared  TYPE * a;
  shared  TYPE * b;
  shared TYPE * res;
  shared TYPE * at;
  shared TYPE * bt;
  shared int* ms = (shared int *) upc_all_alloc(THREADS, sizeof(int));
  shared int* ns = (shared int *) upc_all_alloc(THREADS, sizeof(int));
  read_row_striped_matrix(argv[1], &a, ms, ns);
  read_row_striped_matrix(argv[2], &b, ms, ns);
  int size = *ms;
  
  shared int ** indM = (shared int **) upc_all_alloc(THREADS, size * sizeof(TYPE*));
  for(i=0;i<size;++i){
  	indM[i] = (shared int *) upc_all_alloc(THREADS, size * sizeof(TYPE));
  }
  
  for(i=0;i<size;++i){
  	for(j=0;j<size;++j){
  		indM[i][j] = indM2V(i,j);
  	}
  }
  upc_barrier;
  
  res = (shared TYPE*) upc_all_alloc(THREADS, size*size*sizeof(TYPE));
  at = (shared TYPE*) upc_all_alloc(THREADS, size*size*sizeof(TYPE));
  bt = (shared TYPE*) upc_all_alloc(THREADS, size*size*sizeof(TYPE));
  
  upc_barrier;
  
  
  //przesuwanie macierzy A w lewo
	for(i=1; i<size; i++){
		upc_forall(k=0; k<i; k++;&a[indM[i][0]]){
			TYPE firstValue = a[indM[i][0]];
			for(j=0; j<size-1; j++){
				a[indM[i][j]] = a[indM[i][j+1]];
			}

			a[indM[i][size-1]] = firstValue;

		}
	}
  upc_barrier;  
  //przesuwanie macierzy B w gore
	for(j=1; j<size; j++){
		upc_forall(k=0; k<j; k++;&b[indM[0][j]]){
			TYPE firstValue = b[indM[0][j]];
			for(i=0; i<size-1; i++){
				b[indM2V(i, j)] = b[indM[i+1][j]];
			}
			b[indM[size-1][j]] = firstValue;
		}
	}
  upc_barrier;
	  copyMatrix(size, &a, &at);
	  copyMatrix(size, &b, &bt);
  
	upc_forall( i = 0; i< size*size; i++; &a[i]){
		res[i] = 0;
	}
  
  upc_barrier;
  
  int l;
  int el;
  for(l = 0; l < procSqrt; l++){
	  
	  
	  copyMatrix(size, &at, &a);
	  copyMatrix(size, &bt, &b);
	  
	if(l==0)
	{//ZWYKLE LICZENIE ALE TYLKO W PIERWSZYM KROKU!!
		upc_forall( k = 0; k< size*size; k++; &b[k]){
			res[k] += a[k] * b[k];
		}
	}
	
	upc_barrier;
	  

	//PRZESLANIE CALEGO BLOKU DANYCH DO PROCESOW W LEWO I W GORE
	for(k = 0; k < THREADS; k++){
		int blockLeft = getBlockLeft(k, procSqrt);
		int blockTop = getBlockTop(k, procSqrt);

		int firstElLeft = getFirstElementFromBlock(blockLeft);
		int firstElTop = getFirstElementFromBlock(blockTop);
		int firstElCurrent = getFirstElementFromBlock(k);
		
		upc_forall(i = 0; i<local_rows_global*local_rows_global; i++;&at[mp(i)])
		{
			at[mp(firstElLeft+i)] = a[mp(firstElCurrent + i)];
		}
		
		upc_forall(i = 0; i<local_rows_global*local_rows_global; i++;&bt[mp(i)])
		{
			bt[mp(firstElTop+i)] = b[mp(firstElCurrent + i)];
		}

	}
	
	  upc_barrier;
	
	int k = 0;
	int proc;
	int maxIter = l==procSqrt-1 ? local_rows_global - 1 : local_rows_global; //W OSTATNIM KROKU NIE ROBIMY LICZENIA, BO ZOSTALO ZROBIONE NA SAMYM POCZATKU (GDY l==0)
	for(k = 0; k<maxIter; k++)
	{
		
		upc_forall ( i=0; i<size*size-1; i++; &a[mp(i)]){
			a[mp(i)] = a[mp(i+1)]; //DOSUNIECIE JEDNEJ KOLUMNY A W LEWO
		}
		
		upc_forall ( i=0; i<size*size-local_rows_global; i++; &b[mp(i)]){
			b[mp(i)] = b[mp(i+local_rows_global)]; //DOSUNIECIE JEDNEJ KOLUMNY B W GORE
		}

		upc_barrier;

		
		for(proc = 0 ;proc < THREADS; proc++)
		{
			int firstElCurrent = getFirstElementFromBlock(proc);
			upc_forall(i=0; i<local_rows_global; i++;&a[i]){
				a[mp(firstElCurrent + (local_rows_global)*i + (local_rows_global-1))] = at[mp( firstElCurrent + (local_rows_global)*i + k)];
				//OSTATNIEJ KOLUMNY A NIE MAMY W BLOKACH, MUSIMY POBRAC Z ATEMP
			}
		}
		
		
		for(proc = 0 ;proc < THREADS; proc++)
		{
			int firstElCurrent = getFirstElementFromBlock(proc);
			upc_forall(i=0; i<local_rows_global; i++;&b[i]){
				b[mp(firstElCurrent + i + (local_rows_global*(local_rows_global-1)))] = bt[mp(firstElCurrent + i + k*local_rows_global)]; 
				//OSTATNIEGO WIERSZA B NIE MAMY W BLOKACH, MUSIMY POBRAC Z BTEMP
			}
		}

		upc_barrier;
		
		upc_forall( i = 0; i< size*size; i++; &a[i]){
			//LICZENIE
			res[i] += a[i] * b[i];
		}
            
		upc_barrier;	
				
	}
	
	upc_forall( i = 0; i< size*size; i++; &a[i]){
			//GDY JUZ PRZESUNIETO WSZYSTKIE WARTOSCI, NALEZY JE WPISAC DO at i bt ABY MOZNA BYLO JE WYSLAC W KOLEJNEJ ITERACJI DO POZOSTALYCH PROCESOR
			at[i] = a[i];
			bt[i] = b[i];
		}

	
	
	
	
  }
  
  
  upc_barrier;
  if(MYTHREAD == 0)
  {
	printMatrix(size, res);
  	saveMatrix(size, res);
  }
  
  upc_barrier;
  return 0;
} 
