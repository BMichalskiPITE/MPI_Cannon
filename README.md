# MPI Cannon

## Run

kompilacja i uruchomienie:

wprowadzenie zmiennych srodowiskowych:
 source /opt/nfs/config/source_bupc.sh
 
przez makefile:
make run MATRIX1=matrix.txt MATRIX2=matrix.txt PROC=4

lub recznie:

kompilacja: 
 upcc -gupc -network=mpi main.c -o main

uruchomienie:
 upcrun -n 4 ./main matrix.txt matrix.txt 
w miejsce -n 4 WPISUJEMY ILOSC WATKOW
