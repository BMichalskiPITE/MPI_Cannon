# MPI Cannon

## Run

kompilacja i uruchomienie:

wprowadzenie zmiennych srodowiskowych:
 source /opt/nfs/config/source_bupc.sh
 
kompilacja: 
 upcc -gupc -network=mpi main.c -o main

uruchomienie:
 upcrun -n 4 ./main matrix.txt matrix.txt 
w miejsce -n 4 WPISUJEMY ILOSC WATKOW
