# MPI Cannon

## Comile
```
mpicc main.c -o main
```

##Run
```
mpiexec -n 1 ./main matrix.txt matrix2.txt
mpiexec -n 4 ./main matrix.txt matrix2.txt
mpiexec -n 9 ./main matrix.txt matrix2.txt
mpiexec -n 36 ./main matrix.txt matrix2.txt
```
