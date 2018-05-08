# MPI Cannon

## Run
Makefile uruchamia za nas wszystko, trzeba podać jedynie nazwy plikow z macierzami
```
make run MATRIX1=<nazwa macierzy A> MATRIX2=<nazwa macierzy B>

make run MATRIX1=matrix512_1.txt MATRX2=matrix512_2.txt
```

## Generator macierz
tworzy pliki matrix(rozmiar)_(1|2).txt zgodne z oczwkiwanym formatem programu mpi
```
python3 generate.py <rozmiar>

python3 generate.py 512
```
