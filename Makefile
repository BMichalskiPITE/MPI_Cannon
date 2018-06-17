all:
	make clean; make run;

run:
	upcc -gupc -network=mpi main.c -o main
	upcrun -n $(PROC) ./main $(MATRIX1) $(MATRIX2)
clean:
	rm main
