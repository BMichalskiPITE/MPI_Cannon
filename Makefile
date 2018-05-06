all:
	make clean; make run; make prof

run:
	/opt/nfs/mpich-3.2/bin/mpicc -o main main.c -I/opt/nfs/mpe2-2.4.9b/include -L/opt/nfs/mpe2-2.4.9b/lib -lmpe -lm -lpthread
	/opt/nfs/mpich-3.2/bin/mpiexec  -n 4 ./main $(MATRIX1) $(MATRIX2)
clean:
	rm mpe_logs.slog2 main mpe_logs.clog2
prof:
	/opt/nfs/mpe2-2.4.9b/bin/clog2TOslog2 mpe_logs.clog2 
	/opt/nfs/mpe2-2.4.9b/bin/jumpshot mpe_logs.slog2 
