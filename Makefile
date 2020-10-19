all: 
	make -C src
	make -C external/gsufsort/ DNA=1 TERMINATOR=0 

clean:
	make clean -C src 
	make clean -C external/gsufsort/
