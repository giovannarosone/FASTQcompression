VLIB= -g -O0

MY_CXX_FLAGS= -std=c++11 #-Wall -Wextra -DNDEBUG
MY_CXX_OPT_FLAGS= -O3 -m64 
MY_CXX=g++

M=0
MODE = -DM=$(M) 
B=0
BIN = -DB=$(B)
LFLAGS = -lm -ldl

CXX_FLAGS=$(MY_CXX_FLAGS) $(MY_CXX_OPT_FLAGS) $(LFLAGS) $(MODE) $(BIN)

####

INPUT=../../dataset/reads.100.fastq
OUTPUT=../../output.fq

####


all: main 

main: FASTQcompression.cpp  
	$(MY_CXX) FASTQcompression.cpp -o fastqcompression  $(CXX_FLAGS) 

clean:
	rm -f fastqcompression 

run:
	./fastqcompression -e $(INPUT).bwt -q $(INPUT).bwt.qs -f $(INPUT).bwt -o $(OUTPUT) 
	
valgrind: 
	$(MY_CXX) FASTQcompression.cpp -o fastqcompression-test $(VLIB)
	valgrind --tool=memcheck --leak-check=full --track-origins=yes ./fastqcompression-test -e $(INPUT).bwt -q $(INPUT).bwt.qs -f $(INPUT).bwt -o $(OUTPUT) 

