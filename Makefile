files := main.cc Mesh.cc
CC := g++

all: file1 file2 file

file1 main.cc:
	${CC} -c main.cc -o main.o

file2 file2.cc:
	${CC} -c Mesh.cc -o Mesh.o -I$(eigen)

file main.o Mesh.o:
	${CC} main.o Mesh.o -o main 

clean:
	rm -f main.o Mesh.o