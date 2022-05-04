CC=g++
CFLAGS=-I. -O3
DEPS = utils.h bloomFilter.h assembler.h
OBJ = main.o bloomFilter.o utils.o assembler.o
EXEC = assembler

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

$(EXEC): $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS)

all: clean $(EXEC)

clean:
	rm $(OBJ) $(EXEC)
