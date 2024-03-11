CC = gcc
CFLAGS = -g -Wall -Wextra -pedantic

OBJ = Driver.o SparseMat.o

Driver.ex: $(OBJ)
	$(CC) $(CFLAGS) -o $@ $(OBJ) -lm

%.o: %.c
	$(CC) $(CFLAGS) -c $<

.PHONY: clean

clean:
	rm -f Driver.ex *.o *~
