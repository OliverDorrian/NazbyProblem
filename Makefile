# Variables
CC = gcc
CFLAGS = -fopenmp -O2
LDFLAGS = -fopenmp

target: cInsertion.c fInsertion.c ompcInsertion.c ompfInsertion.c coordReader.h

ci: cInsertion.c coordReader.h
	gcc cInsertion.c coordReader.h -o ci.exe -lm

fi: fInsertion.c coordReader.h
	gcc fInsertion.c coordReader.h -o ci.exe -lm

comp: ompcInsertion.c coordReader.h
	$(CC) $(CFLAGS) -o ompC.exe $^ $(LDFLAGS)  -lm

fomp: ompfInsertion.c coordReader.h
	$(CC) $(CFLAGS) -o ompF.exe $^ $(LDFLAGS)  -lm

icomp: ompcInsertion.c coordReader.h
	$(CC) $(CFLAGS) -o ompCI.exe $^ $(LDFLAGS)  -lm

ifomp: ompfInsertion.c coordReader.h
	$(CC) $(CFLAGS) -o ompFI.exe $^ $(LDFLAGS)  -lm
