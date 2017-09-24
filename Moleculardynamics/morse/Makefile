CC=gcc
CFLAGS=-O3
EXE=main
OBJS= main.o  create_fcc.o veltemp.o  linmomscvelocity.o pbc.o morse.o vverlet.o potentialenergy.o initializevel.o kinenergy.o
LDFLAGS=-lm


.PHONY: all run clean

all: $(EXE)

$(EXE): $(OBJS)
	$(CC) -o $(EXE) $(OBJS) $(LDFLAGS)

run: clean $(EXE)
	./$(EXE)



clean:
	rm -f $(EXE) $(OBJS)

