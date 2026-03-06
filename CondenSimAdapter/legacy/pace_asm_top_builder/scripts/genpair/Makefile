CFLAGS = -O2 -Wall -W -g

all:genPairPACE

genPairPACE: genPairPACE.o  string_op.o file_op.o parse_cg.o geom_op.o numrecp.o 
	gcc -o $@ $^ -lm 

%.o: %.c
	gcc -c $(CFLAGS) $<

depend:
	gcc -MM *.c > $@

clean:
	rm -f *.o

-include depend

