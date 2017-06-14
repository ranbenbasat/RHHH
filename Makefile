CC=gcc -O2 -Wall
all: RandHHH 10RandHHH RandHHH_33 10RandHHH_33 RandHHH2D 10RandHHH2D check1 check1_33 check2
RandHHH: randhhh1D.c randhhh1D.h randmain1.c ulossycount.h ulossycount.c hashtable.h hashtable.c prng.h prng.c Makefile
	$(CC) -DRHHH randhhh1D.c randmain1.c ulossycount.c hashtable.c prng.c -o RandHHH -lm
10RandHHH: randhhh1D.c randhhh1D.h randmain1.c ulossycount.h ulossycount.c hashtable.h hashtable.c prng.h prng.c Makefile
	$(CC) -DRHHH -DVMULT=10 randhhh1D.c randmain1.c ulossycount.c hashtable.c prng.c -o 10RandHHH -lm
RandHHH_33: randhhh1D.c randhhh1D.h randmain1.c ulossycount.h ulossycount.c hashtable.h hashtable.c prng.h prng.c Makefile
	$(CC) -DRHHH -DNUM_MASKS=33 randhhh1D.c randmain1.c ulossycount.c hashtable.c prng.c -o RandHHH_33 -lm
10RandHHH_33: randhhh1D.c randhhh1D.h randmain1.c ulossycount.h ulossycount.c hashtable.h hashtable.c prng.h prng.c Makefile
	$(CC) -DRHHH -DVMULT=10 -DNUM_MASKS=33 randhhh1D.c randmain1.c ulossycount.c hashtable.c prng.c -o 10RandHHH_33 -lm
RandHHH2D: randhhh2D.c main2.c randhhh2D.h ulossycount.c ulossycount.h prng.h prng.c Makefile
	$(CC) -DRANDHHH -DDIMENSION2 randhhh2D.c ulossycount.c prng.c main2.c -o RandHHH2D -lm
10RandHHH2D: randhhh2D.c main2.c randhhh2D.h ulossycount.c ulossycount.h prng.h prng.c Makefile
	$(CC) -DRANDHHH -DVMULT=10 -DDIMENSION2 randhhh2D.c ulossycount.c prng.c main2.c -o 10RandHHH2D -lm
	
check1: check1.cpp Makefile
	g++ -ansi -Wall -O2 check1.cpp -o check1 -lm
check1_33: check1.cpp Makefile
	g++ -ansi -Wall -O2 -DNUM_MASKS=33 check1.cpp -o check1_33 -lm
check2: check2.cpp Makefile
	g++ -ansi -Wall -O2 check2.cpp -o check2 -lm
	