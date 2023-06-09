all:	bec3pc

bec3pc:	bec3p.c++ parameters3.h
	c++ -O3 bec3p.c++ -o bec3pc
