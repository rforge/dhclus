# make
#
CC = gcc
#CFLAGS = -mavx -fpic -O3 -ftree-vectorize -std=gnu99 -fgnu89-inline -w
CFLAGS= -fpic -O3 -std=gnu99 -w
SHFLAG = -shared
HDIR = /usr/share/R/include
LDIR = /usr/lib/R/lib
LIBS = -lR
input = *.c
obj = *.o
icav.so: $(obj)
	$(CC) -shared $(obj) -L$(LDIR) -o $@ $(LIBS)
$(obj):$(input)
	$(CC) -I$(HDIR) -L$(LDIR) $(CFLAGS) -c $(input)
clean:
	rm -rf *.o 
