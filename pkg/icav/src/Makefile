# make
#
CC = gcc
CFLAGS = -march=core2 -mavx -fpic -O3 -ftree-vectorize -std=gnu99 -fgnu89-inline -w
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