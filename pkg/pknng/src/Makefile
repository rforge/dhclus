# make

#gcc
CC = gcc
#gcc flags
CFLAGS = -fpic -O3 -std=gnu99 -w
SHFLAG = -shared
HDIR = /usr/share/R/include
LDIR = /usr/lib/R/lib
#gcc LIBS
LIBS = -lR -lgomp -lpthread
input = *.c
obj = *.o
pknng.so: $(obj)
	$(CC) -shared $(obj) -L$(LDIR) -o $@ $(LIBS)
$(obj):$(input)
	$(CC) -I$(HDIR) -L$(LDIR) $(CFLAGS) -c $(input)
