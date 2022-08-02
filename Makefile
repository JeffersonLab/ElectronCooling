CC = g++
NAME = jspec
TARGET  = $(NAME)

LIBDIR = lib
INCDIR = include
CFLAGS = -O3 -Wall -std=c++11 -fPIC -I$(INCDIR) 
OMPFLAGS = 

#LIBS = -lm -lmuparser
#LINK = -L$(LIBDIR) -s -Wl,-rpath=$(LIBDIR) $(LIBS) 

LINK =  -L$(LIBDIR) -s -Wl,-rpath=$(LIBDIR) -lm -l:libmuparser.so.2 -lgsl -lgslcblas

SRC = $(wildcard src/*.cc)
OBJ = $(SRC:.cc=.o)
DEPS = $(wildcard $(INCDIR)/*.h)

#$(info $$SRC is [${SRC}])
#$(info $$OBJ is [${OBJ}])
#$(info $$CFLAGS is [${CFLAGS}])
#$(info $$OMPFLAGS is [${OMPFLAGS}])

.PHONE: all
all = $(TARGET)

$(TARGET): $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LINK) $(OMPFLAGS)
	
%.o: %.cc $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS) $(OMPFLAGS)

.PHONY: clean
clean:
	rm -f $(OBJ)
