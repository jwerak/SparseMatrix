CC := g++
CFLAGS := -O3

BUILDFILE := matice

SOURCEFILES := *.cpp
OBJFILES := $(patsubst %.cpp,%.o,$(wildcard *.cpp))


all: $(OBJFILES)
	$(CC) $(OBJFILES) -o $(BUILDFILE)

$(OBJFILES): $(SOURCEFILES)
	$(CC) $(CFLAGS) $(SOURCEFILES) -c
	
clean: 
	rm -f $(OBJFILES)

