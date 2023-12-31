CC=gcc
CFFLAGS = -Wall -Wextra -g3 #-Wall -O3
OBJS = main.o tsp.o chrono.o utilities.o convex_hull.o
DEPS = tsp.h
RM = rm -rf
EXE = tsp

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

main: $(OBJS)
	$(CC) -o $@ $^ $(CFLAGS)


### gcc -Wall -Wextra -g3 /Users/elnur/Desktop/OR2/project/main.c -o 
### /Users/elnur/Desktop/OR2/project/output/main 