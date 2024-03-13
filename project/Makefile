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

# # ifeq ($(OS),Darwin)
# setting = 0
# CPLEX_HOME = /Users/elnur/Downloads/CPLEX
# CC = clang -Qunused-arguments
# AR = ar rc
# LIBS = -Wl,-no_compact_unwind -L${CPLEX_HOME}/lib/x86-64_osx/static_pic -L. -lcplex -lm -lpthread 
# INC = -I. -I${CPLEX_HOME}/include/ilcplex
# # endif

# $(EXE): $(OBJS) $(LIBUTILS)
# 	$(CC) $(CFLAGS) -o $(EXE) $(OBJS) $(LIBS)

# $(OBJS) : $(HEADERS)

# $(LIBUTILS): $(OBJS_LIBUTILS)
# 	$(AR) $(LIBUTILS) $(OBJS_LIBUTILS)

# $(LIBUTILS) : $(HEADERS_LIBUTILS)
