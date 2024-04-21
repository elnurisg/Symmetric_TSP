CC = gcc
CFLAGS = -Wall -Wextra -g3 -I/Applications/CPLEX_Studio2211/cplex/include/ilcplex
LDFLAGS = -L/Applications/CPLEX_Studio2211/cplex/lib/x86-64_osx/static_pic -lcplex -lm -lpthread -ldl
RM = rm -rf
EXE = tsp

OBJS = main.o tsp.o heuristics.o metaheuristics.o matheuristics.o chrono.o utilities.o convex_hull.o
DEPS = tsp.h heuristics.h metaheuristics.h matheuristics.h utilities.h convex_hull.h

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

$(EXE): $(OBJS)
	$(CC) -o $@ $^ $(CFLAGS) $(LDFLAGS) 2>&1 | grep -v "libcplex\.a"

clean:
	$(RM) $(EXE) $(OBJS)

