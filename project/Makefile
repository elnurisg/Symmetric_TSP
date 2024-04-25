CC = gcc
CFLAGS = -Wall -Wextra -g3 -I./include -I/Applications/CPLEX_Studio2211/cplex/include/ilcplex
LDFLAGS = -L/Applications/CPLEX_Studio2211/cplex/lib/x86-64_osx/static_pic -lcplex -lm -lpthread -ldl
RM = rm -rf
EXE = tsp

SRCDIR = src
INCDIR = include

OBJS = $(SRCDIR)/main.o $(SRCDIR)/tsp.o $(SRCDIR)/heuristics.o $(SRCDIR)/metaheuristics.o $(SRCDIR)/matheuristics.o $(SRCDIR)/chrono.o $(SRCDIR)/utilities.o $(SRCDIR)/convex_hull.o
DEPS = $(INCDIR)/tsp.h $(INCDIR)/heuristics.h $(INCDIR)/metaheuristics.h $(INCDIR)/matheuristics.h $(INCDIR)/utilities.h $(INCDIR)/convex_hull.h

$(SRCDIR)/%.o: $(SRCDIR)/%.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

$(EXE): $(OBJS)
	$(CC) -o $@ $^ $(CFLAGS) $(LDFLAGS) 2>&1 | grep -v "libcplex\.a"

clean:
	$(RM) $(EXE) $(OBJS)
