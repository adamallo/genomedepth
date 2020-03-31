#Config
ifndef CC
CC=gcc
endif
CFLAGS:= $(CFLAGS) -Wall
PERF=-O3
DBG=-DDEBUG -O0
SRCDIR=./src
BINDIR=./bin
EXECUTABLES=genomedepth

#Internals targets
$(BINDIR)/genomedepth: $(SRCDIR)/main.c
	mkdir -p $(BINDIR)
	$(CC) $(CFLAGS) $(LDFLAGS) $(PERF) $^ -o $@ $(C_LIBS) $(D_LIBS)
	@echo "\ngenomedepth built"
$(BINDIR)/genomedepth_dbg: $(SRCDIR)/main.c
	mkdir -p $(BINDIR)
	$(CC) $(CFLAGS) $(LDFLAGS) $(DBG) -g $^ -o $@ $(C_LIBS) $(D_LIBS)
	@echo "\nDebug version of genomedepth built"
clean:
	@echo "Cleaning debug files directory"
	rm -rf genomedepth_dbg.dSYM

#External targets
genomedepth: $(BINDIR)/genomedepth
genomedepth_debug: $(BINDIR)/genomedepth_dbg
genomedepth_dbg: genomedepth_debug
debug: genomedepth_dbg
all: $(EXECUTABLES)
.PHONY: clean
