# ParaPl - parallel platoon simulation (C11 + OpenMP)
#
#   make            build build/platoon
#   make test       determinism self-test (exact sync must be bit-identical)
#   make bench      full experiment matrix -> results/
#   make clean
#
# Compiler selection:
#   Linux:  gcc -fopenmp
#   macOS:  Apple clang + Homebrew libomp  (brew install libomp)
#           or Homebrew GCC:  make CC=gcc-15
# Override any time with  make CC=<compiler>.

UNAME := $(shell uname -s)

ifeq ($(origin CC),default)
  ifeq ($(UNAME),Darwin)
    CC := clang
  else
    CC := gcc
  endif
endif

# Apple clang needs libomp from Homebrew; GCC (and LLVM clang on Linux) use -fopenmp.
IS_APPLE_CLANG := $(shell $(CC) --version 2>/dev/null | grep -q "Apple clang" && echo 1)
ifeq ($(IS_APPLE_CLANG),1)
  LIBOMP   ?= $(shell brew --prefix libomp 2>/dev/null)
  OMPFLAGS := -Xpreprocessor -fopenmp -I$(LIBOMP)/include
  OMPLIBS  := -L$(LIBOMP)/lib -lomp
else
  OMPFLAGS := -fopenmp
  OMPLIBS  := -fopenmp
endif

# -ffp-contract=off: no fused multiply-add contraction, so the same source
# gives the same floating-point results wherever it is inlined. Together with
# no -ffast-math this keeps sequential and parallel runs bit-identical.
CFLAGS  ?= -O2
CFLAGS  += -std=c11 -D_POSIX_C_SOURCE=200809L -Wall -Wextra -ffp-contract=off $(OMPFLAGS)
LDLIBS  += $(OMPLIBS) -lm

SRC := $(wildcard src/*.c)
OBJ := $(SRC:src/%.c=build/%.o)
BIN := build/platoon

.PHONY: all test bench clean print-cc

all: $(BIN)

$(BIN): $(OBJ)
	$(CC) $(CFLAGS) -o $@ $^ $(LDLIBS)

build/%.o: src/%.c $(wildcard src/*.h) | build
	$(CC) $(CFLAGS) -c -o $@ $<

build:
	mkdir -p build

test: $(BIN)
	./$(BIN) selftest

bench: $(BIN)
	./scripts/run_experiments.sh

clean:
	rm -rf build

print-cc:
	@echo $(CC)
