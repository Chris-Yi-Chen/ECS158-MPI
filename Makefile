# Get all C files in the directory
src := $(sort $(wildcard *.c))

bin := $(src:.c=)

CC := mpicc
CFLAGS := -Wall -Werror -O2

all: $(bin)

# use '-' because some files might be meant to fail
%: %.c
	-$(CC) $(CFLAGS) -o $@ $<

clean:
	-rm -f $(bin)

