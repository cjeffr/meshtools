# Makefile for gridscope C translation

CC = gcc
CFLAGS = -Wall -Wextra -std=c99 -O2 -g
LDFLAGS = -lm

# Source files
SOURCES = gridscope.c
HEADERS = gridscope.h
TEST_SOURCES = test_gridscope.c
OBJECTS = $(SOURCES:.c=.o)
TEST_OBJECTS = $(TEST_SOURCES:.c=.o)

# Targets
TARGET = gridscope
TEST_TARGET = test_gridscope

# Default target
all: $(TARGET) $(TEST_TARGET)

# Build main executable
$(TARGET): $(OBJECTS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

# Build test executable
$(TEST_TARGET): $(TEST_OBJECTS) gridscope.o
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

# Object file compilation
%.o: %.c $(HEADERS)
	$(CC) $(CFLAGS) -c $< -o $@

# Run tests
test: $(TEST_TARGET)
	@echo "Running test suite..."
	@./$(TEST_TARGET)

# Clean up
clean:
	rm -f $(OBJECTS) $(TEST_OBJECTS) $(TARGET) $(TEST_TARGET)
	rm -f *.o *.grd *.14 .tmp01

# Install (if needed)
install: $(TARGET)
	cp $(TARGET) /usr/local/bin/

.PHONY: all test clean install
