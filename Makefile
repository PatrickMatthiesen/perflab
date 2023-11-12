# Student's Makefile for the CS:APP Performance Lab
TEAM = pmat
VERSION = 1
HANDINDIR = 1

CC = gcc
CFLAGS = -Wall -O2 -m32 #-mavx #-floop-block -floop-interchange -floop-strip-mine -floop-parallelize-all -ftree-vectorize -ftree-vectorizer-verbose=2 -fopenmp -fno-strict-aliasing -fno-omit-frame-pointer -g -std=gnu99
LIBS = -lm

OBJS = driver.o kernels.o fcyc.o clock.o

all: driver

driver: $(OBJS) fcyc.h clock.h defs.h config.h
	$(CC) $(CFLAGS) $(OBJS) $(LIBS) -o driver

handin:
	cp kernels.c $(HANDINDIR)/$(TEAM)-$(VERSION)-kernels.c

clean: 
	-rm -f $(OBJS) driver core *~ *.o


