#	Makefile
#	2021-01-30	Markku-Juhani O. Saarinen <mjos@pqshield.com>
#	Copyright (c) 2021, PQShield Ltd.  All rights reserved.

export

CC 		=	gcc
CFLAGS	=	-Wall -Wextra -O3
#-fsanitize=address,undefined -O2 -g
LDLIBS	=	-lfftw3 -lm

XBIN	=	xtest
CSRC	=	$(wildcard *.c)
SSRC	=	$(wildcard *.S)
OBJS	=	$(CSRC:.c=.o) $(SSRC:.S=.o)

$(XBIN):	$(OBJS)
			$(CC) $(CFLAGS) -o $(XBIN) $(OBJS) $(LDLIBS)

%.o:		%.[cS]
			$(CC) $(CFLAGS) -c $^ -o $@

clean:
			rm -f $(OBJS) $(XBIN)
