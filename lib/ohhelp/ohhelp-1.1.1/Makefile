#  File: samplef.mk
#  Version 1.1.1 (2015/10/23)
#  Copyright (C) 2009-2015  Hiroshi Nakashima <h.nakashima@media.kyoto-u.ac.jp>
#                           (ACCMS, Kyoto University)
#  This program can be freely used, redistributed and modified for non-
#  commercial purpose providing that the copyright notice above remains
#  unchanged.

PROGRAM = libohhelp.a

LINKER          = $(FC)

COMMONHDRS      = oh_config.h oh_stats.h
OHHDRS          = ohhelp1.h ohhelp2.h ohhelp3.h oh_part.h
OHOBJS          = ohhelp1.o ohhelp2.o ohhelp3.o

FHDRS           = ohhelp_f.h
FMODS           = oh_type.o oh_mod1.o oh_mod2.o oh_mod3.o oh_mod4p.o

FOBJS           = 

OBJS            = $(FOBJS) $(FMODS) $(OHOBJS)

FFLAGS := $(FFLAGS)
CFLAGS := $(CFLAGS)

$(PROGRAM): $(OBJS)
	ar rc $@ $(OBJS)

$(FOBJS):%.o:   %.F90 $(FMODS) $(FHDRS) $(COMMONHDRS)
	$(FC) $(FFLAGS) -c $< -o $@

$(FMODS):%.o:   %.F90 $(COMMONHDRS)
	$(FC) $(CFLAGS) -c $< -o $@

oh_mod1.o: oh_type.o
oh_mod2.o: oh_mod1.o
oh_mod3.o: oh_mod2.o
oh_mod4p.o: oh_mod3.o

$(OHOBJS):%.o:  %.c $(COMMONHDRS) $(OHHDRS)
	$(CC) $(CFLAGS) -c $< -o $@

clean:;
	-rm *.o *.mod
	-rm $(PROGRAM)
