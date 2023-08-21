# ===============================================================================
# Copyright (C) 2023, Intel Corporation
# Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
# 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
# 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
# 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS “AS IS” AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# ===============================================================================
SHELL=/bin/bash
.PHONY: all clean realclean
.SUFFIXES:
.SUFFIXES: .c .F90 .o

include compiler.mk

BINDIR ?= bin
BASEROOT=$(shell pwd)
versions  = v2b
ompver    = v2b

# all binaries
#bin_v2b       = $(foreach F, $(flags), $(BINDIR)/$(COMPILER)/v2b/$(F)/evp_v2b)
cice_evp1d_standalone    = $(foreach F, $(flags), $(BINDIR)/$(F)/standalone)

# objects from CICE
#ice_kinds_mod      = $(foreach v, $(versions), $(foreach F, $(flags), $(BINDIR)/$(COMPILER)/$(F)/ice_kinds_mod.o))
ice_kinds_mod      = $(foreach F, $(flags), $(BINDIR)/$(F)/ice_kinds_mod.o)
ice_constants      = $(foreach F, $(flags), $(BINDIR)/$(F)/ice_constants.o)
ice_gather_scatter = $(foreach F, $(flags), $(BINDIR)/$(F)/ice_gather_scatter.o)
ice_communicate    = $(foreach F, $(flags), $(BINDIR)/$(F)/ice_communicate.o)
ice_domain         = $(foreach F, $(flags), $(BINDIR)/$(F)/ice_domain.o)
ice_dyn_shared     = $(foreach F, $(flags), $(BINDIR)/$(F)/ice_dyn_shared.o)
# New to be implemented
ice_dyn_evp1d      = $(foreach F, $(flags), $(BINDIR)/$(F)/ice_dyn_evp1d.o)
#Temporary debug 
debug_evp1d        = $(foreach F, $(flags), $(BINDIR)/$(F)/debug_evp1d.o)
#objsde       = $(foreach v, $(versions), $(foreach F, $(flags), $(BINDIR)/$(COMPILER)/$(v)/$(F)/api_fortran_sde.o))
objmyomp           = $(foreach F, $(flags), $(BINDIR)/$(F)/myomp.o)
objnuma            = $(foreach F, $(flags), $(BINDIR)/$(F)/numa.o)
objbench_v2b       = $(foreach F, $(flags), $(BINDIR)/$(F)/bench_v2b.o)
objevp1d_int_v2b   = $(foreach F, $(flags), $(BINDIR)/$(F)/evp1d.o)
objstandalone      = $(foreach F, $(flags), $(BINDIR)/$(F)/ice_dyn_evp1d_standalone.o)
# targets
all: $(cice_evp1d_standalone)

# dependencies
$(objtimer)           : $(ice_kinds_mod)
$(ice_constants)      : $(ice_kinds_mod) 
$(ice_communicate)    : $(ice_kinds_mod)
$(ice_domain)         : $(ice_kinds_mod)
$(ice_dyn_shared)     : $(ice_kinds_mod) $(ice_constants)
$(ice_gather_scatter) : $(ice_kinds_mod)
$(ice_dyn_evp1d)      : $(ice_kinds_mod) $(ice_constants) $(debug_evp1d) $(objbench_v2b) $(ice_communicate) $(ice_gather_scatter) $(ice_domain) $(ice_dyn_shared)
$(debug_evp1d)        : $(ice_kinds_mod)
$(objmyomp)           : $(ice_kinds_mod) $(ice_constants)
$(objbench_v2b)       : $(ice_kinds_mod) $(ice_dyn_shared) $(ice_constant) $(objmyomp)
$(objstandalone)      : $(ice_kinds_mod) $(ice_dyn_evp1d) $(ice_dyn_shared)
clean:
	@rm -fr $(BINDIR)
	@rm -f src/*.mod

realclean:
	@rm -fr $(BINDIR)

#$(objsde): ../src/api_fortran_sde.c
#	test -d $(@D) || mkdir -p $(@D)
#	$(CC) -c $(CFLAGS) $< -o $@

$(ice_kinds_mod): ice_kinds_mod.F90
	test -d $(@D) || mkdir -p $(@D)
	$(COMPILER) -D$(word 3, $(subst 2e,2,$(subst 2d,2,$(subst 2c,2,$(subst 2b,2,$(subst 2a,2,$(subst /, ,$(@D)))))))) $(FCFLAG) $(@D) $(FC$(notdir $(@D))) -c $< -o $@

$(ice_constants): ice_constants.F90
	test -d $(@D) || mkdir -p $(@D)
	$(COMPILER) -D$(word 3, $(subst 2e,2,$(subst 2d,2,$(subst 2c,2,$(subst 2b,2,$(subst 2a,2,$(subst /, ,$(@D)))))))) $(FCFLAG) $(@D) $(FC$(notdir $(@D))) -c $< -o $@

$(ice_gather_scatter): ice_gather_scatter.F90
	test -d $(@D) || mkdir -p $(@D)
	$(COMPILER) -D$(word 3, $(subst 2e,2,$(subst 2d,2,$(subst 2c,2,$(subst 2b,2,$(subst 2a,2,$(subst /, ,$(@D)))))))) $(FCFLAG) $(@D) $(FC$(notdir $(@D))) -c $< -o $@

$(ice_communicate): ice_communicate.F90
	test -d $(@D) || mkdir -p $(@D)
	$(COMPILER) -D$(word 3, $(subst 2e,2,$(subst 2d,2,$(subst 2c,2,$(subst 2b,2,$(subst 2a,2,$(subst /, ,$(@D)))))))) $(FCFLAG) $(@D) $(FC$(notdir $(@D))) -c $< -o $@

$(ice_domain): ice_domain.F90
	test -d $(@D) || mkdir -p $(@D)
	$(COMPILER) -D$(word 3, $(subst 2e,2,$(subst 2d,2,$(subst 2c,2,$(subst 2b,2,$(subst 2a,2,$(subst /, ,$(@D)))))))) $(FCFLAG) $(@D) $(FC$(notdir $(@D))) -c $< -o $@

$(ice_dyn_shared): ice_dyn_shared.F90
	test -d $(@D) || mkdir -p $(@D)
	$(COMPILER) -D$(word 3, $(subst 2e,2,$(subst 2d,2,$(subst 2c,2,$(subst 2b,2,$(subst 2a,2,$(subst /, ,$(@D)))))))) $(FCFLAG) $(@D) $(FC$(notdir $(@D))) -c $< -o $@

$(debug_evp1d): debug_evp1d.F90
	test -d $(@D) || mkdir -p $(@D)
	$(COMPILER) -D$(word 3, $(subst 2e,2,$(subst 2d,2,$(subst 2c,2,$(subst 2b,2,$(subst 2a,2,$(subst /, ,$(@D)))))))) $(FCFLAG) $(@D) $(FC$(notdir $(@D))) -c $< -o $@

$(objmyomp): ../src/myomp.F90
	test -d $(@D) || mkdir -p $(@D)
	$(COMPILER) -D$(word 3, $(subst 2e,2,$(subst 2d,2,$(subst 2c,2,$(subst 2b,2,$(subst 2a,2,$(subst /, ,$(@D)))))))) $(FCFLAG) $(@D) $(FC$(notdir $(@D))) -c $< -o $@

$(objnuma): ../src/numa.F90
	test -d $(@D) || mkdir -p $(@D)
	$(COMPILER) -D$(word 3, $(subst 2e,2,$(subst 2d,2,$(subst 2c,2,$(subst 2b,2,$(subst 2a,2,$(subst /, ,$(@D)))))))) $(FCFLAG) $(@D) $(FC$(notdir $(@D))) -c $< -o $@

$(objbench_v2b): ../src/bench_v2b.F90
	test -d $(@D) || mkdir -p $(@D)
	$(COMPILER) -D$(word 3, $(subst 2e,2,$(subst 2d,2,$(subst 2c,2,$(subst 2b,2,$(subst 2a,2,$(subst /, ,$(@D)))))))) $(FCFLAG) $(@D) $(FC$(notdir $(@D))) -c $< -o $@

#$(objevp_v2b): ../src/evp_v2b.F90
#	test -d $(@D) || mkdir -p $(@D)
#	$(COMPILER) -D$(word 3, $(subst 2e,2,$(subst 2d,2,$(subst 2c,2,$(subst 2b,2,$(subst 2a,2,$(subst /, ,$(@D)))))))) $(FCFLAG) $(@D) $(FC$(notdir $(@D))) -c $< -o $@

$(ice_dyn_evp1d): ice_dyn_evp1d.F90
	test -d $(@D) || mkdir -p $(@D)
	$(COMPILER) -D$(word 3, $(subst 2e,2,$(subst 2d,2,$(subst 2c,2,$(subst 2b,2,$(subst 2a,2,$(subst /, ,$(@D)))))))) $(FCFLAG) $(@D) $(FC$(notdir $(@D))) -c $< -o $@

$(objstandalone): ice_dyn_evp1d_standalone.F90
	test -d $(@D) || mkdir -p $(@D)
	$(COMPILER) -D$(word 3, $(subst 2e,2,$(subst 2d,2,$(subst 2c,2,$(subst 2b,2,$(subst 2a,2,$(subst /, ,$(@D)))))))) $(FCFLAG) $(@D) $(FC$(notdir $(@D))) -c $< -o $@

$(cice_evp1d_standalone): $(ice_kinds_mod) $(ice_constants) $(objmyomp) $(ice_gather_scatter) $(ice_communicate) $(ice_dyn_shared) $(ice_dyn_evp) $(debug_evp1d) $(objstandalone)
	test -d $(@D) || mkdir -p $(@D)
	$(COMPILER) $(FCFLAG) $(@D) $(FC$(notdir $(@D))) -o $@ $(@D)/*.o $(LDFLAGS)
