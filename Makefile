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
versions  = v0 v1 v2a v2b v2c v2d v2e
ompver    = v1 v2a v2b v2c v2d v2e

# all binaries
bin_v0       = $(foreach F, $(flags), $(BINDIR)/$(COMPILER)/v0/$(F)/evp_v0)
bin_v1       = $(foreach F, $(flags), $(BINDIR)/$(COMPILER)/v1/$(F)/evp_v1)
bin_v2a       = $(foreach F, $(flags), $(BINDIR)/$(COMPILER)/v2a/$(F)/evp_v2a)
bin_v2b       = $(foreach F, $(flags), $(BINDIR)/$(COMPILER)/v2b/$(F)/evp_v2b)
bin_v2c       = $(foreach F, $(flags), $(BINDIR)/$(COMPILER)/v2c/$(F)/evp_v2c)
bin_v2d       = $(foreach F, $(flags), $(BINDIR)/$(COMPILER)/v2d/$(F)/evp_v2d)
bin_v2e       = $(foreach F, $(flags), $(BINDIR)/$(COMPILER)/v2e/$(F)/evp_v2e)
allbin       = $(bin_v0) $(bin_v1) $(bin_v2a) $(bin_v2b) $(bin_v2c) $(bin_v2d) #$(bin_v2e)

# objects
objtimer     = $(foreach v, $(versions), $(foreach F, $(flags), $(BINDIR)/$(COMPILER)/$(v)/$(F)/timer.o))
objice       = $(foreach v, $(versions), $(foreach F, $(flags), $(BINDIR)/$(COMPILER)/$(v)/$(F)/ice.o))
objvars      = $(foreach v, $(versions), $(foreach F, $(flags), $(BINDIR)/$(COMPILER)/$(v)/$(F)/vars.o))
objsde       = $(foreach v, $(versions), $(foreach F, $(flags), $(BINDIR)/$(COMPILER)/$(v)/$(F)/api_fortran_sde.o))
objdomp      = $(foreach v, $(ompver),   $(foreach F, $(flags), $(BINDIR)/$(COMPILER)/$(v)/$(F)/domp.o))
objnuma      = $(foreach v, $(ompver),   $(foreach F, $(flags), $(BINDIR)/$(COMPILER)/$(v)/$(F)/numa.o))
objbench_v0  = $(foreach F, $(flags), $(BINDIR)/$(COMPILER)/v0/$(F)/bench_v0.o)
objevp_v0    = $(foreach F, $(flags), $(BINDIR)/$(COMPILER)/v0/$(F)/evp_v0.o)
objbench_v1  = $(foreach F, $(flags), $(BINDIR)/$(COMPILER)/v1/$(F)/bench_v1.o)
objevp_v1    = $(foreach F, $(flags), $(BINDIR)/$(COMPILER)/v1/$(F)/evp_v1.o)
objbench_v2a = $(foreach F, $(flags), $(BINDIR)/$(COMPILER)/v2a/$(F)/bench_v2a.o)
objevp_v2a   = $(foreach F, $(flags), $(BINDIR)/$(COMPILER)/v2a/$(F)/evp_v2a.o)
objbench_v2b = $(foreach F, $(flags), $(BINDIR)/$(COMPILER)/v2b/$(F)/bench_v2b.o)
objevp_v2b   = $(foreach F, $(flags), $(BINDIR)/$(COMPILER)/v2b/$(F)/evp_v2b.o)
objbench_v2c = $(foreach F, $(flags), $(BINDIR)/$(COMPILER)/v2c/$(F)/bench_v2c.o)
objevp_v2c   = $(foreach F, $(flags), $(BINDIR)/$(COMPILER)/v2c/$(F)/evp_v2c.o)
objbench_v2d = $(foreach F, $(flags), $(BINDIR)/$(COMPILER)/v2d/$(F)/bench_v2d.o)
objevp_v2d   = $(foreach F, $(flags), $(BINDIR)/$(COMPILER)/v2d/$(F)/evp_v2d.o)
objbench_v2e = $(foreach F, $(flags), $(BINDIR)/$(COMPILER)/v2e/$(F)/bench_v2e.o)
objevp_v2e   = $(foreach F, $(flags), $(BINDIR)/$(COMPILER)/v2e/$(F)/evp_v2e.o)

# targets
all: $(allbin)

# dependencies
$(objtimer)    : $(objice)
$(objdomp)     : $(objice)
$(objvars)     : $(objice) $(objdomp)
$(objbench_v0) : $(objice) $(objvars)
$(objnuma)     : $(objice) $(objdomp) $(objvars)
#$(objevp_v0)   : $(objice) $(objtimer) $(objvars) $(objbench_v0) $(objsde)
$(objevp_v0)   : $(objice) $(objtimer) $(objvars) $(objbench_v0)
$(objbench_v1) : $(objice) $(objvars) $(objdomp)
$(objevp_v1)   : $(objice) $(objtimer) $(objvars) $(objdomp) $(objbench_v1) $(objnuma)
$(objbench_v2a) : $(objice) $(objvars) $(objdomp)
$(objevp_v2a)   : $(objice) $(objtimer) $(objvars) $(objdomp) $(objbench_v2a) $(objnuma)
#$(objbench_v2b) : $(objice) $(objvars) $(objdomp)
$(objbench_v2b) : $(objice) $(objdomp)
#$(objevp_v2b)   : $(objice) $(objtimer) $(objvars) $(objdomp) $(objbench_v2b) $(objnuma) $(objsde)
$(objevp_v2b)   : $(objice) $(objtimer) $(objvars) $(objdomp) $(objbench_v2b) $(objnuma)
$(objbench_v2c) : $(objice) $(objvars) $(objdomp)
$(objevp_v2c)   : $(objice) $(objtimer) $(objvars) $(objdomp) $(objbench_v2c) $(objnuma)
$(objbench_v2d) : $(objice) $(objvars) $(objdomp)
$(objevp_v2d)   : $(objice) $(objtimer) $(objvars) $(objdomp) $(objbench_v2d) $(objnuma)
$(objbench_v2e) : $(objice) $(objvars) $(objdomp)
$(objevp_v2e)   : $(objice) $(objtimer) $(objvars) $(objdomp) $(objbench_v2e) $(objnuma)

clean:
	@rm -fr $(BINDIR)/$(COMPILER)
	@rm -f src/*.mod

realclean:
	@rm -fr $(BINDIR)

$(objsde): src/api_fortran_sde.c
	test -d $(@D) || mkdir -p $(@D)
	$(CC) -c $(CFLAGS) $< -o $@

$(objtimer): src/timer.F90
	test -d $(@D) || mkdir -p $(@D)
	$(COMPILER) -D$(word 3, $(subst 2e,2,$(subst 2d,2,$(subst 2c,2,$(subst 2b,2,$(subst 2a,2,$(subst /, ,$(@D)))))))) $(FCFLAG) $(@D) $(FC$(notdir $(@D))) -c $< -o $@

$(objice): src/ice.F90
	test -d $(@D) || mkdir -p $(@D)
	$(COMPILER) -D$(word 3, $(subst 2e,2,$(subst 2d,2,$(subst 2c,2,$(subst 2b,2,$(subst 2a,2,$(subst /, ,$(@D)))))))) $(FCFLAG) $(@D) $(FC$(notdir $(@D))) -c $< -o $@

$(objvars): src/vars.F90
	test -d $(@D) || mkdir -p $(@D)
	$(COMPILER) -D$(word 3, $(subst 2e,2,$(subst 2d,2,$(subst 2c,2,$(subst 2b,2,$(subst 2a,2,$(subst /, ,$(@D)))))))) $(FCFLAG) $(@D) $(FC$(notdir $(@D))) -c $< -o $@

$(objdomp): src/domp.F90
	test -d $(@D) || mkdir -p $(@D)
	$(COMPILER) -D$(word 3, $(subst 2e,2,$(subst 2d,2,$(subst 2c,2,$(subst 2b,2,$(subst 2a,2,$(subst /, ,$(@D)))))))) $(FCFLAG) $(@D) $(FC$(notdir $(@D))) -c $< -o $@

$(objnuma): src/numa.F90
	test -d $(@D) || mkdir -p $(@D)
	$(COMPILER) -D$(word 3, $(subst 2e,2,$(subst 2d,2,$(subst 2c,2,$(subst 2b,2,$(subst 2a,2,$(subst /, ,$(@D)))))))) $(FCFLAG) $(@D) $(FC$(notdir $(@D))) -c $< -o $@

$(objbench_v0): src/bench_v0.F90
	test -d $(@D) || mkdir -p $(@D)
	$(COMPILER) -D$(word 3, $(subst 2e,2,$(subst 2d,2,$(subst 2c,2,$(subst 2b,2,$(subst 2a,2,$(subst /, ,$(@D)))))))) $(FCFLAG) $(@D) $(FC$(notdir $(@D))) -c $< -o $@

$(objevp_v0): src/evp_v0.F90
	test -d $(@D) || mkdir -p $(@D)
	$(COMPILER) -D$(word 3, $(subst 2e,2,$(subst 2d,2,$(subst 2c,2,$(subst 2b,2,$(subst 2a,2,$(subst /, ,$(@D)))))))) $(FCFLAG) $(@D) $(FC$(notdir $(@D))) -c $< -o $@

$(objbench_v1): src/bench_v1.F90
	test -d $(@D) || mkdir -p $(@D)
	$(COMPILER) -D$(word 3, $(subst 2e,2,$(subst 2d,2,$(subst 2c,2,$(subst 2b,2,$(subst 2a,2,$(subst /, ,$(@D)))))))) $(FCFLAG) $(@D) $(FC$(notdir $(@D))) -c $< -o $@

$(objevp_v1): src/evp_v1.F90
	test -d $(@D) || mkdir -p $(@D)
	$(COMPILER) -D$(word 3, $(subst 2e,2,$(subst 2d,2,$(subst 2c,2,$(subst 2b,2,$(subst 2a,2,$(subst /, ,$(@D)))))))) $(FCFLAG) $(@D) $(FC$(notdir $(@D))) -c $< -o $@

$(objbench_v2a): src/bench_v2a.F90
	test -d $(@D) || mkdir -p $(@D)
	$(COMPILER) -D$(word 3, $(subst 2e,2,$(subst 2d,2,$(subst 2c,2,$(subst 2b,2,$(subst 2a,2,$(subst /, ,$(@D)))))))) $(FCFLAG) $(@D) $(FC$(notdir $(@D))) -c $< -o $@

$(objevp_v2a): src/evp_v2a.F90
	test -d $(@D) || mkdir -p $(@D)
	$(COMPILER) -D$(word 3, $(subst 2e,2,$(subst 2d,2,$(subst 2c,2,$(subst 2b,2,$(subst 2a,2,$(subst /, ,$(@D)))))))) $(FCFLAG) $(@D) $(FC$(notdir $(@D))) -c $< -o $@

$(objbench_v2b): src/bench_v2b.F90
	test -d $(@D) || mkdir -p $(@D)
	$(COMPILER) -D$(word 3, $(subst 2e,2,$(subst 2d,2,$(subst 2c,2,$(subst 2b,2,$(subst 2a,2,$(subst /, ,$(@D)))))))) $(FCFLAG) $(@D) $(FC$(notdir $(@D))) -c $< -o $@

$(objevp_v2b): src/evp_v2b.F90
	test -d $(@D) || mkdir -p $(@D)
	$(COMPILER) -D$(word 3, $(subst 2e,2,$(subst 2d,2,$(subst 2c,2,$(subst 2b,2,$(subst 2a,2,$(subst /, ,$(@D)))))))) $(FCFLAG) $(@D) $(FC$(notdir $(@D))) -c $< -o $@

$(objbench_v2c): src/bench_v2c.F90
	test -d $(@D) || mkdir -p $(@D)
	$(COMPILER) -D$(word 3, $(subst 2e,2,$(subst 2d,2,$(subst 2c,2,$(subst 2b,2,$(subst 2a,2,$(subst /, ,$(@D)))))))) $(FCFLAG) $(@D) $(FC$(notdir $(@D))) -c $< -o $@

$(objevp_v2c): src/evp_v2c.F90
	test -d $(@D) || mkdir -p $(@D)
	$(COMPILER) -D$(word 3, $(subst 2e,2,$(subst 2d,2,$(subst 2c,2,$(subst 2b,2,$(subst 2a,2,$(subst /, ,$(@D)))))))) $(FCFLAG) $(@D) $(FC$(notdir $(@D))) -c $< -o $@

$(objbench_v2d): src/bench_v2d.F90
	test -d $(@D) || mkdir -p $(@D)
	$(COMPILER) -D$(word 3, $(subst 2e,2,$(subst 2d,2,$(subst 2c,2,$(subst 2b,2,$(subst 2a,2,$(subst /, ,$(@D)))))))) $(FCFLAG) $(@D) $(FC$(notdir $(@D))) -c $< -o $@

$(objevp_v2d): src/evp_v2d.F90
	test -d $(@D) || mkdir -p $(@D)
	$(COMPILER) -D$(word 3, $(subst 2e,2,$(subst 2d,2,$(subst 2c,2,$(subst 2b,2,$(subst 2a,2,$(subst /, ,$(@D)))))))) $(FCFLAG) $(@D) $(FC$(notdir $(@D))) -c $< -o $@

$(objbench_v2e): src/bench_v2e.F90
	test -d $(@D) || mkdir -p $(@D)
	$(COMPILER) -D$(word 3, $(subst 2e,2,$(subst 2d,2,$(subst 2c,2,$(subst 2b,2,$(subst 2a,2,$(subst /, ,$(@D)))))))) $(FCFLAG) $(@D) $(FC$(notdir $(@D))) -c $< -o $@

$(objevp_v2e): src/evp_v2e.F90
	test -d $(@D) || mkdir -p $(@D)
	$(COMPILER) -D$(word 3, $(subst 2e,2,$(subst 2d,2,$(subst 2c,2,$(subst 2b,2,$(subst 2a,2,$(subst /, ,$(@D)))))))) $(FCFLAG) $(@D) $(FC$(notdir $(@D))) -c $< -o $@

$(bin_v0): $(objice) $(objtimer) $(objvars) $(objbench_v0) $(objevp_v0)
	test -d $(@D) || mkdir -p $(@D)
	$(COMPILER) $(FCFLAG) $(@D) $(FC$(notdir $(@D))) -o $@ $(@D)/*.o $(LDFLAGS)

$(bin_v1): $(objice) $(objtimer) $(objvars) $(objdomp) $(objbench_v1) $(objevp_v1) $(objnuma)
	test -d $(@D) || mkdir -p $(@D)
	$(COMPILER) $(FCFLAG) $(@D) $(FC$(notdir $(@D))) -o $@ $(@D)/*.o $(LDFLAGS)

$(bin_v2a): $(objice) $(objtimer) $(objvars) $(objdomp) $(objbench_v2a) $(objevp_v2a) $(objnuma)
	test -d $(@D) || mkdir -p $(@D)
	$(COMPILER) $(FCFLAG) $(@D) $(FC$(notdir $(@D))) -o $@ $(@D)/*.o $(LDFLAGS)

$(bin_v2b): $(objice) $(objtimer) $(objvars) $(objdomp) $(objbench_v2b) $(objevp_v2b) $(objnuma)
	test -d $(@D) || mkdir -p $(@D)
	$(COMPILER) $(FCFLAG) $(@D) $(FC$(notdir $(@D))) -o $@ $(@D)/*.o $(LDFLAGS)

$(bin_v2c): $(objice) $(objtimer) $(objvars) $(objdomp) $(objbench_v2c) $(objevp_v2c) $(objnuma)
	test -d $(@D) || mkdir -p $(@D)
	$(COMPILER) $(FCFLAG) $(@D) $(FC$(notdir $(@D))) -o $@ $(@D)/*.o $(LDFLAGS)

$(bin_v2d): $(objice) $(objtimer) $(objvars) $(objdomp) $(objbench_v2d) $(objevp_v2d) $(objnuma)
	test -d $(@D) || mkdir -p $(@D)
	$(COMPILER) $(FCFLAG) $(@D) $(FC$(notdir $(@D))) -o $@ $(@D)/*.o $(LDFLAGS)

$(bin_v2e): $(objice) $(objtimer) $(objvars) $(objdomp) $(objbench_v2e) $(objevp_v2e) $(objnuma)
	test -d $(@D) || mkdir -p $(@D)
	$(COMPILER) $(FCFLAG) $(@D) $(FC$(notdir $(@D))) -o $@ $(@D)/*.o $(LDFLAGS)
