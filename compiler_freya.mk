# ===============================================================================
# Copyright (C) 2023, Intel Corporation
# Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
# 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
# 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
# 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS “AS IS” AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# ===============================================================================
COMPILER=ftn
CC=CC
FCFLAG=-cpp -module
CFLAGS=-diag-disable=10441
LDFLAGS=

# CPU flags
# f0 (serial) and f1 (multicore but NO SIMD) are used to prove binary identical results across all variants and upstream runs, f2 focus on performance
FCf0:=-g -O0 -traceback -fp-model precise -fpe0 -FR -CB -warn all
FCf1:=-g -O1 -qopenmp -traceback -fp-model precise -fpe0 -FR -CB -no-vec -no-fma -warn all
FCf2:=-qopenmp -O3 -xCORE-AVX512 -qopt-zmm-usage=high -ipo -finline-functions -finline -parallel

# GPU flags (requires COMPILER=ifx, CC=icx)
FCg1:=-D_OPENMP_TARGET -fiopenmp -O2 -fopenmp-targets=spir64="-mllvm -vpo-paropt-enable-64bit-opencl-atomics=true"
FCg1f2018:=-D_OPENMP_TARGET -fiopenmp -O2 -fopenmp-targets=spir64="-mllvm -vpo-paropt-enable-64bit-opencl-atomics=true" -fopenmp-target-do-concurrent

# Common
flags = ${checkflags} f2 # g1 g1f2018
