# CICE-evp1d

The cice-kernel project host files that have been used during the EVP refactorization and performance study conducted by Intel. The underlying refactorization idea of transforming the data representation from 2D into 1D representations were originally developed at the Danish Meteorological Institute and these transformation ideas were tested in an earlier version of the CICE code. This work builds on and extends the outcome of these tests. 

All input data required for running these benchmarks have been generated and provided by DMI in 2022. The input data is based both on the current operational setup at DMI and a setup based on the Regional Arctic Model system (RASM). The variations of the EVP kernel implementations are all based on the upstream CICE implementation which is extracted in the file (src/bench\_v0.F90).
This repository contains the files needed for standalone tests of the 1d evp solver implemented in CICE.

# Source code

The source files are split into three groups:

- *wrapper files* that ensure declaration of arrays, initialization, IO handling, timers, OpenMP init etc. This is write-once code solely useful for supporting current driver files with sufficient infrastructure.

- *driver files* that are used to run the EVP benchmarks, i.e. they call code that will ensure proper runtime initialization based on the namelist settings and they have timers around the benchmark and exposes how EVP is supposed to be called. There is one driver file per EVP implementation and some driver files include both a CPU and a GPU version. The driver files all have the prefix bench in their filename.

- *evp files* these are the files that contain the EVP code. They are all based on the upstream version found in bench-v0.F90. There is one file for each EVP implementation and some of the files contain both a CPU and a GPU version. The name of these files matches the name of the driver files.

# Building

```
Link the correct compiler options to compiler.mk
Currently two sets are provided
compiler_Freya.mk (hpc at DMI). Note that version 2e do not compile., Till Rasmussen, DMI.
compiler_oneapi.mk (original) provided by Intel
make
ls bin/${COMPILER}/v*/*/evp_v? bin/${COMPILER}/v*/*/evp_v?? will list all the binaries that have been generated.
```

# Running

It requires input files in order to run this code. Please consult DMI in case you are interested in running the kernel.
As for the namelist parameters there are a couple of them that are tied to the testcase at hand, namely:

```
&kernel_nml
 ndte      = 1000           ! number of outer convergence iterations (typical value is 1000)
 testscale = 1              ! allowing one to weak-scale the given testcase, v0 is 4x increase for each +1 in testscale, v1-v2 is 2x increase for each +1 in testscale.
 binoutput = .true.         ! whether or not one wishes to dump results emerging from stress and stepu, respectively
 lindividual = .true.       ! whether or not one wishes to run both stress and stepu alone (i.e. withour barriers between) on top of running EVP to get individual measurements of the two stages too.
/

```

# Testing for correctness

The baseline version (v0) does not support SIMD instructions and we cannot expect binary identical results once we start using these. Thus, a necessary condition for retaining binary identical results is not to use SIMD instructions. This is also reflected in the md5sums for f2 below whereas f0 (serial) and f1 (multicore but no SIMD) retain binary identical results:

```
for f in f0 f1 f2
do
  for vv in 0 1 2a 2b 2c 2d 2e
  do
    ~/check-out/cice-evp1d/bin/ifort/v${vv}/${f}/evp_v${vv} > v${vv}_${f}.log
    thread=$(grep "Running " v${vv}_${f}.log)
    echo "Testing results for flag $f for v${vv} : ${thread}"
  done
  md5sum output_stress_v*1d.bin
  md5sum output_stepu_v*1d.bin
done
```

```
Testing results for flag f0 for v0 :
Testing results for flag f0 for v1 : Running without openMP
Testing results for flag f0 for v2a : Running without openMP
Testing results for flag f0 for v2b : Running without openMP
Testing results for flag f0 for v2c : Running without openMP
Testing results for flag f0 for v2d : Running without openMP
Testing results for flag f0 for v2e : Running without openMP
cb12f98a34640d65ece95770fa3e7d68  output_stress_v0_1d.bin
cb12f98a34640d65ece95770fa3e7d68  output_stress_v1_1d.bin
cb12f98a34640d65ece95770fa3e7d68  output_stress_v2a_1d.bin
cb12f98a34640d65ece95770fa3e7d68  output_stress_v2b_1d.bin
cb12f98a34640d65ece95770fa3e7d68  output_stress_v2c_1d.bin
cb12f98a34640d65ece95770fa3e7d68  output_stress_v2d_1d.bin
cb12f98a34640d65ece95770fa3e7d68  output_stress_v2e_1d.bin
3004476cb73ceed3dd0797a503b2f69c  output_stepu_v0_1d.bin
3004476cb73ceed3dd0797a503b2f69c  output_stepu_v1_1d.bin
3004476cb73ceed3dd0797a503b2f69c  output_stepu_v2a_1d.bin
3004476cb73ceed3dd0797a503b2f69c  output_stepu_v2b_1d.bin
3004476cb73ceed3dd0797a503b2f69c  output_stepu_v2c_1d.bin
3004476cb73ceed3dd0797a503b2f69c  output_stepu_v2d_1d.bin
3004476cb73ceed3dd0797a503b2f69c  output_stepu_v2e_1d.bin
Testing results for flag f1 for v0 :
Testing results for flag f1 for v1 : Running openMP with     8 threads
Testing results for flag f1 for v2a : Running openMP with     8 threads
Testing results for flag f1 for v2b : Running openMP with     8 threads
Testing results for flag f1 for v2c : Running openMP with     8 threads
Testing results for flag f1 for v2d : Running openMP with     8 threads
Testing results for flag f1 for v2e : Running openMP with     8 threads
cb12f98a34640d65ece95770fa3e7d68  output_stress_v0_1d.bin
cb12f98a34640d65ece95770fa3e7d68  output_stress_v1_1d.bin
cb12f98a34640d65ece95770fa3e7d68  output_stress_v2a_1d.bin
cb12f98a34640d65ece95770fa3e7d68  output_stress_v2b_1d.bin
cb12f98a34640d65ece95770fa3e7d68  output_stress_v2c_1d.bin
cb12f98a34640d65ece95770fa3e7d68  output_stress_v2d_1d.bin
cb12f98a34640d65ece95770fa3e7d68  output_stress_v2e_1d.bin
3004476cb73ceed3dd0797a503b2f69c  output_stepu_v0_1d.bin
3004476cb73ceed3dd0797a503b2f69c  output_stepu_v1_1d.bin
3004476cb73ceed3dd0797a503b2f69c  output_stepu_v2a_1d.bin
3004476cb73ceed3dd0797a503b2f69c  output_stepu_v2b_1d.bin
3004476cb73ceed3dd0797a503b2f69c  output_stepu_v2c_1d.bin
3004476cb73ceed3dd0797a503b2f69c  output_stepu_v2d_1d.bin
3004476cb73ceed3dd0797a503b2f69c  output_stepu_v2e_1d.bin
Testing results for flag f2 for v0 :
Testing results for flag f2 for v1 : Running openMP with     8 threads
Testing results for flag f2 for v2a : Running openMP with     8 threads
Testing results for flag f2 for v2b : Running openMP with     8 threads
Testing results for flag f2 for v2c : Running openMP with     8 threads
Testing results for flag f2 for v2d : Running openMP with     8 threads
Testing results for flag f2 for v2e : Running openMP with     8 threads
16d6ae3b093e7ef4d6d3dcea2984def8  output_stress_v0_1d.bin
2b6b8e4a297dff9d2e457de158442ba9  output_stress_v1_1d.bin
0fdc28b5ab574a19d41a43c4dc8542fe  output_stress_v2a_1d.bin
0fdc28b5ab574a19d41a43c4dc8542fe  output_stress_v2b_1d.bin
0fdc28b5ab574a19d41a43c4dc8542fe  output_stress_v2c_1d.bin
ff333fe6c22cf4a0d0124c3cd494280f  output_stress_v2d_1d.bin
0fdc28b5ab574a19d41a43c4dc8542fe  output_stress_v2e_1d.bin
848ed7323e1df9cc4b2149362cd0a92a  output_stepu_v0_1d.bin
98731a62d9db1306b3f94f9b21caf3d0  output_stepu_v1_1d.bin
86ef0a968dfad37ae3c5713f6ad89de3  output_stepu_v2a_1d.bin
43d22ac175030c2958688900b8faf5ce  output_stepu_v2b_1d.bin
43d22ac175030c2958688900b8faf5ce  output_stepu_v2c_1d.bin
e795c11d6841e41cb888ad9f5784efcf  output_stepu_v2d_1d.bin
323ae6467f6665f5ba8af944e5781ffd  output_stepu_v2e_1d.bin
```

