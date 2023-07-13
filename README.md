# CICE-evp1d

The cice-kernel project host files that have been used during the EVP refactorization and performance study conducted by Intel. The underlying refactorization idea of transforming the data representation from 2D into 1D representations were originally developed at the Danish Meteorological Institute and these transformation ideas were tested in an earlier version of the CICE code. This work builds on and extends the outcome of these tests. 

All input data required for running these benchmarks have been generated and provided by DMI in 2022. The input data is based both on the current operational setup at DMI and a setup based on the Regional Arctic Model system (RASM). The variations of the EVP kernel implementations are all based on the upstream CICE implementation which is extracted in the file (src/bench\_v0.F90).
This repository contains the files needed for standalone tests of the 1d evp solver implemented in CICE.

# Building

```
make
ls bin/ifort/v*/*/evp_v?? # will list all the binaries that have been generated.
```

# Running

It requires input files in order to run this code. Please consult DMI in case you are interested in running the kernel.

# Testing for correctness

We expect binary identical results as long as we confine ourselves not to use SIMD operations. The baseline version (v0) does not support SIMD instructions and we cannot expect binary identical results once we start using SIMD instructions. This is also reflected in the md5sums for f2 below whereas f0 (serial) and f1 (multicore but no SIMD) retain binary identical results:

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

