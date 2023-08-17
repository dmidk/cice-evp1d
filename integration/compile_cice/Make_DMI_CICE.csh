#!/bin/csh -f

# Based on CICE compile script.
#if (! $?CLEANUP) set CLEANUP = 0
set CLEANUP = 0

# compiledir / sourcecode
setenv ICE_BUILDDIR ${cwd}    # compiledir
# set to location of cice code
setenv ICE_SANDBOX /lustre/research/tar/test_evp/test_implementation/cice_main/CICE
# default commandline input
set compiler = intel
set ICE_BLDDEBUG = false

# Usage / input args
if ( `echo $argv | grep -c "\-h"` ) then
  echo "Usage: $0 {compiler} {DEBUG}"
  echo "\nDefault: $0 ${compiler} ${ICE_BLDDEBUG}"
  exit 2
endif

if ($#argv > 0) set compiler = $1
if ($#argv > 1) set ICE_BLDDEBUG = $2

setenv ICE_BLDDEBUG $ICE_BLDDEBUG

#### USER INPUT
setenv ICE_THREADED   true    # use openMP otherwise false
setenv ICE_IOTYPE     netcdf  # set to none if netcdf library is unavailable
setenv ICE_COMMDIR    serial     # mpi or serial
setenv ICE_CONSTOPT   cice
setenv ICE_MACHINE_BLDTHRDS 1 # 1/2
#setenv ICE_BLDDEBUG   true   # true: build debug flags
#### END USER UNPUT

#---------------------------------------------
#-- Machine specific / compiler
#---------------------------------------------
echo "_____________________"
echo " ${compiler}"
echo "_____________________"

# start modules:
source /opt/modules/default/init/tcsh

module rm PrgEnv-intel
module rm PrgEnv-gnu
module rm PrgEnv-cray
switch ( ${compiler} )
  case intel:
  case gnu:
  case cray:
    breaksw
  default:
    echo "ERROR: Compiler not supported : "${compiler}
    exit 1
    breaksw
endsw
module add PrgEnv-${compiler}
setenv ARCH DMI.${compiler}

module add cray-netcdf
echo "-----------------------"
echo "COMPILER: ${compiler}"
echo "MACHINE:  `uname -mn`"
echo "ARCH:     ${ARCH}"
echo "-----------------------"

module list
#---------------------------------------------
#-- Settings / CPPFLAGS etc.
#---------------------------------------------
setenv ICE_DRVOPT cice
setenv ICE_DRVDIR standalone
if (${ICE_BLDDEBUG} == 'true') then
   setenv ICE_OBJDIR  ${ICE_BUILDDIR}/compile/${compiler}_${ICE_DRVDIR}_${ICE_DRVOPT}_debug/
else
   setenv ICE_OBJDIR  ${ICE_BUILDDIR}/compile/${compiler}_${ICE_DRVDIR}_${ICE_DRVOPT}/
endif
cp Macros.${ARCH} ${ICE_SANDBOX}/
setenv MACFILE ${ICE_SANDBOX}/Macros.${ARCH}
if (! -e ${MACFILE}) then
  echo "ERROR: Macro file not found: "${MACFILE}
  exit 1
endif

# make log-file. ONLY stdout if empty
setenv ICE_BLDLOG_FILE
setenv ICE_BLDLOG_FILE cice.${compiler}.buildlog

#-----------------
# cpp flags for compilation
#DMI TEST and DMI_v1 are related to output for cases
# vert_therm_bug fixes a division by small numbers when transfer from ice/snow volume to thickness 
#setenv ICE_CPPDEFS "-DFORTRANUNDERSCORE -DDMI_TEST -DDMI_V1 -Dvert_therm_bug"
setenv ICE_CPPDEFS "-DFORTRANUNDERSCORE -Dvert_therm_bug -Dintegrate"
#setenv ICE_CPPDEFS "-DFORTRANUNDERSCORE -Dvert_therm_bug"
# -DDMI_TEST -DDMI_V1"
if (${ICE_IOTYPE} == 'netcdf') then
  set IODIR = io_netcdf
  setenv ICE_CPPDEFS "${ICE_CPPDEFS} -DUSE_NETCDF"
else if (${ICE_IOTYPE} == 'pio') then
  set IODIR = io_pio
  setenv ICE_CPPDEFS "${ICE_CPPDEFS} -DUSE_NETCDF"
else if (${ICE_IOTYPE} == 'netcdf_bin') then
  set IODIR = io_netcdf_bin
  setenv ICE_CPPDEFS "${ICE_CPPDEFS} -DUSE_NETCDF"
else
  set IODIR = io_binary
endif

echo "cleaning objdir"
rm -r -f ${ICE_OBJDIR}
mkdir -p ${ICE_OBJDIR}
cd ${ICE_OBJDIR}
echo "Goto objdir: "${ICE_OBJDIR}
### List of source code directories (in order of importance).
#${ICE_SANDBOX}/cicecore/drivers/${ICE_DRVDIR}/${ICE_DRVOPT} later versions than 6.0
cat >! Filepath << EOF
${ICE_SANDBOX}/cicecore/drivers/${ICE_DRVDIR}/${ICE_DRVOPT}
${ICE_SANDBOX}/cicecore/cicedynB/dynamics
${ICE_SANDBOX}/cicecore/cicedynB/general
${ICE_SANDBOX}/cicecore/cicedynB/analysis
${ICE_SANDBOX}/cicecore/cicedynB/infrastructure
${ICE_SANDBOX}/cicecore/cicedynB/infrastructure/io/${IODIR}
${ICE_SANDBOX}/cicecore/cicedynB/infrastructure/comm/${ICE_COMMDIR}
${ICE_SANDBOX}/cicecore/shared
${ICE_SANDBOX}/cicecore/shared/constants/${ICE_CONSTOPT}
${ICE_SANDBOX}/icepack/columnphysics
${ICE_SANDBOX}/icepack/columnphysics/constants/${ICE_CONSTOPT}
EOF

#---------------------------------------------
#-- Build / make
#---------------------------------------------

# make log-file. ONLY stdout if empty
setenv ICE_BLDLOG_FILE
setenv ICE_BLDLOG_FILE cice.${compiler}.buildlog

# makdep: Dependicies
# Use Makefile inline compile: define SCC and ICE_CASEDIR (path to makdep.c)
setenv SCC gcc
setenv ICE_CASEDIR ${ICE_SANDBOX}/configuration/scripts
setenv OBJS_DEPGEN $ICE_CASEDIR/makdep.c

#  # ... or Remember to comment out OBJS_DEPGEN in Makefile
#  echo "building makdep using GNU compiler"
#  gcc -o makdep ${ICE_SANDBOX}/configuration/scripts/makdep.c || exit 2

echo "gmake clean"
gmake VPFILE=Filepath CPPDEFS="${ICE_CPPDEFS}" \
      -f  ${ICE_SANDBOX}/configuration/scripts/Makefile \
          OBJS_DEPGEN=${OBJS_DEPGEN} MACFILE=${MACFILE} clean | tee ${ICE_BLDLOG_FILE}

echo "building cice"
gmake -j ${ICE_MACHINE_BLDTHRDS} VPFILE=Filepath EXEC=${ICE_OBJDIR}/cice_${compiler}_${ICE_DRVOPT} \
      -f ${ICE_SANDBOX}/configuration/scripts/Makefile \
         OBJS_DEPGEN=${OBJS_DEPGEN} MACFILE=${MACFILE} | tee ${ICE_BLDLOG_FILE} 

if ($status != 0) then
  echo "COMPILE FAILED: ${0}"
  exit 99
else
  echo "*** COMPILE SUCCESSFUL ***"
endif

if ($CLEANUP == 1) then
  \rm ${ICE_OBJDIR}/{*.o,*.mod,*.d,makdep,Filepath}
endif

