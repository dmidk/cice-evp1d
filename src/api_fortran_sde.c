/*!===============================================================================
  ! Copyright (C) 2023, Intel Corporation
  ! Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
  ! 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
  ! 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
  ! 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
  ! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS “AS IS” AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
  !===============================================================================*/
#include<stdio.h>

#define SSC_MARK( MARK_ID ) \
        __asm__ __volatile__ ( \
        "push %rbx \n\t" \
        "movl $"#MARK_ID", %ebx \n\t" \
        ".byte 0x64 \n\t" \
        ".byte 0x67 \n\t" \
        ".byte 0x90 \n\t" \
        "pop %rbx" \
        );

void fortran_ssc_start() {
       printf ("Entered SSC START\n");
       SSC_MARK(0x111);
}

void fortran_ssc_stop() {
       printf ("Entered SSC STOP\n");
       SSC_MARK(0x222);
}
