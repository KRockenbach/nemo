#!/usr/bin/awk -f

##################################################################################
#
# MIT License
#
# Copyright (c) 2025 Kevin Rockenbach, Agnieszka Golicz
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#
##################################################################################


#=========================================================================================================
#title: separate_motifs.awk
#description: separates out individual pfms from combined JASPAR pfm file
#author: Kevin Rockenbach
#email: kevin.rockenbach@ag.uni-giessen.de
#date: 2025-08-28
#version: 1.0.0
#usage: awk -f separate_motifs.awk <family level combined pfm file>
#=========================================================================================================


# This script takes in a collection of PWMs from JASPAR and separates them out
# output directory should be piped into awk, concatenated PWM file is second input

BEGIN{
    OFS="\t"
    outdir="."
    split("", TF_arr) # initialize empty array
    skip="False"
}

{
    if (NR==1 && NF==1) { # out dir path
        outdir=$0
    } else {
        if ($0 ~ "^>") {
            gsub(/ /, "\t", $0)
            num_a=split($2,a,".")
            if (num_a==4){
                TF=a[3]"."a[4]
            } else {
                TF=a[3]
            }
            motif=a[1]
            version=a[2]
            for (key in TF_arr){
                if (key == TF && version < TF_arr[TF]){
                    skip="True"
                    next
                }
            }
            TF_vers[TF]=version
            outfile=outdir"/"TF".pwm"
            printf("echo '# %s %s.%s' > %s\n", TF, motif, version, outfile) | "bash"
        }

        if (skip == "True") {
            n = NR + 4
            skip = "False"
        }

        if (NR < n) {next}


        if ($0 !~ "^>") {
            gsub(/^\s*/, "", $0)
            gsub(/\s\s*/, "\t", $0)
            printf("echo '%s' >> %s\n", $0, outfile) | "bash"
        }
    }
}

END { close("bash") }
