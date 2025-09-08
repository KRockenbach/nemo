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
#title: get_gene_stats.awk
#description: extracts gene statistics from GFF file for calculation of Xpresso halflife features (gene ID, total gene length, number of CDS features, total exon length)
#author: Kevin Rockenbach
#email: kevin.rockenbach@ag.uni-giessen.de
#date: 2025-08-28
#version: 1.0.0
#usage: awk -f get_gene_stats.awk <GFF> > output.tsv
#=========================================================================================================


# based on https://github.com/nathanweeks/scripts/blob/master/intron-length.awk


BEGIN {
    FS = OFS = "\t"
    processed = 0
}

/^#/ || NF == 0 { next } # skip comments or blank lines

$3 == "mRNA" || $3 == "transcript" {

    if (processed > 0){ # print stats of previous mRNA before processing next
        print id, total_gene_len, num_cds, total_exon_len
    }

    total_exon_len = num_cds = 0
    total_gene_len = $5 - $4 + 1
    split($9, id_split1, ";"a)
    split(id_split1[1], id_split2, "=")
    split(id_split2[2], id_split3, ".")
    id  = id_split3[1]
    processed += 1
}

$3 == "exon" {
    # for length calculation strandedness doesn't matter
    exon_len = $5 - $4 + 1 # end - start + 1
    total_exon_len += exon_len
}

$3 == "CDS" {
    num_cds += 1
}


END { # print stats of last processed mRNA
    print id, total_gene_len, num_cds, total_exon_len
}
