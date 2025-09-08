#!/usr/bin/Rscript --vanilla

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


rsq <- read.table("../results/rsq_df.tsv", header=T, sep="\t")
sub <- rsq[(rsq$partitioning=="graphpart" | rsq$partitioning=="random"),]
sub <- sub[sub$test_organism=="Bnapus",]
sub <- sub[sub$train_organism=="Bnapus",]


png("fig_S43.png", height=10, width=17, res=1200, units="cm")

par(mar=c(5,5,0,0))
plot(x=c(0:6), y=c(rep(0.2,4), rep(0.8,3)), axes=F, type="n", 
     ylab=expression("Model Performance (r "^2*")"),
     xlab="Input Masking", cex.lab=1.2)
axis(1, labels = c("CDS masked", "Clear sequence", "Non-coding masked"), at=c(1,3,5))
axis(2)
mask_types <- c("masked", "clear", "only_cds")
part_types <- c("graphpart", "random")
for(m in 1:3){
  masking <- mask_types[m]
  at <- 1 + ((m-1)*2)
  at_p <- at - 0.3
  at_r <- at + 0.3
  gp <- sub$rsq[sub$masking==masking & sub$partitioning=="graphpart"]
  rand <- sub$rsq[sub$masking==masking & sub$partitioning=="random"]
  boxplot(gp, at=at_p, add=T, col="purple", axes=F)
  boxplot(rand, at=at_r, add=T, col="cyan", axes=F)
  text(x=at, y=0.3, labels=bquote(Delta*"Q"[2]~"="~.(round(median(rand)-median(gp),3))))
}
legend("topright", title="Partitioning", fill=c("purple", "cyan"), legend=c("Gene family-guided", "Random"))

dev.off()