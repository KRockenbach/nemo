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


source("plot_utils.R")

organism <- "Bnapus"

families <- c("MYBrelated", "Homeobox", "C2C2dof", "WRKY",
              "TCP", "bZIP", "MYB", "MADS", "bHLH", "Trihelix",
              "NAC", "G2like", "HSF", "AP2EREBP", "C2C2gata")

plot_lines <-function(pdf, legend=F, ...){
  x <- pdf[,"position"]
  min <- min(as.matrix(pdf[,families]))
  max <- max(as.matrix(pdf[,families]))
  plot(x, c(min,rep(max,(length(x)-1))), type="n", ...)
  cols <- rainbow(15)
  for (i in 1:15){
    fam <- families[i]
    lines(x=x, y=pdf[,fam], col=cols[i], lwd=2) 
  }
  if(legend){
    legend(x=(550), y=max, legend=families, fill = cols, ncol=1, cex=0.7, bty="n") 
  }
}

pdf <- read.table(paste0("../results/nemo/",organism,
                         "/masked_graphpart/attribs/promoter_TF_importance.tsv"),
                  header=T)
tdf <- read.table(paste0("../results/nemo/",organism,
                         "/masked_graphpart/attribs/terminator_TF_importance.tsv"),
                  header=T)
tdf$position <- tdf$position + 1 # shift position, because TTS is part of upstream sequence
# trim to +- 500 bp
pdf <- pdf[501:1501,]
tdf <- tdf[501:1501,]
pdf <- pdf[,c("position","mean", "sd", families)]
tdf <- tdf[,c("position","mean", "sd", families)]

######################
png("fig_S6.png", width=17, height=15, units="cm", res=1200)
par(mfrow=c(3,1), cex.axis=0.9)
padded <- c(rep(0,10000),pdf$mean,rep(0,10000))
freq <- c(1:length(padded))/length(padded)
plot(y=Mod(fft(padded)), x=freq, type="l", ylim=c(0,1), ylab="Amplitude", xlab="Frequency", axes=F)#, xlim=c(0,round(length(padded)/5,0)))
axis(1, at=(c(1:14)*1/15), labels=c("1/15", "2/15", "3/15", "4/15", "5/15", "6/15", "7/15",
                                    "8/15", "9/15", "10/15", "11/15", "12/15", "13/15", "14/15"))
axis(2)
legend("top", fill=transparent("seagreen",0.7), border=NA, legend="Filtered Frequency Ranges")
box()
transformed <- fft(padded)
period <- 15
for (i in c(1,2,13,14)){#1:(period-1)){
  f <- i/period
  at <- round(length(padded)*f,0)
  bw <- abs((i-(period/2))/(period/2))*100
  fl <- floor(at-bw)/length(padded)
  fh <- ceiling(at+bw)/length(padded)
  transformed[floor(at-bw):ceiling(at+bw)] <- 0
  polygon(x=c(fl,fl,fh,fh),
          y=c(-1,2,2,-1), border=NA, col=transparent("seagreen",0.7))
}
imp <- (pdf$mean)
plot(y=imp, x=c(-500:500), col="grey", type="l", ylab="Importance", xlab="Position Relative to TSS")
lines(y=(Re(fft(transformed, inverse=T))[10001:11001])/length(padded),
      x=c(-500:500), col="blue")
legend("topleft", lty=1, col=c("grey", "blue"), legend=c("Raw Importance", "Frequency Filter Applied"),
       bty="n")
N=45  #75
P=1
imp_filtered <- sgolayfilt(imp, p = P, n = N)
plot(imp, x=c(-500:500), col="grey", type="l",  ylab="Importance", xlab="Position Relative to TSS")
lines(y=imp_filtered, x=c(-500:500), col="blue")
legend("topleft", lty=1, col=c("grey", "blue"), legend=c("Raw Importance", "Savitzky-Golay Filter Applied"),
       bty="n")
dev.off()
########################

for (organism in c("Bnapus", "Athaliana")){
  pdf <- read.table(paste0("../results/nemo/",organism,
                           "/masked_graphpart/attribs/promoter_TF_importance.tsv"),
                    header=T)
  tdf <- read.table(paste0("../results/nemo/",organism,
                           "/masked_graphpart/attribs/terminator_TF_importance.tsv"),
                    header=T)
  tdf$position <- tdf$position + 1 # shift position, because TTS is part of upstream sequence
  # trim to +- 500 bp
  pdf <- pdf[501:1501,]
  tdf <- tdf[501:1501,]
  pdf <- pdf[,c("position","mean", "sd", families)]
  tdf <- tdf[,c("position","mean", "sd", families)]
  ##############################
  if (organism == "Bnapus"){
    figname <- "fig_S7.png"
  } else {
    figname <- "fig_S8.png"
  }
  png(figname, width=17, height=20, units="cm", res=1200)
  par(mfrow=c(5,1), mar=c(2,6,0.1,6), xpd=T)
  plot_lines(pdf, ylab="Raw Importance", xlab="")
  
  pdf.proc <- pdf
  tdf.proc <- tdf
  ## savitzky-golay filtering
  for (c in 2:ncol(pdf)){
    R=1
    N=45 #75
    P=1
    pdf.proc[,c] <- sgolayfilt(pdf.proc[,c], p = P, n = N)
    tdf.proc[,c] <- sgolayfilt(tdf.proc[,c], p = P, n = N)
  }
  plot_lines(pdf.proc, ylab="Savitzky-Golay Filtered", xlab="")
  ### mean normalize
  for (f in families){
    concat <- c(pdf.proc[,f], tdf.proc[,f])
    concat <- mean_norm(concat) #!
    pdf.proc[,f] <- concat[1:length(pdf.proc[,f])]
    tdf.proc[,f] <- concat[(length(pdf.proc[,f])+1):length(concat)]
  }
  plot_lines(pdf.proc, ylab="Mean Noramlized", xlab="", legend=T)
  #### mean normalize (across sequence) mean (across samples)
  concat_mean <- c(pdf.proc$mean, tdf.proc$mean)
  concat_mean <- mean_norm(concat_mean) #!
  pdf.proc$mean <- concat_mean[1:length(pdf.proc$mean)]
  tdf.proc$mean <- concat_mean[(length(pdf.proc$mean)+1):length(concat_mean)]
  
  ## standardize
  for (f in families){
    pdf.proc[,f] <- (pdf.proc[,f] - pdf.proc$mean)
    tdf.proc[,f] <- (tdf.proc[,f] - tdf.proc$mean)
  }
  plot_lines(pdf.proc, ylab="Overall Mean Subtracted\nfrom Group Mean", xlab="")
  for (f in families){
    pdf.proc[,f] <- (pdf.proc[,f]/pdf.proc$sd)
    tdf.proc[,f] <- (tdf.proc[,f]/tdf.proc$sd)
  }
  par(mar=c(5,6,0.1,6))
  plot_lines(pdf.proc, ylab="Divided By Overall\nStandard Deviation", xlab="Position Relative to TSS")
  dev.off()
  
}

