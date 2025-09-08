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

png("fig_S4.png", width=17, height=8, units="cm", res=1200)

par(mfrow=c(1,2), mar=c(0,0,0,0), oma=c(0,0,0,0))

fig_lab_cex=1.5

# A

df <- read.table("../results/nemo/Bnapus/masked_graphpart_Bn/logs/lr_log_412_53543_t_0_v_1.csv",
                 header=T, sep="\t")

steps_per_epoch = 412

epoch <- 1
total_step <- c(1)
for (i in 2:nrow(df)){
  if (df$step[i] < df$step[i-1]){
    epoch <- epoch + 1
  }
  total_step[i] <- (steps_per_epoch * (epoch-1)) + (df$step[i] + 1)
}

df$total_step <- total_step

par(mar=c(2.5,2.9,2,2.9), xpd=T, tck=-0.03, mgp=c(1.05,0.4,0))

plot(lr~total_step, data=df, type="l", axes=F, 
     ylab="", xlab="", col="purple3", ylim=c(0.00,0.08), lwd=1.5, cex.lab=0.9)
text(x=0, y=0.08, labels=c("Learning Rate"), adj=0.7, cex=0.9)
title(xlab="Epoch", cex.lab=0.9)#, line=1)
axis(side=1, at=seq(0,(412*6),412), labels=c(0:6), cex.axis=0.7, lwd=1.5)
axis(side=2, col="purple3", at=seq(0.00,0.07,0.01), las=2, cex.axis=0.7, lwd=1.5)


par(new=T)
plot(beta_1~total_step, data=df, type="l", axes=F, 
     ylab="", xlab="", col="green3", ylim=c(0.89,0.938), lwd=1.5)
text(x=412*6.3, y=0.937, labels=expression("Momentum ( "*beta[1]*" )"), cex=0.9)
axis(side=4, col="green3", at=seq(0.89, 0.93, 0.01), las=2, cex.axis=0.8, lwd=1.5)


fig_label(expression(bold("A")), cex=fig_lab_cex)


# B

df <- read.table("../results/nemo/Bnapus/masked_graphpart_Bn/logs/trainlog_t_0_v_1.csv",
                 header=T, sep="\t")

df$epoch <- df$epoch + 1

epoch <- c(1:6)
train_median <- c()
train_mad <- c()
val_median <- c()
val_mad <- c()
for (e in epoch){
  sub <- df[df$epoch==e,]
  train_median <- c(train_median, median(sub$loss))
  train_mad <- c(train_mad, mad(sub$loss))
  val_median <- c(val_median, median(sub$val_loss))
  val_mad <- c(val_mad, mad(sub$val_loss))
}

par(mar=c(2.5,3,0.5,0), xpd=T, tck=-0.02, mgp=c(1.45,0.35,0), cex.lab=0.8, cex.axis=0.7)
new_df <- data.frame(epoch, train_median, train_mad, val_median, val_mad)
plot(val_median~epoch, data=new_df, type="n", 
     ylim=c(0.19, 0.45), xlab="", ylab="Huber Loss", axes=F)
title(xlab="Epoch", cex.lab=0.9, line=1.05)
axis(side=1, lwd=1.5)
axis(side=2, las=2, lwd=1.5)
polygon(x=c(new_df$epoch, rev(new_df$epoch)), 
        y=c(new_df$val_median+new_df$val_mad, rev(new_df$val_median-new_df$val_mad)),
        col=transparent("lightskyblue", 0.6), border=NA)
polygon(x=c(new_df$epoch, rev(new_df$epoch)), 
        y=c(new_df$train_median+new_df$train_mad, rev(new_df$train_median-new_df$train_mad)),
        col=transparent("goldenrod1", 0.6), border=NA)
lines(x=new_df$epoch, y=new_df$val_median, col="blue", lwd=1.5)
lines(x=new_df$epoch, y=new_df$train_median, col="darkorange", lwd=1.5)

polygon(x=c(2.4,2.4,2.6,2.6), y=c(0.42,0.43,0.43,0.42), 
        col=transparent("goldenrod1", 0.6), border=NA)
lines(x=c(2.4,2.6), y=c(0.425,0.425), col="darkorange", lwd=1.5)
polygon(x=c(2.4,2.4,2.6,2.6), y=c(0.405,0.415,0.415, 0.405), 
        col=transparent("lightskyblue", 0.6), border=NA)
lines(x=c(2.4,2.6), y=c(0.41,0.41), col="blue", lwd=1.5)
text(x=2.7, y=c(0.425,0.41), labels=c("Training Loss", "Validation Loss"), 
     adj=0, cex=0.8)
text(x=3.2, y=0.39, labels="Median \u00B1 Median Absolute Deviation", 
     adj=0.5, col="grey40", cex=0.7)

fig_label(expression(bold("B")), cex=fig_lab_cex)
dev.off()



