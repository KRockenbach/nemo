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


source("./plot_utils.R")

png("fig_2.png", width=17, height=17, units="cm", res=1200)

lab_cex=1
ax_cex=0.9
text_cex=1.1
fig_lab_cex=1.5

layout(
  matrix(c(1,1,2,2,
           4,0,7,0,
           3,5,6,8), ncol=4, byrow=TRUE), 
  widths=c(3,0.3,3,0.3), 
  heights=c(3,0.3,3)
)


# 2A
df = read.table("../results/rsq_df.tsv", header=T, sep="\t")
sub=df[df$partitioning=="graphpart_Bn",]
sub$model[is.na(sub$valid_fold)] <- "nemo90"
par(mar=c(5,4,3,0), mgp=c(1.4,0.4,0), tck=-0.03, xpd=T)
plot(x=c(0.5:5.5), y=c(min(sub$rsq), rep(max(sub$rsq),5)), 
     ylab=expression("Performance on"~italic("B. napus")~"Test Set (r"^2*")"),
     xlab="", cex.lab=lab_cex, cex.axis = ax_cex, type="n", axes=F)
model_names <- c("nemo90", "nemo", "nemo_rand_prom", "nemo_rand_term", "xpresso")
for (m in 1:5){
  model <- model_names[m]
  boxplot(sub$rsq[sub$model==model], at=m, add=T, axes=F)
}

axis(side=2, at=seq(0.05, 0.55, 0.05))
text(x=c(1:5), y=0.58, labels=rep("N = 10", 4), col="grey70")
for (m in 1:5){
  model <- model_names[m]
  med_rsq <- as.character(round(median(sub$rsq[sub$model==model]),2))
  text(x=m, y=0.3, labels=bquote("Q"[2]~"="~.(med_rsq)), cex=0.75, col="grey30")
}

text(x=1, y=-0.0375, labels=expression(italic("n")*"emo"["90"]), cex=0.8)
text(x=2, y=-0.0375, labels=expression(italic("n")*"emo"["80"]), cex=0.8)
Lines <- list(bquote(italic("n")*"emo"["80"]),bquote("with"),bquote("Randomized"),bquote("Promoter"))
text(x=3, y=c(0, -0.025,-0.05,-0.075), labels=do.call(expression, Lines), cex=0.8)
Lines <- list(bquote(italic("n")*"emo"["80"]),bquote("with"),bquote("Randomized"),bquote("Terminator"))
text(x=4, y=c(0, -0.025,-0.05,-0.075), labels=do.call(expression, Lines), cex=0.8)
text(x=5, y=-0.0375, labels="Xpresso", cex=0.8)

fig_label(expression(bold("A")), cex=fig_lab_cex)




# 2B

df = read.table("../results/rsq_df.tsv", header=T, sep="\t")

sub=df[df$masking=="masked" & df$partitioning=="graphpart" & df$test_organism %in% c("Athaliana", "BnapusDS"), ]
cols=c("seagreen2","yellow")


par(mar=c(3,3,3,0), mgp=c(1.1,0.6,0.2), xpd=T)

plot(x=c(1.2,1.8,3.2,3.8), y=rep(0.5, 4), ylim=c(0.3,0.7), xlim=c(0.8,4.2), axes=F, xlab="", ylab="", type="n")
tr_At_te_At <- sub$rsq[sub$train_organism=="Athaliana" & sub$test_organism=="Athaliana"]
tr_At_te_Bn <- sub$rsq[sub$train_organism=="Athaliana" & sub$test_organism=="BnapusDS"]
tr_Bn_te_At <- sub$rsq[sub$train_organism=="Bnapus" & sub$test_organism=="Athaliana"]
tr_Bn_te_Bn <- sub$rsq[sub$train_organism=="Bnapus" & sub$test_organism=="BnapusDS"]
boxplot(tr_At_te_At, col="seagreen2", at=1.2, boxwex=0.5, border="darkgreen", add=T, frame=F, axes=F)
boxplot(tr_At_te_Bn, col="yellow", at=1.8, boxwex=0.5, border="gold4", add=T, frame=F, axes=F)
boxplot(tr_Bn_te_At, col="seagreen2", at=3.2, boxwex=0.5, border="darkgreen", add=T, frame=F, axes=F)
boxplot(tr_Bn_te_Bn, col="yellow", at=3.8, boxwex=0.5, border="gold4", add=T, frame=F, axes=F)
Lines <- list(bquote("Q"[2]~"="~.(round(median(tr_At_te_At),2))),
              bquote("Q"[2]~"="~.(round(median(tr_At_te_Bn),2))),
              bquote("Q"[2]~"="~.(round(median(tr_Bn_te_At),2))),
              bquote("Q"[2]~"="~.(round(median(tr_Bn_te_Bn),2))))
text(x=c(1.2,1.8,3.2,3.8), y=0.33, labels=do.call(expression, Lines), cex=0.75, col="grey30")

axis(1, at=c(1.5, 3.5),
     labels=c(expression(italic("A. thaliana")),
              expression(italic("B. napus"))),
     cex.axis=1, tick=F, cex.lab=1.3)
axis(2, tck=-0.03, cex.axis=ax_cex, line=-0.5)


title(ylab=expression('Performance on Test Set (r ' ^ 2*")"), 
      xlab="Training Organism", cex.lab=lab_cex, line=1.5)

legend("top", legend=c(expression(italic("A. thaliana")), expression(italic("B. napus"))), 
       fill=c("seagreen2", "yellow"), cex=0.8, title="Test Organism", border=c("darkgreen", "gold4"),
       bty="n")

orgs <- c("Athaliana", "Bnapus")


print(shapiro.test(sub$rsq[sub$train_organism=="Athaliana" & sub$test_organism=="Athaliana"]))
print(shapiro.test(sub$rsq[sub$train_organism=="Athaliana" & sub$test_organism=="BnapusDS"]))
print(shapiro.test(sub$rsq[sub$train_organism=="Bnapus" & sub$test_organism=="Athaliana"]))
print(shapiro.test(sub$rsq[sub$train_organism=="Bnapus" & sub$test_organism=="BnapusDS"]))
## shapito test is significant for models trained on Athaliana and tested on Bnapus!
## --> wilcoxon test!


p_vals <- c()
for (tn in 1:2){
  tno <- orgs[tn]
  p_vals[tn] <- wilcox.test(sub$rsq[sub$train_organism==tno & sub$test_organism=="Athaliana"],
                            sub$rsq[sub$train_organism==tno & sub$test_organism=="BnapusDS"])$p.val
  print(tno)
  print(wilcox.test(sub$rsq[sub$train_organism==tno & sub$test_organism=="Athaliana"],
                    sub$rsq[sub$train_organism==tno & sub$test_organism=="BnapusDS"]))
}


labs <- c()
for (p in p_vals){
  if (p <= 0.001){
    labs <- c(labs,"***")
  } else if (p <= 0.01) {
    labs <- c(labs,"**")
  } else if (p <= 0.05) {
    labs <- c(labs,"*")
  } else {
    labs <- c(labs,"ns")
  }
}
text(x=c(1.5, 3.5), y=c(0.515,0.62), labels=labs, 
     col="grey40", cex=0.85)
text(x=c(1.2,1.8,3.2,3.8), y=rep(0.7,4), labels=rep("N = 10", 4), 
     col="grey70")#, cex=0.85)

fig_label(expression(bold("B")), cex=fig_lab_cex)




#  2C

# Bnapus concatenated
path <- "../results/nemo/Bnapus/masked_graphpart/nemoB_preds/concat_preds.tsv"
df <- read.table(path, sep="\t", header=T)
# Gaussian KDE
df$Density <- get_density(df$Actual, df$Predicted, n = 100, h = c(1, 1))
# scale density to range (0,1)
df$Density <- df$Density/max(df$Density)
# scale to integer range [1,1000]
df$Density <- 1 + floor(df$Density*999)
# sort df by density, so that densest points get drawn last
df <- df[order(df$Density),]

par(mar=c(3,3,0,0), mgp=c(1.4,0.4,0), tck=-0.02)
plot(df$Actual, df$Predicted, col=viridis(1000)[df$Density], pch=19,
     xlab = expression("Observed Expression [log"[10]*"(TPM + 0.1)]"),
     ylab = expression("Predicted Expression [log"[10]*"(TPM + 0.1)]"),
     cex.lab=lab_cex, axes=F, xlim=c(-1,4), ylim=c(min(df$Predicted),4))
for (i in 1:100){
  bottom = 2+((i-1)/100)
  top = 2+((i-1)/100)+(1/100)
  polygon(x=c(4,4,4.1,4.1), y=c(bottom,top,top,bottom), 
          col=viridis(100)[i], border=NA)
}
lines(x=c(3.95,4.1), y=c(2,2))
text(x=3.91, y=2, labels="Low", adj=1, cex=0.95, col="grey50")
lines(x=c(3.95,4.1), y=c(3,3))
text(x=3.91, y=3, labels="High", adj=1, cex=0.95, col="grey50")
text(x=4.1, y=3.6, labels="Kernel\nDensity\nEstimation", cex=0.95, adj=1)
box()
axis(side=1, cex.axis=0.85, tck=-0.02)
axis(side=2, cex.axis=0.85, tck=-0.02)

text(x=1.5, y=3.8, labels=expression(italic("B. napus")), cex=1.2)

r <- cor(x=df$Actual, y=df$Predicted, method="pearson")
rsq <- r^2
print(rsq)
n.lab <- paste("N = ", as.character(length(df$Actual)), sep="")
text(x=c(-0.5), y=c(3.8), labels=bquote("r"^2*" = "~.(format(round(rsq,3), nsmall=3))), 
     col="black", cex=text_cex)
text(x=3.6, y=-1, labels=n.lab, 
     col="grey70", cex=text_cex)


#  2C upper margin
dx <- density(df$Actual, cut=F)
par(mar=c(0,3,0,0))
plot(x=c(-1,dx$x), y=c(0,dx$y)-0.1, type="l", ylab="", xlab="", axes=F, xlim=c(-1,4))
polygon(x=c(-1,dx$x), y=c(0,dx$y)-0.1, col="grey90", border=NA)
lines(x=c(-1,dx$x), y=c(0,dx$y)-0.1)
fig_label(expression(bold("C")), cex=fig_lab_cex)


# 2C right margin
dy <- density(df$Predicted, cut=F)
par(mar=c(3,0,0,0))
plot(x=c(0,dy$y)-0.1, y=c(min(df$Predicted),dy$x), type="l", xlab="", ylab="", axes=F, ylim=c(min(dy$x),4))
polygon(x=c(0,dy$y)-0.1, y=c(min(dy$x),dy$x), col="grey90", border=NA)
lines(x=c(0,dy$y)-0.1, y=c(min(df$Predicted),dy$x))




# 2 D

# Athaliana concatenated
# TRAINED AND TESTED ON ATHALIANA

path <- "../results/nemo/Athaliana/masked_graphpart/nemoA_preds/concat_preds.tsv"
df <- read.table(path, sep="\t", header=T)
# Gaussian KDE
df$Density <- get_density(df$Actual, df$Predicted, n = 100, h = c(1, 1))
# scale density to range (0,1)
df$Density <- df$Density/max(df$Density)
# scale to integer range [1,1000]
df$Density <- 1 + floor(df$Density*999)
# sort df by density, so that densest points get drawn last
df <- df[order(df$Density),]

par(mar=c(3,3,0,0), mgp=c(1.4,0.4,0), tck=-0.02)
plot(df$Actual, df$Predicted, col=viridis(1000)[df$Density], pch=19,
     xlab = expression("Observed Expression [log"[10]*"(TPM + 0.1)]"),
     ylab = expression("Predicted Expression [log"[10]*"(TPM + 0.1)]"),
     cex.lab=lab_cex, axes=F, xlim=c(-1,4), ylim=c(min(df$Predicted),4))
for (i in 1:100){
  bottom = 2+((i-1)/100)
  top = 2+((i-1)/100)+(1/100)
  polygon(x=c(4,4,4.1,4.1), y=c(bottom,top,top,bottom), 
          col=viridis(100)[i], border=NA)
}
lines(x=c(3.95,4.1), y=c(2,2))
text(x=3.91, y=2, labels="Low", adj=1, cex=0.95, col="grey50")
lines(x=c(3.95,4.1), y=c(3,3))
text(x=3.91, y=3, labels="High", adj=1, cex=0.95, col="grey50")
text(x=4.1, y=3.6, labels="Kernel\nDensity\nEstimation", cex=0.95, adj=1)
box()
axis(side=1, cex.axis=0.85, tck=-0.02)
axis(side=2, cex.axis=0.85, tck=-0.02)

text(x=1.5, y=3.8, labels=expression(italic("A. thaliana")), cex=1.2)

r <- cor(x=df$Actual, y=df$Predicted, method="pearson")
rsq <- r^2

n.lab <- paste("N = ", as.character(length(df$Actual)), sep="")
text(x=c(-0.5), y=c(3.8), labels=bquote("r"^2*" = "~.(format(round(rsq,3), nsmall=3))), 
     col="black", cex=text_cex)
text(x=3.6, y=-1, labels=n.lab, 
     col="grey70", cex=text_cex)


#  2D upper margin
dx <- density(df$Actual, cut=F)
par(mar=c(0,3,0,0))
plot(x=c(-1,dx$x), y=c(0,dx$y)-0.1, type="l", ylab="", xlab="", axes=F, xlim=c(-1,4))
polygon(x=c(-1,dx$x), y=c(0,dx$y)-0.1, col="grey90", border=NA)
lines(x=c(-1,dx$x), y=c(0,dx$y)-0.1)
fig_label(expression(bold("D")), cex=fig_lab_cex)


# 2D right margin
dy <- density(df$Predicted, cut=F)
par(mar=c(3,0,0,0))
plot(x=c(0,dy$y)-0.1, y=c(min(df$Predicted),dy$x), type="l", xlab="", ylab="", axes=F, ylim=c(min(dy$x),4))
polygon(x=c(0,dy$y)-0.1, y=c(min(dy$x),dy$x), col="grey90", border=NA)
lines(x=c(0,dy$y)-0.1, y=c(min(df$Predicted),dy$x))

dev.off()

  
