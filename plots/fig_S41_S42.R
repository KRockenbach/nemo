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


png("fig_S41.png", width=17, height=10, units="cm", res=1200)
lab_cex=0.8
ax_cex=0.7
text_cex=0.8
fig_lab_cex=1

layout(matrix(c(1,2,5,
                3,4,5), byrow=T, ncol=3),
       widths=c(3,3,2),
       heights=c(3,3))

par(mar=c(2.5,3,0.1,0.1))

for (org in c("Bnapus", "Athaliana")){
  init <- strsplit(org, split="", fixed=T)[[1]][1]
  df <- read.table(paste0("../results/nemo/",org,"/masked_graphpart/nemo",init,"_preds/concat_preds.tsv"), header=T, sep="\t")
  err <- (df$Predicted - df$Actual)
  exp <- read.table(paste0("../data/", org, "/derived_data/expr_matrix.tsv"), header=T, sep="\t", row.names=1)
  exp[is.na(exp)] <- 0.0 
  exp <- exp[df$ID,]
  genes <- rownames(exp)
  exp <- exp[,!(colnames(exp) == "Gene")]
  med_exp <- df$Actual
  dense <- get_density(med_exp, err, n = c(300,100))
  new_df <- data.frame(med_exp, err, dense)
  new_df$dense <- new_df$dense/max(new_df$dense)
  # scale to integer range [1,1000]
  new_df$dense <- 1 + floor(new_df$dense*999)
  new_df <- new_df[order(new_df$dense),]
  cols <- wes_palette("Zissou1", 1000, type = "continuous")
  plot(x=new_df$med_exp, y=new_df$err, type="n", col=cols[new_df$dense],
       ylab = "", xlab = "", cex.lab=lab_cex, cex.axis=ax_cex, axes=F, xlim=c(-1,4))
  
  if (org == "Athaliana"){
    ylim=c(-3,2)
  } else {
    ylim=c(-4,4)
  }
  axis(1, tck=-0.03, at=seq(-1, 4, 1), labels=rep("", length(seq(-1, 4, 1))), cex.axis=ax_cex, line=0)
  axis(2, tck=-0.03, at=seq(ylim[1], ylim[2], 1), labels=rep("", length(seq(ylim[1], ylim[2], 1))), cex.axis=ax_cex, line=0)
  axis(1, lwd=0, at=seq(-1, 4, 1), cex.axis=ax_cex, line=-0.7)
  axis(2, lwd=0, at=seq(ylim[1], ylim[2], 1), cex.axis=ax_cex, line=-0.5)
  title(ylab = expression("Predicted − Observed ("~epsilon~")"),
        xlab = expression("Observed Median Expression ["~log10("TPM + 0.1")~"]"), cex.lab=lab_cex, line=1.2)

  cutoff <- median(abs(err))
  polygon(x=c(-2,-2,5,5), y=c(-5,5,5,-5), col="darkolivegreen1", border=NA) #lightseagreen
  polygon(x=c(-2,-2,5,5), y=c(-5,-cutoff,-cutoff,-5), col="darkslategray1", border=NA)
  polygon(x=c(-2,-2,5,5), y=c(5,cutoff,cutoff,5), col="lemonchiffon", border=NA) #red4

  points(x=new_df$med_exp, y=new_df$err, pch=19, col=cols[new_df$dense])
  abline(lm(err~med_exp, data=new_df), col="black", lwd=2)
  print(cor.test(new_df$med_exp, new_df$err))
  r <- round(cor(new_df$med_exp, new_df$err, method="pearson"),3)
  if (org == "Athaliana"){
    legend("topright", lwd=c(2,NA), col="black", 
           legend=c("linear regression", expression(italic("A. thaliana"))), bty="n",
           text.col=c("black", "grey60"))
    legend("bottomleft", pch=NA, legend=bquote(rho["P"]~"="~.(as.character(round(r,3)))~"***"), bty="n", text.col="black") 
  } else {
    legend("topright", lwd=c(2,NA), col="black", 
           legend=c("linear regression", expression(italic("B. napus"))), bty="n",
           text.col=c("black", "grey60"))
    legend("bottomleft", pch=NA, legend=bquote(rho["P"]~"="~.(as.character(round(r,3)))~"***"), bty="n", text.col="black") 
  }
  if (org == "Athaliana"){
    fig_label(expression(bold("B")), cex=fig_lab_cex) 
  } else {
    fig_label(expression(bold("A")), cex=fig_lab_cex) 
  }  
  
}



#######################
#######################
par(mar=c(2.5,3,0.1,0.1))

for (org in c("Bnapus", "Athaliana")){
  init <- strsplit(org, split="", fixed=T)[[1]][1]
  df <- read.table(paste0("../results/nemo/",org,"/masked_graphpart/nemo",init,"_preds/concat_preds.tsv"), header=T, sep="\t")
  err <- (df$Predicted - df$Actual)
  exp <- read.table(paste0("../data/", org, "/derived_data/expr_matrix.tsv"), header=T, sep="\t", row.names=1)
  exp[is.na(exp)] <- 0.0 
  exp <- exp[df$ID,]
  genes <- rownames(exp)
  exp <- exp[,!(colnames(exp) == "Gene")]
  tau <- apply(exp, MARGIN=1, FUN=get_tau)
  dense <- get_density(tau, err, n = c(300,100))
  new_df <- data.frame(tau, err, dense)
  new_df$dense <- new_df$dense/max(new_df$dense)
  # scale to integer range [1,1000]
  new_df$dense <- 1 + floor(new_df$dense*999)
  new_df <- new_df[order(new_df$dense),]
  cols <- wes_palette("Zissou1", 1000, type = "continuous")
  plot(x=new_df$tau, y=new_df$err, type="n", col=cols[new_df$dense],
       ylab = "", xlab = "", cex.lab=lab_cex, cex.axis=ax_cex, axes=F)
  
  if (org == "Athaliana"){
    ylim=c(-3,2)
  } else {
    ylim=c(-4,4)
  }
  axis(1, tck=-0.03, at=seq(0.3, 1, 0.1), labels=rep("", length(seq(0.3, 1, 0.1))), cex.axis=ax_cex, line=0)
  axis(2, tck=-0.03, at=seq(ylim[1], ylim[2], 1), labels=rep("", length(seq(ylim[1], ylim[2], 1))), cex.axis=ax_cex, line=0)
  axis(1, lwd=0, at=seq(0.3, 1, 0.1), cex.axis=ax_cex, line=-0.7)
  axis(2, lwd=0, at=seq(ylim[1], ylim[2], 1), cex.axis=ax_cex, line=-0.5)
  title(ylab = expression("Predicted − Observed ("~epsilon~")"),
        xlab = expression("Tissue-Specificity ("~tau~")"), cex.lab=lab_cex, line=1.2)
  ###################
  cutoff <- median(abs(err))
  polygon(x=c(0,0,2,2), y=c(-5,5,5,-5), col="darkolivegreen1", border=NA) #lightseagreen
  polygon(x=c(0,0,2,2), y=c(-5,-cutoff,-cutoff,-5), col="darkslategray1", border=NA)
  polygon(x=c(0,0,2,2), y=c(5,cutoff,cutoff,5), col="lemonchiffon", border=NA) #red4
  #####################
  points(x=new_df$tau, y=new_df$err, pch=19, col=cols[new_df$dense])
  abline(lm(err~tau, data=new_df), col="black", lwd=2)
  print(cor.test(new_df$tau, new_df$err))
  r <- round(cor(new_df$tau, new_df$err, method="pearson"),3)
  
  if (org == "Athaliana"){
    legend("topleft", lwd=c(2,NA), col="black", 
           legend=c("linear regression", expression(italic("A. thaliana"))), bty="n",
           text.col=c("black", "grey60"))
    legend("bottomleft", pch=NA, legend=bquote(rho["P"]~"="~.(as.character(round(r,3)))~"***"), bty="n", text.col="black") 
  } else {
    legend("topleft", lwd=c(2,NA), col="black", 
           legend=c("linear regression", expression(italic("B. napus"))), bty="n",
           text.col=c("black", "grey60"))
    legend("bottomleft", pch=NA, legend=bquote(rho["P"]~"="~.(as.character(round(r,3)))~"***"), bty="n", text.col="black") 
  }
  if (org == "Athaliana"){
    fig_label(expression(bold("D")), cex=fig_lab_cex) 
  } else {
    fig_label(expression(bold("C")), cex=fig_lab_cex) 
  }  
  
}

par(mar=c(2.5,0.1,0.1,0.1))
plot(x=c(1:100), y=c(1:100), type="n", ylab="", xlab="", axes=F)
for (i in 1:1000){
  h = (50/1000)
  bottom=25+((i-1)*h)
  top=25+(i*h)
  left=23
  right=37
  polygon(x=c(left, left, right, right), y=c(bottom,top,top,bottom), col=cols[i], border=NA)
}
lines(x=c((left-0.1),(right+0.1)), y=c(25,25), col="black")
lines(x=c((left-0.1),(right+0.1)), y=c(75,75), col="black")
text(x=30, y=78, labels="High", col="grey30", cex=0.7)
text(x=30, y=22, labels="Low", col="grey30", cex=0.7)
text(x=30, y=85, labels="Kernel\nDensity\nEstimation", col="black", cex=0.8)

left=62
right=87
polygon(x=c(left,left,right,right), y=c(70,90,90,70), col="lemonchiffon", border=NA)
polygon(x=c(left,left,right,right), y=c(40,60,60,40), col="darkolivegreen1", border=NA)
polygon(x=c(left,left,right,right), y=c(10,30,30,10), col="darkslategray1", border=NA)
text(x=rep(((left+right)/2),3), y=c(80,50,20), labels=c("Over-\npredicted", "Well-\npredicted", "Under-\npredicted"), col="grey30", cex=0.7)


dev.off()













###########################################
###########################################
###########################################



png("fig_S42.png", width=17, height=9, units="cm", res=1200)

lab_cex=0.8
ax_cex=0.7
text_cex=0.8
fig_lab_cex=1.5

layout(
  matrix(c(1,2,3), ncol=3, byrow=TRUE), 
  widths=c(2.5,0.65,0.9), 
  heights=c(3)
)

org <- "Athaliana"
###############################
#A

init <- strsplit(org, split="", fixed=T)[[1]][1]
df <- read.table(paste0("../results/nemo/",org,"/masked_graphpart/nemo",init,"_preds/concat_preds.tsv"), header=T, sep="\t")

x <- df$Actual
y <- df$Predicted

err <- (y - x)
abs_err <- abs(err)
cutoff <- median(abs_err)



par(mar=c(3,3,2,0), mgp=c(1.4,0.4,0), tck=-0.03, xpd=F)
plot(x, y, type="n", xlab="", ylab="", axes=F)
axis(side=1, cex.axis=ax_cex, tck=-0.02)
axis(side=2, cex.axis=ax_cex, tck=-0.02)
title(xlab=expression("Observed Expression"~"[log"[10]*"(TPM + 0.1)]"),
      ylab=expression("Predicted Expression"~"[log"[10]*"(TPM + 0.1)]"), cex.lab=lab_cex)


polygon(x=c(-2,-2,5,5), y=c(-2,5,5,-2), col="darkolivegreen1", border=NA) #lightseagreen


x_under <- c(-2, 5, 5)
y_under <- c(((-2)-cutoff), (5-cutoff), ((-2)-cutoff))
polygon(x=x_under, y=y_under, col="darkslategray1", border=NA)

x_over <- c(-2, -2, 5)
y_over <- c(((-2)+cutoff), (5+cutoff), (5+cutoff))
polygon(x=x_over, y=y_over, col="lemonchiffon", border=NA) #red4


dense <- get_density(x, y, n = 100, h = c(1, 1))
# scale to [0,1]
dense <- dense/max(dense)
# get colors
cols <- grey_scale(dense)
x <- x[order(dense)]
y <- y[order(dense)]
cols <- cols[order(dense)]
points(x, y, col=cols, pch=20)

# create grey_scale legend
#####################
x_max <- max(x)
x_min <- min(x)
y_max <- max(y)
y_min <- min(y)

low_end <- y_min + (0.4*(y_max-y_min))
high_end <- y_min + (0.6*(y_max-y_min))
left <- x_min + (0.99*(x_max-x_min))
width <- (0.01*(x_max-x_min))
right <- left + width
height_increment <- (0.2*(y_max-y_min))/1000

cols <- grey_scale(seq(0,1,1/999))
for (i in 1:1000){
  bottom <- low_end + ((i-1)*height_increment)
  top <- bottom + height_increment
  polygon(x=c(left, left, right, right), y=c(bottom, top, top, bottom), col=cols[i], border=NA)
}

text(x=((left+right)/2), y=(high_end+(200*height_increment)), labels=expression(bold("KDE")), cex=0.6)
lines(x=c((left-width), (right+width)), y=c(low_end, low_end), lwd=0.5)
lines(x=c((left-width), (right+width)), y=c(high_end, high_end), lwd=0.5)
text(x=(left-(4*width)), y=low_end, labels="Low", cex=0.6)
text(x=(left-(4*width)), y=high_end, labels="High", cex=0.6)
######################


x_under <- c(-2, 5, 5)
y_under <- c(((-2)-cutoff), (5-cutoff), ((-2)-cutoff))
polygon(x=x_under, y=y_under, border="blue", lwd=1.5)

x_over <- c(-2, -2, 5)
y_over <- c(((-2)+cutoff), (5+cutoff), (5+cutoff))
polygon(x=x_over, y=y_over, border="gold", lwd=1.5)



abline(a=0, b=1, col="green", lwd=1.5)
abline(a=(-cutoff), b=1, col="blue", lwd=1.5)
abline(a=cutoff, b=1, col="gold", lwd=1.5)


legend("bottomright", legend=c(expression("y = x" + "Median( "*abs(""~epsilon~"")*" )"),
                               "y = x",
                               expression("y = x" - "Median( "*abs(""~epsilon~"")*" )")),
       box.col=transparent("white",0.2),
       bg=transparent("white",0.2),
       box.lwd = 0,
       lty=c(1,1,1),
       col=c("gold", "green", "blue"),
       seg.len=0.9, cex=0.9, lwd=1.5)
legend("topleft", legend=expression(italic("A. thaliana")), pch=NA, bty="n", text.col="grey40")


box()
fig_label(expression(bold("A")), cex=fig_lab_cex) 



par(mar=c(1,0,0,0.5))

plot(x=c(30:60), y=c(rep(-20,15),rep(120,16)), type="n", axes=F, xlab="", ylab="")

oID <- read.table(paste0("../results/nemo/",org,"/masked_graphpart/IDs/prediction/overpredicted.lst"), 
                  header=F, sep="\t")[,1]
wID <- read.table(paste0("../results/nemo/",org,"/masked_graphpart/IDs/prediction/well-predicted.lst"), 
                  header=F, sep="\t")[,1]
uID <- read.table(paste0("../results/nemo/",org,"/masked_graphpart/IDs/prediction/underpredicted.lst"), 
                  header=F, sep="\t")[,1]


ALLw <- ((length(wID) + length(uID)) / (length(oID) + length(wID) + length(uID))) * 100

ALLu <- (length(uID) / (length(oID) + length(wID) + length(uID))) * 100



#all
polygon(x=c(40,40,50,50), y=c(1,100,100,1), col="lemonchiffon")
polygon(x=c(40,40,50,50), y=c(1,ALLw,ALLw,1), col="darkolivegreen1")
polygon(x=c(40,40,50,50), y=c(1,ALLu,ALLu,1), col="darkslategray1")

text(adj=c(0,0.5), x=51, y=(((100+ALLw)/2)), 
     labels=paste(as.character(round((100-ALLw),1)), "%", sep=""), cex=text_cex)
text(adj=c(0,0.5), x=51, y=(((ALLw+ALLu)/2)), 
     labels=paste(as.character(round((ALLw-ALLu),1)), "%", sep=""), cex=text_cex)
text(adj=c(0,0.5), x=51, y=(((ALLu)/2)), 
     labels=paste(as.character(round((ALLu),1)), "%", sep=""), cex=text_cex)


text(adj=0.5, x=45, y=110, 
     labels=c("All\nGenes"), 
     cex=text_cex)


polygon(x=c(30,30,60,60), y=c(-25,-2,-2,-25))
text(x=c(35,35,35), y=c(-6, -13, -20), 
     labels=c("Overpredicted", "Well-Predicted", "Underpredicted"),
     adj=0, cex=0.85)
points(x=32, y=-6, pch=22, col="black", bg="lemonchiffon", cex=1.5)
points(x=32, y=-13, pch=22, col="black", bg="darkolivegreen1", cex=1.5)
points(x=32, y=-20, pch=22, col="black", bg="darkslategray1", cex=1.5)


#######################################################################
#######################################################################
# B
#######################################################################
#######################################################################

par(mar=c(3,3,2,0.1), mgp=c(1.4,0.4,0), tck=-0.03, xpd=F)
d <- density(abs(err))
plot(d$x, d$y, type="l", cex.lab=lab_cex, cex.axis=ax_cex,
     ylab="Density", xlab=expression("Absolute Error ("~abs(""~epsilon~"")~")"), xlim=c(0,5), ylim=c(0,1.7))
polygon(c(d$x[d$x <= median(abs(err))], median(abs(err))), 
        c(d$y[d$x <= median(abs(err))], 0),
        col="darkolivegreen1", border="grey30")
polygon(c(median(abs(err)), d$x[d$x >= median(abs(err))]), 
        c(0, d$y[d$x >= median(abs(err))]),
        col="grey95", border="grey30")
lines(d$x, d$y, lwd=1.2)
text(x=1.8, y=1.2, labels="Well-Predicted",
     cex=0.75)
text(x=2.7, y=0.8, labels=bquote("Median(" ~ abs(""~epsilon~"") ~ ") = "*.(as.character(round(cutoff,3)))),
     cex=0.75)
text(x=2.7, y=0.3, labels="Poorly Predicted",
     cex=0.75)
lines(x=c(0.15,1.6), 
      y=c(1,1.17), 
      lty=2, col=transparent("grey30", 0.6))
lines(x=c(median(abs(err)),2.5), 
      y=c(0.5,0.77), 
      lty=2, col=transparent("grey30", 0.6))
lines(x=c(1,2.5), 
      y=c(0.15,0.27), 
      lty=2, col=transparent("grey30", 0.6))

legend("topright", legend=expression(italic("A. thaliana")), text.col="grey40", bty="n")
fig_label(expression(bold("B")), cex=fig_lab_cex) 



dev.off() 
