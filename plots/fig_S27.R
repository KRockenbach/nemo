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


png("supp_figs/fig_S27.png", width=17, height=10, units="cm", res=1200)
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
  df <- read.table(paste0("../result_subset/nemo/",org,"/masked_graphpart/nemo",init,"_preds/concat_preds.tsv"), header=T, sep="\t")
  err <- (df$Predicted - df$Actual)
  exp <- read.table(paste0("../data_subset/", org, "/derived_data/expr_matrix.tsv"), header=T, sep="\t", row.names=1)
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
  df <- read.table(paste0("../result_subset/nemo/",org,"/masked_graphpart/nemo",init,"_preds/concat_preds.tsv"), header=T, sep="\t")
  err <- (df$Predicted - df$Actual)
  exp <- read.table(paste0("../data_subset/", org, "/derived_data/expr_matrix.tsv"), header=T, sep="\t", row.names=1)
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

