source("./plot_utils.R")

fig_lab_cex=1.5

zs11_preds <- read.table("../validation/ZS11_predictions.tsv", sep="\t", header=T)
colnames(zs11_preds) <- c("Gene", "ZS11_preds")
expr_preds <- read.table("../validation/Express617_predictions.tsv", sep="\t", header=T)
colnames(expr_preds) <- c("Gene", "E617_preds")
df <- merge(expr_preds, zs11_preds, by="Gene")

zs11_TPM <- read.table("../validation/ZS11_EAGLE_TPM.genes.tsv", sep="\t", header=T)
colnames(zs11_TPM) <- c("Gene", "ZS11_Leaf_rep1","ZS11_Leaf_rep2","ZSA3","ZSA3_2","ZSA3_3")
ZS11_Woolfenden_samples <- c("ZS11_Leaf_rep1","ZS11_Leaf_rep2")
ZS11_Zhang_samples <- c("ZSA3","ZSA3_2","ZSA3_3")
expr_TPM <- read.table("../validation/Express_EAGLE_TPM.genes.tsv", sep="\t", header=T)
colnames(expr_TPM) <- c("Gene", "E617_Leaf_rep1", "E617_Leaf_rep2", "E617_ExL21.TPM.lst")
E617_Woolfenden_samples <- c("E617_Leaf_rep1", "E617_Leaf_rep2", "E617_ExL21.TPM.lst")
df <- merge(df, zs11_TPM, by="Gene")
df <- merge(df, expr_TPM, by="Gene")

# exclude non-expressed genes (0 across the board)
df <- df[apply(df[c(ZS11_Woolfenden_samples,E617_Woolfenden_samples,ZS11_Zhang_samples)],MARGIN=1, FUN=sum) > 0,]


median_zs11_w <- as.vector(apply(df[,ZS11_Woolfenden_samples], MARGIN=1, FUN=median))
median_e617_w <- as.vector(apply(df[,E617_Woolfenden_samples], MARGIN=1, FUN=median))
median_zs11_z <- as.vector(apply(df[,ZS11_Zhang_samples], MARGIN=1, FUN=median))


median_zs11_w <- log10(median_zs11_w + 0.1)
median_zs11_z <- log10(median_zs11_z + 0.1)
median_e617_w <- log10(median_e617_w + 0.1)


plot_scatter <- function(x,y,genotype="ZS11", plotlabel=expression(bold("C"))){
  text_cex=0.7
  lab_cex=0.8
  df <- data.frame("Actual"=x,"Predicted"=y)
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
  text(x=3.91, y=2, labels="Low", adj=1, cex=text_cex, col="grey50")
  lines(x=c(3.95,4.1), y=c(3,3))
  text(x=3.91, y=3, labels="High", adj=1, cex=text_cex, col="grey50")
  text(x=4.1, y=3.6, labels="Kernel\nDensity\nEstimation", cex=text_cex, adj=1, col="grey40")
  box()
  axis(side=1, cex.axis=0.7, tck=-0.02)
  axis(side=2, cex.axis=0.7, tck=-0.02)
  
  text(x=-1, y=3.25, labels=genotype, cex=0.7, adj=0)
  
  r <- cor(x=df$Actual, y=df$Predicted, method="pearson")
  rsq <- r^2
  print(rsq)
  n.lab <- paste("N = ", as.character(length(df$Actual)), sep="")
  text(x=c(-1), y=c(3.8), labels=bquote("r"^2*" = "~.(format(round(rsq,3), nsmall=3))), 
       col="black", cex=text_cex, adj=0)
  text(x=3.6, y=-1, labels=n.lab, 
       col="grey70", cex=text_cex)
  
  
  #  2C upper margin
  dx <- density(df$Actual, cut=F)
  par(mar=c(0,3,0,0))
  plot(x=c(-1,dx$x), y=c(0,dx$y)-0.1, type="l", ylab="", xlab="", axes=F, xlim=c(-1,4))
  polygon(x=c(-1,dx$x), y=c(0,dx$y)-0.1, col="grey90", border=NA)
  lines(x=c(-1,dx$x), y=c(0,dx$y)-0.1)
  fig_label(plotlabel, cex=fig_lab_cex)
  
  
  # 2C right margin
  dy <- density(df$Predicted, cut=F)
  par(mar=c(3,0,0,0))
  plot(x=c(0,dy$y)-0.1, y=c(min(df$Predicted),dy$x), type="l", xlab="", ylab="", axes=F, ylim=c(min(dy$x),4))
  polygon(x=c(0,dy$y)-0.1, y=c(min(dy$x),dy$x), col="grey90", border=NA)
  lines(x=c(0,dy$y)-0.1, y=c(min(df$Predicted),dy$x))
}







png("supp_figs/fig_S8.png", width=17, height=15, units="cm", res=1000)

lab_cex=1
ax_cex=0.9
text_cex=1.1
fig_lab_cex=1.5

layout(
  matrix(c(2,0,5,0,
           1,3,4,6,
           0,0,8,0,
           0,0,7,9), ncol=4, byrow=TRUE), 
  widths=c(3,0.3,3,0.3), 
  heights=c(0.3,3,0.3,3)
)


plot_scatter(median_e617_w,df$E617_preds,genotype="Express617 leaves\nWoolfenden et al.", plotlabel=expression(bold("A")))
plot_scatter(median_zs11_w,df$ZS11_preds,genotype="ZS11 leaves\nWoolfenden et al.", plotlabel=expression(bold("B")))
plot_scatter(median_zs11_z,df$ZS11_preds,genotype="ZS11 leaves\nZhang et al.", plotlabel=expression(bold("C")))

dev.off()

