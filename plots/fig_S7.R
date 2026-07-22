

seq <- c("prom_up", "prom_down", "term_up", "term_down")
rsq_df <- as.data.frame(matrix(ncol=4, nrow=10)) 
colnames(rsq_df) <- seq

for (s in seq){
  column <- c()
  for (t in c(0:9)){
    df <- read.table(paste0("../result_subset/nemo/Bnapus/masked_graphpart/randomization/rand_",s,"_predictions.t_",as.character(t),".txt"), header=T, sep="\t")
    r <- cor(df$Actual_Median, df$Predicted_Median, method="pearson")
    column <- c(column, r**2)
  }
  rsq_df[,s] <- column
}

png("supp_figs/fig_S7.png", width=17, height=17, res=800, unit="cm")
par(mar=c(5,5,2,0.5), xpd=T)
cols=turbo(4)
boxplot(rsq_df, at=c(0.6,1.6,3.4,4.4), axes=F, xlab="Randomized Input Region", ylab=expression("Model Performance (r²)"), cex.lab=1.2, xlim=c(-0.5,5.5), col=cols, border="grey40")
axis(2)
text(x=c(-0.2,2.4,2.6,5.2), y=apply(rsq_df, MARGIN = 2, FUN = median), 
     labels = as.character(round(apply(rsq_df, MARGIN = 2, FUN = median),2)))
text(x=c(0.5,1.7,3.3,4.5), y=0.02, labels = c("Upstream\nPromoter",
                                              "Downstream\nPromoter",
                                              "Upstream\nTerminator",
                                              "Downstream\nTerminator"),
     cex=0.8)
lines(x=c(-0.5,2.1), y=rep(0.52,2), lty=3, lwd=1.5, col=cols[1])
lines(x=c(2.9,5.5), y=rep(0.52,2), lty=3, lwd=1.5, col=cols[4])

lines(x=c(-0.2,1.1), y=rep(0.52,2), lty=1, lwd=3, col=cols[1])
lines(x=c(1.1,2.1), y=rep(0.52,2), lty=1, lwd=3, col=cols[2])

lines(x=c(2.9,3.9), y=rep(0.52,2), lty=1, lwd=3, col=cols[3])
lines(x=c(3.9,5.2), y=rep(0.52,2), lty=1, lwd=3, col=cols[4])


lines(x=c(1.1,1.1), y=c(0.52,0.53), lty=1, lwd=2)
lines(x=c(3.9,3.9), y=c(0.52,0.53), lty=1, lwd=2)

lines(x=c(1.1,1.2), y=c(0.53,0.53), lty=1, lwd=2)
lines(x=c(3.84,3.96), y=c(0.53,0.53), lty=1, lwd=2)

polygon(x=c(1.2,1.2,1.23), y=c(0.528,0.532,0.53), col="black")

text(x=c(1.1,3.9), y=0.54, labels=c("TSS", "TTS"), cex=1)


dev.off()
