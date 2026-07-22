png("supp_figs/fig_S6.png", width=17, height=10, units="cm", res=1200)
par(mar=c(4,4.5,0.2,0.1))

param_df <- read.table("../data_subset/plot_data/model_params.tsv", header=F, sep="\t")
colnames(param_df) <- c("modelname", "num_params")

# 2A
df = read.table("../result_subset/rsq_df.tsv", header=T, sep="\t")
sub=df[df$partitioning=="graphpart_Bn",]
sub$model[is.na(sub$valid_fold)] <- "nemo90"
model_names <- param_df$modelname
nm <- length(model_names)
medium_rsq <- c()
for (m in 1:nm){
  model <- model_names[m]
  rsq <- sub$rsq[sub$model==model & sub$tpm_type=="median"]
  medium_rsq <- c(medium_rsq,round(median(rsq),2))
}

param_df <- cbind(param_df, medium_rsq)

plot((param_df$num_params)/1e6, param_df$medium_rsq,
     ylab=expression("Model Performance (r²)"),
     xlab="Million Model Parameters",
     pch=16, cex=1, bty="n", ylim=c(0.35, 0.55), xlim=c(0,25), cex.axis=0.8, cex.lab=0.8)
text(x=((param_df$num_params)/1e6), y=((param_df$medium_rsq)+0.01), 
     labels=c("Basenji-5K", "Xpresso", "Xpresso without halflife features", 
              expression(italic(n)*"emo"), "PhytoExpr-ensemble", "PhytoExpr-transformer"), 
     adj=0.1, cex=0.5)

dev.off()
