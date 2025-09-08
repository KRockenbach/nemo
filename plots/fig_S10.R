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


png("fig_S10.png", height=20, width=17, units="cm", res=1200)
par(mfrow=c(2,1), mar=c(4.2,4.2,0.1,0.1))

for (organism in c("Bnapus", "Athaliana")){
  idf <- read.table(paste0("../data/plot_data/", organism, "_promoter_TF_inserted.tsv"), header=T, sep="\t")
  
  
  auc <- function(x){
    x <- x[501:1501]
    area <- 0
    for (i in 1:length(x)){
      if (i == 1){
        prev <- x[i]
        next
      } else {
        current <- x[i]
        # is previous > current, third term becomes negative
        # negative values result in negarive area
        area <- area + prev + ((current-prev)/2)
      }
    }
    return(area)
  }
  
  idf <- idf[,-1]
  area <- apply(as.matrix(idf), MARGIN=2, FUN=auc)
  area <- area[order(names(area))]
  
  initial <- strsplit(organism, split="")[[1]][1]
  exp_path <- paste0("../results/nemo/", organism, "/masked_graphpart/nemo", initial, "_preds/concat_preds.tsv")
  exp_df <- read.table(exp_path, header=T, sep="\t")
  
  families <- c("MYBrelated", "Homeobox", "C2C2dof", "WRKY",
                "TCP", "bZIP", "MYB", "MADS", "bHLH", "Trihelix",
                "NAC", "G2like", "HSF", "AP2EREBP", "C2C2gata")
  upstream <- c("MYBrelated", "Homeobox", "C2C2dof", "WRKY",
                "TCP", "bZIP", "MYB", "MADS")
  both <- c("bHLH", "Trihelix", "NAC", "G2like")
  downstream <- c("HSF", "AP2EREBP", "C2C2gata")
  
  
  exp <- c()
  index <- 1
  for (f in families){
    if (organism == "Bnapus"){
      ID_path <- paste0("../data/Bnapus/parent_data/TF_ids/",f,".at.cov.ids")
    } else {
      ID_path <- paste0("../data/Athaliana/parent_data/TF_ids/",f,".all.narrowPeak.n.intersect.ids")
    }
    IDs <- read.table(ID_path, sep="\t", header=F)[,1]
    for (i in 1:length(IDs)){
      IDs[i] <- strsplit(IDs[i], split=".", fixed=T)[[1]][1]
    }
    exp[index] <- mean(exp_df$Actual[exp_df$ID %in% IDs])
    index <- index + 1
  }
  
  names(exp) <- families
  exp <- exp[order(names(exp))]
  
  
  cols <- c("brown","red","tomato", "tan1", "gold", "khaki", "lawngreen",
            "aquamarine", "seagreen", "cornflowerblue", "blue", "darkslategray", "plum", "purple", "magenta")

  names(cols) <- families
  
  plot(exp, area, col=cols, pch=19, cex=1, type="n",
       ylab=expression("AUC( "*Delta*"E"*" )"),
       xlab="Average Expression Level")
  
  for (f in families){
    if (f %in% upstream){
      char <- 24
    } else if (f %in% downstream){
      char <- 25
    } else {
      char <- 23
    }
    points(exp[f], area[f], col="black", bg=cols[f], pch=char, cex=1)
  }
  
  legend("topleft", 
         legend=c(expression(bold("Family")), families[1:8], expression(bold("TF binding (DAP-Seq)")), "mostly upstream", "mostly downstream", 
                             NA, families[9:15], NA, NA, "both", NA), 
         pch=c(NA, rep(19,8), NA, 24, 25, NA, rep(19,7), NA, NA, 23, NA), 
         col=c(NA, cols[families[1:8]], NA, NA, NA, NA, cols[families[9:15]], NA, NA, NA, NA),
         pt.bg="black", ncol=2, pt.cex=1, cex=0.65, bty="n")
  
  abline(lm(area~exp))
  print(cor.test(area,exp, method="pearson"))
  ct <- cor.test(area,exp, method="pearson")
  r <- ct$estimate
  p <- ct$p.value
  
  if (organism == "Bnapus"){
    org_lab <- "B. napus"
  } else {
    org_lab <- "A. thaliana"
  }
  legend("bottomright", bty="n", legend=bquote(italic(.(org_lab))), text.col="grey40")
  Lines <- list(bquote(rho["p"]*" = "*.(round(r,3))),bquote("p-value = "*.(round(p,3))))
  legend("top", bty="n", legend=do.call(expression,Lines)) 

}


dev.off()

