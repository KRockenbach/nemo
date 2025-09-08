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



organisms <- c("Bnapus", "Athaliana")


families <- c("MYBrelated", "Homeobox", "C2C2dof", "WRKY",
              "TCP", "bZIP", "MYB", "MADS", "bHLH", "Trihelix",
              "NAC", "G2like", "HSF", "AP2EREBP", "C2C2gata")

families_full <- c("MYB-related", "Homeobox", "C2C2-Dof", "WRKY",
                   "TCP", "bZIP", "MYB", "MADS", "bHLH", "Trihelix",
                   "NAC", "G2-like", "HSF", "AP2/EREBP", "C2C2-GATA")


fig_lab_cex=1.5


plot_heatmap <- function(heat_mat, organism, palette=NA){
  par(mar=c(3,0,0,3))
  plot(x=c(0:99), y=c(0:99), type="n", axes=F, xlab="", ylab="")
  l <- 100/(ncol(heat_mat))
  if (is.na(palette)){
    col_range <- scico(n=1000, palette="berlin")
  } else {
    col_range <- viridis(n=1000, option=palette)
  }
  for (r in 1:nrow(heat_mat)){
    for (c in 1:nrow(heat_mat)){
      bottom <- (r-1)*l
      top <- r*l
      left <- (c-1)*l
      right <- c*l
      if (is.na(palette)){
        scaled <- (heat_mat[(16-r),c] + 1)/2
      } else {
        scaled <- heat_mat[(16-r),c] 
      }
      polygon(y=c(bottom,top,top,bottom), x=c(left,left,right,right), col=col_range[1+(ceiling(scaled*999))])
      threshold <- 0.45
      if(abs(heat_mat[(16-r),c]) >= threshold){
        text_col="black"
      } else {
        text_col="white"
      }
      text(x=(left+(l/2)), y=(bottom+(l/2)), labels=as.character(round(heat_mat[(16-r),c],2)), col=text_col, cex=0.35)
    }
  }
  axis(1, at=(l/2 + c(0:14)*l), labels=families_full, tick=F, las=2, line=-1.1, cex.axis=0.5)
  axis(4, at=(l/2 + c(0:14)*l), labels=rev(families_full), tick=F, las=2, line=-1.1, cex.axis=0.5)
}


plot_heat_key <- function(){
  plot(x=c(1:100), y=c(1:100), type="n", axes=F, xlab="", ylab="")
  col_range <- scico(n=1000, palette="berlin")
  for(i in 1:1000){
    l <- 50/1000
    bottom <- 25+((i-1)*l)
    top <- 25+(i*l)
    left <- 30
    right <- 40
    polygon(x=c(left,left,right,right), y=c(bottom,top,top,bottom), border=NA, col=col_range[i])
  }
  text(x=36.5, y=78, labels=expression(bold(rho["S"])), cex=1.1, col="grey40")
  lines(x=c(30, 41), y=c(75,75), col="grey40")
  lines(x=c(30, 41), y=c(50,50), col="grey40")
  lines(x=c(30, 41), y=c(25,25), col="grey40")
  text(x=43, y=75, adj=-1, labels="1", cex=0.8, col="grey40")
  text(x=43, y=50, adj=-1, labels="0", cex=0.8, col="grey40")
  text(x=43, y=25, adj=-1, labels="-1", cex=0.8, col="grey40")
}


plot_clust <- function(D, organism){
  C <- hclust(D, method = "median") # WPGMC
  par(xpd=T)
  swap <- C["order"]
  print(as.matrix(D))
  heat_mat <- 1-(as.matrix(D)[swap[[1]],swap[[1]]]) # 1 minus distance
  print(heat_mat)
  par(mar=c(3,0,0,3))
  plot(x=c(0:99), y=c(0:99), type="n", axes=F, xlab="", ylab="")
  l <- 100/(ncol(heat_mat))
  col_range <- magma(1000)
  for (r in 1:nrow(heat_mat)){
    for (c in 1:nrow(heat_mat)){
      bottom <- (r-1)*l
      top <- r*l
      left <- (c-1)*l
      right <- c*l
      polygon(y=c(bottom,top,top,bottom), x=c(left,left,right,right), col=col_range[1+(ceiling(heat_mat[r,(16-c)]*999))])
      threshold <- 0.45
      if(heat_mat[r,(16-c)] >= threshold){
        text_col="black"
      } else {
        text_col="white"
      }
      text(x=(left+(l/2)), y=(bottom+(l/2)), labels=as.character(round(heat_mat[r,(16-c)],2)), col=text_col, cex=0.35)
    }
  }
  axis(1, at=(l/2 + c(0:14)*l), labels=rev(C$labels[C$order]), tick=F, las=2, line=-1.1, cex.axis=0.5)
  axis(4, at=(l/2 + c(0:14)*l), labels=C$labels[C$order], tick=F, las=2, line=-1.1, cex.axis=0.5)
  
}

plot_clust_key <- function(){
  plot(x=c(1:100), y=c(1:100), type="n", axes=F, xlab="", ylab="")
  col_range <- magma(1000)
  for(i in 1:1000){
    l <- 50/1000
    bottom <- 25+((i-1)*l)
    top <- 25+(i*l)
    left <- 30
    right <- 40
    polygon(x=c(left,left,right,right), y=c(bottom,top,top,bottom), border=NA, col=col_range[i])
  }
  text(x=27, y=82, labels=expression(bold("TFBS\nCo-occurence\n(Szymkiewicz-\nSimpson\nCoefficient)")), cex=0.7, col="grey40", adj=0.3)
  lines(x=c(30, 41), y=c(75,75), col="grey40")
  lines(x=c(30, 41), y=c(50,50), col="grey40")
  lines(x=c(30, 41), y=c(25,25), col="grey40")
  text(x=43, y=75, adj=-1, labels="1", cex=0.8, col="grey40")
  text(x=43, y=50, adj=-1, labels="0.5", cex=0.8, col="grey40")
  text(x=43, y=25, adj=-1, labels="0", cex=0.8, col="grey40")
}


plot_dendro <- function(D, horiz=F, reverse=F){
  C <- hclust(D, method = "average") # UPGMA
  hcd <- as.dendrogram(C)
  # Default plot
  if (reverse){
    plot(rev(hcd), type = "rectangle", ylab = "", axes=F, horiz=horiz, edge.root=F, leaflab="none")
  } else {
    plot(hcd, type = "rectangle", ylab = "", axes=F, horiz=horiz, edge.root=F, leaflab="none") 
  }
}



png(paste0("fig_S9.png"), height=15, width=17, res=1200, units="cm")
layout(matrix(c(8,6,0,12,10,14,
                7,5,0,11,9,14,
                4,2,0,0,0,0,
                3,1,13,0,0,0
                ), ncol=6, byrow=T),
       heights=c(0.5,4,0.5,4),
       widths=c(0.5,4,1,0.5,4,1))

lab_cex=0.75
ax_cex=0.6
main_cex=0.85
fig_lab_cex=1


###################################

organism <- "Athaliana"
genes <- c()
class <- c()
for (f in 1:length(families)){
  TF <- families[f]
  TF_name <- families_full[f]
  df <- read.table(paste0("../data/Athaliana/parent_data/TF_ids/",TF,".all.narrowPeak.n.intersect.ids"), header=F)
  genes <- c(genes,df[,1])
  class <- c(class,rep(TF_name, nrow(df)))
}
df <- data.frame(genes, class)


# calculate distance matrix with Szymkiewicz-Simpson coefficient
S <- matrix(ncol=length(families), nrow=length(families))
colnames(S) <- families_full
rownames(S) <- families_full
for (fA in 1:length(families)){
  famA <- families_full[fA]
  for (fB in 1:length(families)){
    famB <- families_full[fB]
    A <- df$genes[df$class==famA]
    B <- df$genes[df$class==famB]
    inter_len <- sum(A %in% B)
    print(length(A))
    print(length(B))
    if (length(A) > length(B)){
      SSC <- inter_len/length(B)
    } else {
      SSC <- inter_len/length(A)
    }
    if (is.na(SSC)){
      SSC <- 0
    }
    S[fA,fB] <- SSC
  }
}
cooc_mat <- S


plot_heatmap(cooc_mat, organism = "Athaliana", palette="inferno")
par(mar=c(0,0,0.5,3))
plot(c(1:100),c(1:100), type="n", axes=F, ylab="", xlab="")
par(mar=c(3,0.5,0,0))
plot(c(1:100),c(1:100), type="n", axes=F, ylab="", xlab="")
if (organism=="Athaliana"){
  par(mar=c(0,0,0,0))
  plot(c(1:100),c(1:100), type="n", axes=F, ylab="", xlab="")
  fig_label(bquote(atop(bold("C"),"")), cex=fig_lab_cex, col="black")
  fig_label(bquote(atop("",italic("A. thaliana"))), cex=0.7, col="grey40")
}
####################################
  

for (organism in c("Bnapus","Athaliana")){
  imp <- read.table(paste0("../results/nemo/",organism,
                           "/masked_graphpart/attribs/promoter_TF_importance.tsv"),
                    header=T)
  imp <- imp[501:1501,]
  imp <- imp[,c("position","mean", "sd", families)]
  ## savitzky-golay filtering
  for (c in 2:ncol(imp)){
    R=1
    N=75
    P=1
    imp[,c] <- sgolayfilt(imp[,c], p = P, n = N)
  }
  for (f in families){
    imp[,f] <- mean_norm(imp[,f])
  }
  imp$mean <- mean_norm(imp$mean)
  
  ## standardize
  for (f in families){
    imp[,f] <- (imp[,f] - imp$mean)
  }
  for (f in families){
    imp[,f] <- (imp[,f]/imp$sd)
  }
  
  #####################
  
  dap <- data.table::transpose(read.table("../data/plot_data/DAPseq/DAPSeq.TF.promoter.tsv", header=F), make.names = 1)
  dap <- dap[,families]
  position <- as.numeric(rownames(dap))
  dap <- cbind(position, dap)
  cor_mat <- matrix(ncol=15, nrow=15)
  for (c in 1:15){
    for (r in 1:15){
      cor_mat[r,c] <- cor(dap[,(r+1)], imp[,(c+3)], method="spearman") 
    }
  }
  colnames(cor_mat) <- rownames(cor_mat) <- families_full
  plot_heatmap(heat_mat=cor_mat, organism)
  par(mar=c(0,0,0,3))
  plot(c(1:100),c(1:100), type="n", axes=F, ylab="", xlab="")
  text(x=50, y=50, labels="Normalized Importance", adj=0.5)
  par(mar=c(3,0,0,0))
  plot(c(1:100),c(1:100), type="n", axes=F, ylab="", xlab="")
  text(x=50, y=50, labels="DAP Signal", adj=0.5, srt=90)
  par(mar=c(0,0,0,0))
  plot(c(1:100),c(1:100), type="n", axes=F, ylab="", xlab="")
  if (organism=="Athaliana"){
    fig_label(bquote(atop(bold("B"),"")), cex=fig_lab_cex, col="black")
    fig_label(bquote(atop("",italic("A. thaliana"))), cex=0.8, col="grey40")
  } else {
    fig_label(bquote(atop(bold("A"),"")), cex=fig_lab_cex, col="black")
    fig_label(bquote(atop("",italic("B. napus"))), cex=0.8, col="grey40")
  }
}  

plot_clust_key()
plot_heat_key()

dev.off()


