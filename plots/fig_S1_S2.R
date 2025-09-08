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

plot_clust <- function(D, organism){
  C <- hclust(D, method = "median") # WPGMC
  par(xpd=T)
  swap <- C["order"]
  heat_mat <- 1-(as.matrix(D)[swap[[1]],swap[[1]]]) # 1 minus distance
  par(mar=c(18,0,0,18))
  plot(x=c(0:99), y=c(0:99), type="n", axes=F, xlab="", ylab="")
  l <- 100/(ncol(heat_mat))
  col_range <- plasma(1000)
  for (r in 1:nrow(heat_mat)){
    for (c in 1:nrow(heat_mat)){
      bottom <- (r-1)*l
      top <- r*l
      left <- (c-1)*l
      right <- c*l
      polygon(y=c(bottom,top,top,bottom), x=c(left,left,right,right), col=col_range[1+(ceiling(heat_mat[r,(ncol(heat_mat)+1-c)]*999))], border=NA)
      threshold <- 0.45
      if(heat_mat[r,(ncol(heat_mat)+1-c)] >= threshold){
        text_col="black"
      } else {
        text_col="white"
      }
    }
  }
  axis(1, at=(l/2 + c(0:(ncol(heat_mat)-1))*l), labels=rev(C$labels[C$order]), tick=F, line=-1.1, cex.axis=0.4, las=2)
  axis(4, at=(l/2 + c(0:(ncol(heat_mat)-1))*l), labels=C$labels[C$order], tick=F, line=-1.1, cex.axis=0.4, las=2)
  
  for (i in 0:999){
    delta <- (200 - 110)/1000
    left=110+(i*delta)
    right=left+delta
    polygon(x=c(left,left,right,right), y=((-1)*c(40,80,80,40)), col=col_range[i+1], border=NA)
  }
  polygon(x=c(110,110,200,200), y=((-1)*c(40,80,80,40)))
  x <- density(heat_mat, from=0, to=1)$x
  x <- (x*90)+110
  y <- density(heat_mat, from=0, to=1)$y
  y <- (((y-min(y))/(max(y)-min(y)))*20)-80
  lines(x,y, col="cyan")
  x <- seq(from=-1, to=1, by=0.2)
  for (x_i in x){
    x_shift <- (((x_i+1)/2)*90)+110
    lines(x=c(x_shift, x_shift), y=c(-80,-81))
    text(x=x_shift, y=-83, labels=as.character(x_i), cex=0.5)
  }
  x_shift <- (((x+1)/2)*90)+110
  text(x=mean(x_shift), y=-86, labels=expression(rho["Spearman"]), cex=0.7)
  text(x=mean(x_shift), y=-35, labels=bquote(mu~"="~.(as.character(round(mean(heat_mat), 2)))), cex=0.6)
}

plot_clust_key <- function(){
  plot(x=c(1:100), y=c(1:100), type="n", axes=F, xlab="", ylab="")
  col_range <- plasma(1000)
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



for (organism in c("Bnapus", "Athaliana")){
  df <- read.table(paste0("../data/", organism, "/derived_data/expr_matrix.tsv"), header=T, row.names=1, sep="\t")
  colnames(df) <- str_replace_all(string=colnames(df), pattern=fixed(".."), replacement=("_"))
  df[is.na(df)] <- 0.0
  medians <- apply(as.matrix(df),MARGIN=1,FUN=median)
  df <- df[order(medians, decreasing = T)[1:5000],] # select 1000 genes with highest median expression
  mat <- matrix(ncol=ncol(df), nrow=ncol(df))
  for (col in 1:ncol(df)){
    for (row in 1:ncol(df)){
      mat[row,col] <- cor(df[,col], df[,row], method="spearman")
    }
  }
  mat[is.na(mat)] <- 0.0 # due to zero standard deviation
  colnames(mat) <- rownames(mat) <- colnames(df)
  if (organism == "Bnapus"){
    figname <- "fig_S1.png"
  } else {
    figname <- "fig_S2.png"
  }
  png(figname, width=17, height=17, res=1200, units="cm")

  D <- as.dist((1-mat)/2)
  
  plot_clust(D, organism)

  dev.off()
}

