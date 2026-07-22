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


# processing steps:
########################
# importance:
# 1) importance for each gene is calculated as the sum of absolute attributions across bases at each position
# 2) global mean and sd across all genes is calculated at each position (taking masked segments into account)
# 3) mean across each group of genes is calculated (taking masked segments into account)
# 4) savitzky-golay filter is applied to global mean, global sd, and group means
# 5) global mean is mean-normalized (mean across both sequences is subtracted from each position)
# 6) group means are mean-normalized (mean across both sequences is subtracted from each position)
# 7) normalized global mean is subtracted from normalized group means
# 8) difference is divided by global standard deviation
# 9) color values are scaled globally for each sequence 
    #(min-max scaling to the range [0,1], 
    #where the 0 represents the minumum across groups and positions 
    #and 1 represents the maximum across groups and positions)


# DAP-seq:
# Min-max scaling to range [0,1] within each group (by row) across both sequences! 



organisms <- c("Bnapus", "Athaliana")

families <- c("MYBrelated", "Homeobox", "C2C2dof", "WRKY",
              "TCP", "bZIP", "MYB", "MADS", "bHLH", "Trihelix",
              "NAC", "G2like", "HSF", "AP2EREBP", "C2C2gata")
families_full <- c("MYB-related", "Homeobox", "C2C2-Dof", "WRKY",
                   "TCP", "bZIP", "MYB", "MADS", "bHLH", "Trihelix",
                   "NAC", "G2-like", "HSF", "AP2/EREBP", "C2C2-GATA")


fig_lab_cex=1.5



plot_info <- function(organism, numbers=NULL, descriptor=NULL, label=expression(bold("A")), text_cex=0.6, descriptor_height=110){
  par(mar=c(0.1,0,0.8,0), tck=-0.03, mgp=c(1.5,0.4,0))
  plot(x=c(-460:460), y=c(0:-920), type="n", 
       axes=F, xlab="", ylab="", ylim=c(-1500,0))
  fig_label(label, cex=fig_lab_cex)
  if (organism=="Athaliana"){
    text(x=-300, y=-750, labels=expression(italic("A. thaliana")), adj=0, cex=1, col="grey40", srt=90)
  } else {
    text(x=-300, y=-750, labels=expression(italic("B. napus")), adj=0, cex=1, col="grey40", srt=90)
  }
  if (!is.null(descriptor)){
    text(x=250, y=descriptor_height, labels=descriptor, adj=0, cex=0.6, col="grey50") 
  }
  for (f in 1:length(families)){
    y1 <- (-1)*(f-1)*100
    y2 <- (-1)*(f)*100
    if (!is.null(numbers)){
      text(x=250, y=(y1+y2)/2, labels=as.character(numbers[f]), adj=0, cex=text_cex)
      text(x=200, y=110, labels="Family", adj=1, cex=0.6, col="grey50")
      text(x=200, y=(y1+y2)/2, labels=families_full[f], adj=1, cex=text_cex)
    } else {
      text(x=400, y=110, labels="Family", adj=1, cex=0.6, col="grey50")
      text(x=400, y=(y1+y2)/2, labels=families_full[f], adj=1, cex=text_cex)
    }
  }
}


plot_legend <- function(palette="inferno", descriptor="DAP-Seq Peak Coverage\nScaled to [0,1]\nWithin Family", 
                        legend_labs=c("0","1"), middle="0.5"){
  par(mar=c(1.5,0.5,1.5,0.5))
  plot(x=c(-60:60), y=c(0:120), type="n", 
       axes=F, xlab="", ylab="", ylim=c(0,120), cex.lab=lab_cex)
  if (palette=="custom"){
    cols <- get_colors(seq(-1,1,2/999), palette="custom") 
  } else {
    cols <- get_colors(seq(0,1,1/999), palette=palette)
  }
  #legend width
  lw <- 15
  #legend height
  lh <- 0.05
  for (i in 1:1000){
    polygon(x=((c(-1,-1,1,1)*lw)-5), y=c(25+((i-1)*lh), 25+(i*lh), 25+(i*lh), 25+((i-1)*lh)), col=cols[i], border=NA)
  }
  lines(x=(c(-lw-1,lw+1)-5), y=c(25,25))
  lines(x=(c(-lw-1,lw+1)-5), y=c(75,75))
  text(x=c(lw,lw), y=c(25,75), labels=legend_labs, cex=0.5, adj=0, col="grey40")
  lines(x=(c(-lw-1,-lw+3)-5), y=c(50,50), col="black")
  lines(x=(c(lw-3,lw+1)-5), y=c(50,50), col="black")
  text(x=lw, y=50, labels=middle, cex=0.5, adj=0, col="grey40")
  par(xpd=T)
  text(x=0, y=110, labels=descriptor, cex=0.55)
}

plot_header <- function(seq="Promoter"){
  if (seq == "Promoter"){
    par(mar=c(0.1,0.1,0.3,0.3), cex.axis=0.6, xpd=T)
  } else {
    par(mar=c(0.1,0.3,0.3,0.1), cex.axis=0.6, xpd=T)
  }
  plot(x=seq(-500,500,100), y=c(0:10), type="n", axes=F, xlab="", ylab="")
  text(x=0, y=5, labels=bquote(bold(.(seq))), cex=main_cex)
}


plot_heatmap <- function(df, seq="Promoter", palette="inferno", midline_col="white", border_col=transparent("grey66", alpha=0.3), header=F, axis=F){
  if (seq == "Promoter"){
    ref_point <- "TSS"
  } else {
    ref_point <- "TTS"
  }
  if (seq == "Promoter"){
    par(mar=c(0.1,0.1,0.8,0.3), cex.axis=0.6, xpd=T)
  } else {
    par(mar=c(0.1,0.3,0.8,0.1), cex.axis=0.6, xpd=T)
  }
  plot(x=c(-460:460), y=c(0:-920), type="n", 
       axes=F, xlab="", ylab="", ylim=c(-1500,0), cex.lab=lab_cex)
  
  ###############
  for (f in 1:15){
    fam <- families[f]
    y1 <- (-1)*(f-1)*100
    y2 <- (-1)*f*100
    cols <- get_colors(df[,fam], palette=palette)
    step <- round((1000/nrow(df)),0)
    if (seq == "Promoter"){indices = c(1:nrow(df))}else{indices = c(nrow(df):1)} # reverse plotting order for terminator, so TTS part of upstream
    for (idx in indices){
      x1 <- (-500) + ((idx-1)*step)
      x2 <- (-500) + (idx*step)
      polygon(x=c(x1,x1,x2,x2), y=c(y1,y2,y2,y1), border=NA, col=cols[idx])
    }
    for (segment_border in seq(-400,400,100)){
      lines(x=c(segment_border,segment_border), y=c(y1,y2), col=transparent(midline_col, alpha=0.5), lwd=0.5)
    }
    polygon(x=c(-500,-500,500,500), y=c(y1,y2,y2,y1), lwd=1, border=border_col)
  } 
  
  lines(x=c(0,0), y=c(-1500,0), col=transparent(midline_col, alpha=0.6), lwd=2)
  if (axis){
    axis(side=1, at=seq(-500,500,100), las=2, cex.axis=ax_cex)
  }
  if (header){
    title(main=seq, line=0.02, cex.main=main_cex) 
  }
}


plot_axlab <- function(seq="Promoter"){
  if (seq == "Promoter"){
    par(mar=c(0,0.1,2,0.3), cex.axis=0.6, xpd=T)
    ref_point="TSS"
  } else {
    par(mar=c(0,0.3,2,0.1), cex.axis=0.6, xpd=T)
    ref_point="TTS"
  }
  plot(x=seq(-500,500,100), y=c(0:10), type="n", axes=F, xlab="", ylab="")
  text(x=0, y=9, labels=paste0("Position Relative to ", ref_point), cex=lab_cex)
}

#############################################


pdf("fig_5.pdf", height=18, width=17, pointsize=30)
par(xpd=T)
layout(matrix(c(0,1,2,0,
                3,4,5,6,
                7,8,9,10,
                11,12,13,14,
                15,16,17,18,
                19,20,21,22,
                0,23,24,0), ncol=4, byrow=T),
       heights=c(0.2, rep(1,5), 0.5),
       widths=c(1.5,2.5,2.5,0.8))


seqs <- c("promoter", "terminator")


lab_cex=0.75
ax_cex=0.6
main_cex=1


ppath <- "../data_subset/plot_data/DAPseq/DAPSeq.TF.promoter.tsv"
tpath <- "../data_subset/plot_data/DAPseq/DAPSeq.TF.terminator.tsv"

# load and transpose df
load_and_transpose <- function(path){
  df <- read.table(path, header=F)
  col_vec <- df[,1]
  df <- df[,-1]
  df <- t(df)
  df <- as.data.frame(df)
  colnames(df) <- col_vec
  rownames(df) <- c(1:nrow(df))
  return(df)
}
pdf <- load_and_transpose(ppath)
tdf <- load_and_transpose(tpath)
tdf <- tdf[c(nrow(tdf):1),] # positions of terminator df are reversed

pdf <- pdf[,families]
tdf <- tdf[,families]

max_cov <- c()
for (f in families){
  max_cov <- c(max_cov, max(c(pdf[,f],tdf[,f])))
}

# scale data within families
for (f in families){
  norm.cat <- min_max.positive.scale(c(pdf[,f], tdf[,f]))
  pdf[,f] <- norm.cat[1:(length(norm.cat)/2)]
  tdf[,f] <- norm.cat[((length(norm.cat)/2)+1):length(norm.cat)]
}

for (f in families){
  R=1
  N=101
  P=1
  pdf[,f] <- sgolayfilt(pdf[,f], p = P, n = N)
  tdf[,f] <- sgolayfilt(tdf[,f], p = P, n = N)
}

plot_header("Promoter")
plot_header("Terminator")
plot_info("Athaliana", numbers=max_cov, descriptor="Max. Cov.", label=expression(bold("A")), descriptor_height=113)
plot_heatmap(pdf, seq="Promoter", palette="inferno", midline_col="white", border_col="black")
plot_heatmap(tdf, seq="Terminator", palette="inferno", midline_col="white", border_col="black")
plot_legend(palette="inferno", descriptor="DAP-Seq Peak Coverage\nScaled to [0,1]\nWithin Family", legend_labs=c("0","1"))
#######################################

for (organism in c("Bnapus", "Athaliana")){
  if (organism == "Athaliana"){
    all_genes <- read.table("../result_subset/nemo/Athaliana/masked_graphpart/attribs/gene_names.lst", header=F)
    counts <- c()
    for (f in families){
      group_genes <- read.table(paste0("../data_subset/Athaliana/parent_data/TF_ids/", f, ".all.narrowPeak.n.intersect.ids"), header=F)
      counts <- c(counts,sum(group_genes[,1] %in% all_genes[,1]))
    }
  } else {
    all_genes <- read.table("../result_subset/nemo/Bnapus/masked_graphpart/attribs/gene_names.lst", header=F)
    counts <- c()
    for (f in families){
      group_genes <- read.table(paste0("../data_subset/Bnapus/parent_data/TF_ids/", f, ".at.cov.ids"), header=F)
      for (i in 1:nrow(group_genes)){
        group_genes[i,1] <- str_replace(group_genes[i,1], pattern=fixed(".1"), replacement="") 
      }
      counts <- c(counts,sum(group_genes[,1] %in% all_genes[,1]))
    }
  }
  pdf <- read.table(paste0("../result_subset/nemo/",organism,
                           "/masked_graphpart/attribs/promoter_TF_importance.tsv"),
                    header=T)
  tdf <- read.table(paste0("../result_subset/nemo/",organism,
                           "/masked_graphpart/attribs/terminator_TF_importance.tsv"),
                    header=T)
  tdf$position <- tdf$position + 1 # shift position, because TTS is part of upstream sequence
  # trim to +- 500 bp
  pdf <- pdf[501:1501,]
  tdf <- tdf[501:1501,]
  pdf <- pdf[,c("position","mean", "sd", families)]
  tdf <- tdf[,c("position","mean", "sd", families)]
  
  # savitzky-golay filtering
  for (c in 2:ncol(pdf)){
    R=1
    N=101 #75
    P=1
    pdf[,c] <- sgolayfilt(pdf[,c], p = P, n = N)
    tdf[,c] <- sgolayfilt(tdf[,c], p = P, n = N)
  }
  # mean normalize
  for (f in families){
    concat <- c(pdf[,f], tdf[,f])
    concat <- mean_norm(concat) #!
    pdf[,f] <- concat[1:length(pdf[,f])]
    tdf[,f] <- concat[(length(pdf[,f])+1):length(concat)]
    
  }
  # mean normalize (across sequence) mean (across samples)
  concat_mean <- c(pdf$mean, tdf$mean)
  concat_mean <- mean_norm(concat_mean) #!
  pdf$mean <- concat_mean[1:length(pdf$mean)]
  tdf$mean <- concat_mean[(length(pdf$mean)+1):length(concat_mean)]
  
  
  # standardize
  for (f in families){
    pdf[,f] <- (pdf[,f] - pdf$mean)
    tdf[,f] <- (tdf[,f] - tdf$mean)
  }
  
  for (f in families){
    pdf[,f] <- (pdf[,f]/pdf$sd)
    tdf[,f] <- (tdf[,f]/tdf$sd)
  }
  
  # scale globally
  norm.cat <- min_max.positive.scale(as.matrix(rbind(pdf[,families], tdf[,families])))
  pdf[,families] <- norm.cat[1:(nrow(norm.cat)/2),]
  tdf[,families] <- norm.cat[((nrow(norm.cat)/2)+1):nrow(norm.cat),]
  
  if (organism=="Athaliana"){
    label=expression(bold("C"))
    plot_info("Athaliana", numbers=counts, descriptor="# Genes", label=label)
  } else {
    label=expression(bold("B"))
    plot_info("Bnapus", numbers=counts, descriptor="# Genes", label=label)
  }
  
  plot_heatmap(pdf, seq="Promoter", palette="inferno", midline_col="white", border_col="black")
  plot_heatmap(tdf, seq="Terminator", palette="inferno", midline_col="white", border_col="black")
  plot_legend(palette="inferno", descriptor="Normalized Importance\nGlobally Scaled to [0,1]", legend_labs=c("0","1"))
}


#######################################
for (o in 1:2){
    organism <- organisms[o]
    if (organism == "Athaliana"){
      axis=T
    } else {
      axis=F
    }

    df <- read.table(paste0("../result_subset/nemo/",organism,
                             "/masked_graphpart/modisco/seqlet_coordinates.tsv"),
                      header=T)

    bin_start <- seq(-500,400,100)
    pdf <- as.data.frame(bin_start)
    tdf <- as.data.frame(bin_start)
    for (family in families){
        fam_col = c()
        for (bs in bin_start){
          be <- bs + 100
          fam_col <- c(fam_col, nrow(df[df$seq=="promoter" &
                                        df$family==family &
                                        bs <= df$start &
                                        df$end <= be,])) # seqlet contained in bin
        }
        pdf <- cbind(pdf, fam_col)
        colnames(pdf)[ncol(pdf)] <- family
        fam_col = c()
        for (bs in bin_start){
          be <- bs + 100
          fam_col <- c(fam_col, nrow(df[df$seq=="terminator" &
                                          df$family==family &
                                          bs <= df$start &
                                          df$end <= be,]))
        }
        tdf <- cbind(tdf, fam_col)
        colnames(tdf)[ncol(tdf)] <- family
    }

    seqlet_num <- c()
    for (f in families){
      seqlet_num <- c(seqlet_num, sum(c(pdf[,f],tdf[,f])))
    }

    # scale data within families
    for (f in families){
      norm.cat <- min_max.positive.scale(c(pdf[,f], tdf[,f]))
      pdf[,f] <- norm.cat[1:(length(norm.cat)/2)]
      tdf[,f] <- norm.cat[((length(norm.cat)/2)+1):length(norm.cat)]
    }


    if (organism=="Athaliana"){
      label=expression(bold("E"))
    } else {
      label=expression(bold("D"))
    }

    plot_info(organism, numbers=seqlet_num, descriptor="# Seqlets", label=label)
    plot_heatmap(pdf, seq="Promoter", palette="inferno", midline_col="white", border_col=transparent("grey66", alpha=0.3), axis=axis)
    plot_heatmap(tdf, seq="Terminator", palette="inferno", midline_col="white", border_col=transparent("grey66", alpha=0.3), axis=axis)
    plot_legend(palette="inferno", descriptor="Number of Seqlets\nScaled to [0,1]\nWithin Family", legend_labs=c("0","1"))

}

plot_axlab("Promoter")
plot_axlab("Terminator")

dev.off()

