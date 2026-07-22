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



plot_info <- function(organism, numbers=NULL, descriptor=NULL, label=expression(bold("A"))){
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
    text(x=250, y=110, labels=descriptor, adj=0, cex=0.55, col="grey50") 
  }
  for (f in 1:length(families)){
    y1 <- (-1)*(f-1)*100
    y2 <- (-1)*(f)*100
    if (!is.null(numbers)){
      text(x=250, y=(y1+y2)/2, labels=as.character(numbers[f]), adj=0, cex=0.6)
      text(x=200, y=110, labels="Family", adj=1, cex=0.55, col="grey50")
      text(x=200, y=(y1+y2)/2, labels=families_full[f], adj=1, cex=0.6)
    } else {
      text(x=400, y=110, labels="Family", adj=1, cex=0.55, col="grey50")
      text(x=400, y=(y1+y2)/2, labels=families_full[f], adj=1, cex=0.6)
    }
  }
}


plot_legend <- function(palette="inferno", descriptor="Relative DAP_seq\nPeak Coverage\n(Scaled within Family)", 
                        legend_labs=c("0","1"), middle="0.5", min_max=NULL){
  par(mar=c(1.5,0.5,1.5,0.5))
  plot(x=c(-60:60), y=c(0:120), type="n", 
       axes=F, xlab="", ylab="", ylim=c(0,120), cex.lab=lab_cex)
  if (palette=="custom"){
    if (is.null(min_max)){
      max <- 1
      min <- -1
    } else {
      max <- max(abs(min_max))
      min <- -max
    }
    cols <- get_colors(seq(min,max,(2*max)/999), palette="custom") 
  } else {
    cols <- get_colors(seq(0,1,1/999), palette=palette)
  }
  #legend width
  lw <- 15
  #legend height
  lh <- 0.05
  for (i in 1:1000){
    polygon(x=c(-1,-1,1,1)*lw, y=c(25+((i-1)*lh), 25+(i*lh), 25+(i*lh), 25+((i-1)*lh)), col=cols[i], border=NA)
  }
  lines(x=c(-lw-1,lw+1), y=c(25,25))
  lines(x=c(-lw-1,lw+1), y=c(75,75))
  text(x=c(0,0), y=c(18,82), labels=legend_labs, cex=0.55, adj=0.5, col="grey40")
  lines(x=c(-lw-1,-lw+3), y=c(50,50), col=transparent("grey20", alpha=0.5))
  lines(x=c(lw-3,lw+1), y=c(50,50), col=transparent("grey20", alpha=0.5))
  text(x=lw+5, y=50, labels=middle, cex=0.55, adj=0, col="grey40")
  par(xpd=T)
  text(x=0, y=100, labels=descriptor, cex=0.55)
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
      lines(x=c(segment_border,segment_border), y=c(y1,y2), col=transparent(midline_col, alpha=0.2), lwd=0.5)
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


get_overall_min_max <- function(){
  # overall boundaries used as reverence for scaling data across subplots in supplementary figure
  max_vec <- c()
  min_vec <- c()
  for (organism in c("Athaliana", "Bnapus")){
    for (family in families){
      df <- read.table(paste0("../result_subset/nemo/",organism,"/masked_graphpart/motif_insertion/", family, "_medium.tsv"), header=T)
      df$diff <- (df$mutated_pred - df$baseline)
      for (type in c("TFBS", "scrambled", "random", "reversed")){
        x <- c()
        for (position in c(-500:499)){
          x <- c(x, median(df$diff[df$motif_type==type & df$mutated_sequence_type=="promoter" & df$relative_idx==position]))
        }
        R=1
        N=101
        P=1
        x <- sgolayfilt(x, p = P, n = N)
        max_vec <- c(max_vec, max(x))
        min_vec <- c(min_vec, min(x))
      }
    }
  }
  return(c(min(min_vec), max(max_vec)))
}
overall_min_max <- get_overall_min_max()


get_mutation_data <- function(organism, type="TFBS", control=NULL){
  pdf <- data.frame(position=c(-500:499))
  tdf <- data.frame(position=c(-500:499))
  for (family in families){
    df <- read.table(paste0("../result_subset/nemo/",organism,"/masked_graphpart/motif_insertion/", family, "_medium.tsv"), header=T)
    if (!is.null(control)){
      df_contr <- df[df$motif_type==control,] ### 
    }
    df <- df[df$motif_type==type,]
    if (!is.null(control)){
      df$diff <- (df$mutated_pred - df_contr$mutated_pred) 
    } else {
      df$diff <- (df$mutated_pred - df$baseline) 
    }
    pdiff <- c()
    tdiff <- c()
    for (p in pdf$position){
      pdiff <- c(pdiff, median(df$diff[df$mutated_sequence_type=="promoter" & df$relative_idx==p]))
      tdiff <- c(tdiff, median(df$diff[df$mutated_sequence_type=="terminator" & df$relative_idx==p]))
    }
    pdf <- cbind(pdf, pdiff)
    colnames(pdf)[ncol(pdf)] <- family
    tdf <- cbind(tdf, tdiff)
    colnames(tdf)[ncol(tdf)] <- family
  }
  return(list(pdf,tdf))
}


##############
get_effects <- function(){
  
  spec_vec <- c()
  fam_vec <- c()
  motif_vec <- c()
  val_vec <- c()
  
  for (organism in c("Bnapus","Athaliana")){
    for (type in c("TFBS", "reversed", "scrambled", "random")){
      data <- get_mutation_data(organism, type=type)
      pdf <- data[[1]]
      tdf <- data[[2]]
      axis=F
      
      for (f in families){
        R=1
        N=101
        P=1
        pdf[,f] <- sgolayfilt(pdf[,f], p = P, n = N)
        tdf[,f] <- sgolayfilt(tdf[,f], p = P, n = N)
        # keep track of average effects
        spec_vec <- c(spec_vec, organism)
        fam_vec <- c(fam_vec, f)
        motif_vec <- c(motif_vec, type)
        val_vec <- c(val_vec, round(max(c(pdf[,f], tdf[,f])),3))
      }
    }
  }
  effect_df <- data.frame("species"=spec_vec, "family"=fam_vec, "motif_type"=motif_vec, "value"=val_vec)
  effect_matrix <- matrix(ncol=4, nrow=15)
  colnames(effect_matrix) <- c("TFBS", "Reversed", "Scrambled", "Random")
  rownames(effect_matrix) <- families_full
  c <- 0
  for (type in c("TFBS", "reversed", "scrambled", "random")){
    c <- c+1
    r <- 0
    for (f in families){
      r <- r+1
      Bn_ref <- effect_df$value[effect_df$species=="Bnapus" & effect_df$family==f & effect_df$motif_type=="TFBS"]
      At_ref <- effect_df$value[effect_df$species=="Athaliana" & effect_df$family==f & effect_df$motif_type=="TFBS"]
      Bn_val <- round(100*((effect_df$value[effect_df$species=="Bnapus" & effect_df$family==f & effect_df$motif_type==type])/Bn_ref),0)
      At_val <- round(100*((effect_df$value[effect_df$species=="Athaliana" & effect_df$family==f & effect_df$motif_type==type])/At_ref),0)
      effect_matrix[r,c] <- paste0(as.character(Bn_val),"% (", At_val, "%)")
    }
  }
  return(effect_matrix)
}

relative_positive_effect_sizes <- get_effects()
#############


plot_subfig <- function( type, labels=c( expression(bold("A")), expression(bold("B")) ), descriptor="Reversed Motifs", min_max=NULL, add_distr=F){
  for (organism in c("Bnapus","Athaliana")){   
    data <- get_mutation_data(organism, type=type)
    pdf <- data[[1]]
    tdf <- data[[2]]
    axis=F
    
    if (organism=="Athaliana"){
      label=labels[2]
      if (type=="random"){
        axis=T
      }
    } else {
      label=labels[1]
      axis=F
    }
    
    for (f in families){
      R=1
      N=101
      P=1
      pdf[,f] <- sgolayfilt(pdf[,f], p = P, n = N)
      tdf[,f] <- sgolayfilt(tdf[,f], p = P, n = N)
      # keep track of average effects
      spec_vec <- c(spec_vec, organism)
      fam_vec <- c(fam_vec, f)
      motif_vec <- c(motif_vec, type)
      val_vec <- c(val_vec, mean(abs(pdf[,f])))
    }
    
    total_max <- round(max(as.matrix(rbind(pdf[,families],tdf[,families]))),3)
    total_min <- round(min(as.matrix(rbind(pdf[,families],tdf[,families]))),3)
    
    # scale with respect to overall absolute maximum
    norm.cat <- min_max.zero_center.scale(as.matrix(rbind(pdf[,families], tdf[,families])), min_max=min_max)
    pdf[,families] <- norm.cat[1:(nrow(norm.cat)/2),]
    tdf[,families] <- norm.cat[((nrow(norm.cat)/2)+1):nrow(norm.cat),]
    
    legend_range <- as.character(max(abs(c(total_min,total_max))))
    legend_labs <- c(paste0("-",as.character(legend_range)), paste0("+",as.character(legend_range)))
    
    plot_info(organism, label=label)
    plot_heatmap(pdf, seq="Promoter", palette="custom", midline_col="black", border_col=transparent("grey66", alpha=0.3), axis=axis)
    plot_heatmap(tdf, seq="Terminator", palette="custom", midline_col="black", border_col=transparent("grey66", alpha=0.3), axis=axis)
    plot_legend(palette="custom", descriptor="(Mutated \u2212 Wild-Type)", legend_labs=legend_labs, middle="0", min_max=c(total_min/max(abs(overall_min_max)), total_max/max(abs(overall_min_max))))
    if (add_distr){
      norm.cat <- min_max.positive.scale(as.matrix(rbind(pdf[,families], tdf[,families])))
      x <- c()
      for (f in families){
        x <- c(norm.cat[,f])
      }
      x <- (x * 50) + 25
      dens <- density(x, from=25, to=75)
      lines(x=(((-dens$y)*100)-15), y=dens$x)
    }
    text(x=0, y=126, labels=expression(Delta*" Prediction"), cex=0.55)
    text(x=0, y=113, labels=descriptor, cex=0.55)
    
  }
}




pdf("supp_figs/fig_S24.pdf", height=20, width=17, pointsize=27)
par(xpd=T)
layout(matrix(c(0,1,2,0,
                3,4,5,6,
                7,8,9,10,
                11,12,13,14,
                15,16,17,18,
                19,20,21,22,
                23,24,25,26,
                27,28,29,30,
                31,32,33,34,
                0,35,36,0), ncol=4, byrow=T),
       heights=c(0.2, rep(1,8), 0.5),
       widths=c(1.5,2.5,2.5,0.8))

seqs <- c("promoter", "terminator")

lab_cex=0.75
ax_cex=0.6
main_cex=1

plot_header("Promoter")
plot_header("Terminator")

plot_subfig( type="TFBS", labels=c( expression(bold("A")), expression(bold("B")) ), descriptor=expression(bold("TFBS Motifs")), min_max=overall_min_max, add_distr=T)
plot_subfig( type="reversed", labels=c( expression(bold("C")), expression(bold("D")) ), descriptor=expression(bold("Reversed Motifs")), min_max=overall_min_max, add_distr=T)
plot_subfig( type="scrambled", labels=c( expression(bold("E")), expression(bold("F")) ), descriptor=expression(bold("Scrambled Motifs")), min_max=overall_min_max, add_distr=T)
plot_subfig( type="random", labels=c( expression(bold("G")), expression(bold("H")) ), descriptor=expression(bold("Random Motifs")), min_max=overall_min_max, add_distr=T)
plot_axlab("Promoter")
plot_axlab("Terminator")

dev.off()



######################################
#######################################


pdf("fig_6.pdf", height=12, width=17, pointsize=30)
par(xpd=T)
layout(matrix(c(0,1,2,0,
                3,4,5,6,
                7,8,9,10,
                11,12,13,14,
                15,16,17,18,
                0,19,20,0), ncol=4, byrow=T),
       heights=c(0.2, rep(1,4), 0.5),
       widths=c(1.5,2.5,2.5,0.8))


seqs <- c("promoter", "terminator")

lab_cex=0.75
ax_cex=0.6
main_cex=1

plot_header("Promoter")
plot_header("Terminator")


plot_subfig <- function( type, labels=c( expression(bold("A")), expression(bold("B")) ), descriptor="Reversed Motifs" ){
  for (organism in c("Bnapus","Athaliana")){   
    data <- get_mutation_data(organism, control=type)
    pdf <- data[[1]]
    tdf <- data[[2]]
    axis=F
    
    if (organism=="Athaliana"){
      label=labels[2]
      if (type=="random"){
        axis=T
      }
    } else {
      label=labels[1]
      axis=F
    }
    
    for (f in families){
      R=1
      N=101
      P=1
      pdf[,f] <- sgolayfilt(pdf[,f], p = P, n = N)
      tdf[,f] <- sgolayfilt(tdf[,f], p = P, n = N)
    }
    total_max <- round(max(as.matrix(rbind(pdf[,families],tdf[,families]))),3)
    total_min <- round(min(as.matrix(rbind(pdf[,families],tdf[,families]))),3)
    
    print(type)
    print("TCP max")
    print(round(max(as.matrix(rbind(pdf[,"TCP"],tdf[,"TCP"]))),3))
    print("AP2 max")
    print(round(max(as.matrix(rbind(pdf[,"AP2EREBP"],tdf[,"AP2EREBP"]))),3))
    
    # scale with respect to overall absolute maximum
    norm.cat <- min_max.zero_center.scale(as.matrix(rbind(pdf[,families], tdf[,families])))
    pdf[,families] <- norm.cat[1:(nrow(norm.cat)/2),]
    tdf[,families] <- norm.cat[((nrow(norm.cat)/2)+1):nrow(norm.cat),]
    
    legend_range <- as.character(max(abs(c(total_min,total_max))))
    legend_labs <- c(paste0("-",as.character(legend_range)), paste0("+",as.character(legend_range)))
    if (type == "random"){
      contr_lab <- "Random"
    } else {
      contr_lab <- "Scrambled"
    }  
    plot_info(organism, label=label)
    plot_heatmap(pdf, seq="Promoter", palette="custom", midline_col="black", border_col=transparent("grey66", alpha=0.3), axis=axis)
    plot_heatmap(tdf, seq="Terminator", palette="custom", midline_col="black", border_col=transparent("grey66", alpha=0.3), axis=axis)
    plot_legend(palette="custom", descriptor=paste0("(TFBS \u2212 ",contr_lab,")"), legend_labs=legend_labs, middle="0")
    text(x=0, y=113, labels=expression(Delta*" Prediction"), cex=0.55)
    #text(x=0, y=113, labels=descriptor, cex=0.55)
    
  }
}

plot_subfig( type="scrambled", labels=c( expression(bold("A")), expression(bold("B")) ), descriptor=expression(bold("Scrambled Motifs")) )
plot_subfig( type="random", labels=c( expression(bold("C")), expression(bold("D")) ), descriptor=expression(bold("Random Motifs")) )
plot_axlab("Promoter")
plot_axlab("Terminator")

dev.off()


#############################################################################
#############################################################################


plot_boxes <- function(organism="Bnapus"){
  par(mar=c(3,2,0,0), oma=c(0,5,1,0), xpd=T)
  colors <- c("palegoldenrod", "yellow", "lawngreen", "turquoise")
  colors <- rev(colors)
  y_range=c(-0.1,0.1)
  for (start in c(-500,-250,0,250)){
    stop <- (start+249)
    interval=c(start:stop)
    offset=0
    idx=1
    for (type in c("TFBS", "reversed", "scrambled", "random")){
      d <- matrix(ncol=length(families), nrow=40000)
      colnames(d) <- families
      for (f in families){
        df <- read.table(paste0("../result_subset/nemo/",organism,"/masked_graphpart/motif_insertion/", f, "_medium.tsv"), header=T)
        df <- df[(df$motif_type==type & df$mutated_sequence_type=="promoter"),]
        df$diff <- (df$mutated_pred - df$baseline)
        response <- c()
        for (p in interval){
          response <- c(response, median(df$diff[df$relative_idx==p]))
          #response <- c(response, df$diff[df$relative_idx==p])
        }
        d[,f] <- response
      }
      d <- d[,rev(colnames(d))]
      if (type == "TFBS"){
        boxplot(d, boxwex=0.15, at = (seq(1,(2*length(families)),2)+offset), xlim=c(0.9,31), col=colors[idx], outcol=colors[idx], axes=F, outline=F, horizontal = T)
        lines(y=c(1,30), x=c(0,0), lty=1, col="grey50", lwd=1.3)
        boxplot(d, boxwex=0.15, at = (seq(1,(2*length(families)),2)+offset), col=colors[idx], outcol=colors[idx], axes=F, outline=F, add=T, horizontal = T)
        text(x=0, y=31, adj=0.5, labels=bquote(bold("["*.(start)*", "*.(stop)*"]")), srt=0, col="black", cex=0.8)        } else {
          boxplot(d, boxwex=0.15, add=T, at = (seq(1,(2*length(families)),2)+offset), col=colors[idx], outcol=colors[idx], axes=F, outline=F, horizontal = T) 
        }
      offset <- offset+0.2
      idx <- idx+1
    }
    if (start == -500){
      axis(2, at=(seq(1,(2*length(families)),2)+0.3), labels=rev(families_full), las=2)
    }
    axis(1)
  }
}

############################


pdf("supp_figs/fig_S25.pdf", height=20, width=17, pointsize=30)
layout(matrix(c(1,2,3,4,6,
                5,5,5,5,0), ncol=5, byrow=T),
       heights=c(8,0.2),
       widths=c(2,2,2,2,1))
plot_boxes("Bnapus")
par(mar=c(0,0,0,0))

plot(x=c(1:100), y=(1:100), type="n", axes=F, ylab="", xlab="")
text(x=50,  y=50, labels=expression(Delta*" Prediction (Mutated - Wild-Type)"), cex=1, srt=0, adj=0.5)

plot(x=c(1:100), y=(1:100), type="n", axes=F, ylab="", xlab="")
colors <- rev(c("palegoldenrod", "yellow", "lawngreen", "turquoise"))
points(x=rep(50,4), y=c(30,40,50,60), col=colors, cex=2, pch=15)
text(x=rep(50,4),  y=c(33,43,53,63), labels=c("TFBS", "Reversed\nTFBS", "Scrambled\nTFBS", "Random"), cex=0.9)
text(x=50,  y=70, adj=0.5, labels=expression(bold(atop("Inserted","Motif"))), cex=0.9)

dev.off()


#################################

pdf("supp_figs/fig_S26.pdf", height=20, width=17, pointsize=30)
layout(matrix(c(1,2,3,4,6,
                5,5,5,5,0), ncol=5, byrow=T),
       heights=c(8,0.2),
       widths=c(2,2,2,2,1))
plot_boxes("Athaliana")
par(mar=c(0,0,0,0))

plot(x=c(1:100), y=(1:100), type="n", axes=F, ylab="", xlab="")
text(x=50,  y=50, labels=expression(Delta*" Prediction (Mutated - Wild-Type)"), cex=1, srt=0, adj=0.5)

plot(x=c(1:100), y=(1:100), type="n", axes=F, ylab="", xlab="")
colors <- rev(c("palegoldenrod", "yellow", "lawngreen", "turquoise"))
points(x=rep(50,4), y=c(30,40,50,60), col=colors, cex=2, pch=15)
text(x=rep(50,4),  y=c(33,43,53,63), labels=c("TFBS", "Reversed\nTFBS", "Scrambled\nTFBS", "Random"), cex=0.9)
text(x=50,  y=70, adj=0.5, labels=expression(bold(atop("Inserted","Motif"))), cex=0.9)

dev.off()
#########################