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

# function that takes vector of values from 0 to one and translates them into gray scale colors
grey_scale <- function(x){
  idx <- 1+(((1-x)**(3/2))*999)
  cols <- grey.colors(n=1000, start=0.0, end=0.95)[idx]
  return(cols)
}
 
png("fig_5.png", width=17, height=15, units="cm", res=1200)

lab_cex=0.8
ax_cex=0.7
text_cex=0.8
fig_lab_cex=1.5

layout(
  matrix(c(1,2,
           3,4), ncol=2, byrow=TRUE), 
  widths=c(3,1.9), 
  heights=c(3,3)
)

org <- "Bnapus"
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

#se_file <- "../data/plot_data/Superenhancers/cognategenes_SEs_chr_all_v3-v1.tr.tsv"
seO_file <- "../data/plot_data/Superenhancers/SE_overpredicted_v3-v1.tsv"
seW_file <- "../data/plot_data/Superenhancers/SE_wellpredicted_v3-v1.tsv"
seU_file <- "../data/plot_data/Superenhancers/SE_underpredicted_v3-v1.tsv"
me_file <- "../data/plot_data/H3K27me3/Marked.genes.Fig4.Ex.ids"

#eID <- read.table(se_file, header=F, sep="\t") 
#eID <- eID[,2]
eoID <- read.table(seO_file, header=F, sep="\t") 
eoID <- eoID[,2]
ewID <- read.table(seW_file, header=F, sep="\t") 
ewID <- ewID[,2]
euID <- read.table(seU_file, header=F, sep="\t") 
euID <- euID[,2]
eID <- c(eoID,ewID,euID)

mID <- read.table(me_file, header=F, sep="\t") 
mID <- gsub(pattern=".1", x=mID[,1], replacement="", fixed=T)
# exclude genes that have been filtered or are not included in expression set
eID <- eID[eID %in% df$ID]
mID <- mID[mID %in% df$ID]
epreds <- df[df$ID %in% eID,]
mpreds <- df[df$ID %in% mID,]


mpreds$cols <- rep(transparent("magenta",0.3),nrow(mpreds))
mpreds$pch <- rep(24, nrow(mpreds))
epreds$cols <- rep(transparent("lightseagreen",0.3),nrow(epreds))
epreds$pch <- rep(25, nrow(epreds))
em_preds <- rbind(mpreds, epreds)
# randomly shuffle plotting order
em_preds <- em_preds[sample(nrow(em_preds)),]

points(x=em_preds$Actual,y=em_preds$Predicted, bg=em_preds$cols,
       col=NULL, pch=em_preds$pch, cex=1.2, lwd=1)

abline(a=0, b=1, col="green", lwd=1.5)
abline(a=(-cutoff), b=1, col="blue", lwd=1.5)
abline(a=cutoff, b=1, col="gold", lwd=1.5)


legend("topleft", legend=c(expression("y = x" + "Median( "*abs(""~epsilon~"")*" )"),
                           "y = x",
                           expression("y = x" - "Median( "*abs(""~epsilon~"")*" )")),
       box.col=transparent("white",0.2),
       bg=transparent("white",0.2),
       box.lwd = 0,
       lty=c(1,1,1),
       col=c("gold", "green", "blue"),
       seg.len=0.9, cex=0.9, lwd=1.5)
legend("bottomright", legend=expression(italic("B. napus")), pch=NA, bty="n", text.col="grey40")


box()
fig_label(expression(bold("A")), cex=fig_lab_cex) 



par(mar=c(1,0,0,0.5))

plot(x=c(1:100, 1:41), y=c(-20:120), type="n", axes=F, xlab="", ylab="")

oID <- read.table(paste0("../results/nemo/",org,"/masked_graphpart/IDs/prediction/overpredicted.lst"), 
                  header=F, sep="\t")[,1]
wID <- read.table(paste0("../results/nemo/",org,"/masked_graphpart/IDs/prediction/well-predicted.lst"), 
                  header=F, sep="\t")[,1]
uID <- read.table(paste0("../results/nemo/",org,"/masked_graphpart/IDs/prediction/underpredicted.lst"), 
                  header=F, sep="\t")[,1]





ew <- ((sum(eID %in% wID) + sum(eID %in% uID)) / (length(eID))) * 100
mw <- ((sum(mID %in% wID) + sum(mID %in% uID)) / (length(mID))) * 100
ALLw <- ((length(wID) + length(uID)) / (length(oID) + length(wID) + length(uID))) * 100
eu <- (sum(eID %in% uID) / (length(eID))) * 100
mu <- (sum(mID %in% uID) / (length(mID))) * 100
ALLu <- (length(uID) / (length(oID) + length(wID) + length(uID))) * 100


EO <- sum(eID %in% oID)
EW <- sum(eID %in% wID)
EU <- sum(eID %in% uID)

MO <- sum(mID %in% oID)
MW <- sum(mID %in% wID)
MU <- sum(mID %in% uID)

contingency_table <- matrix(c(EO,length(oID),
         EW,length(wID),
         EU,length(uID)),
         byrow=T, ncol=2)
NE <- sum(contingency_table[,1])
print(NE)
Nall <- sum(contingency_table[,2])
print(Nall)

colnames(contingency_table) <- c("Enhancer", "All")
rownames(contingency_table) <- c("Overpredicted", "Well-predicted", "Underpredicted")

#
print(chisq.test(contingency_table))

contingency_table <- matrix(c(MO,length(oID),
                              MW,length(wID),
                              MU,length(uID)),
                            byrow=T, ncol=2)
NM <- sum(contingency_table[,1])
print(NM)

colnames(contingency_table) <- c("H3K27me3", "All")
rownames(contingency_table) <- c("Overpredicted", "Well-predicted", "Underpredicted")

#
print(chisq.test(contingency_table)) 


if (org == "Athaliana"){
  ew_lst <- eID[(eID %in% wID)]
  write.table(as.data.frame(ew_lst), file="enhancer.well-predicted.txt", row.names=F,col.names=F,quote=F,sep="\t")
  eo_lst <- eID[(eID %in% oID)]
  write.table(as.data.frame(eo_lst), file="enhancer.overpredicted.txt", row.names=F,col.names=F,quote=F,sep="\t")
  eu_lst <- eID[(eID %in% uID)]
  write.table(as.data.frame(eu_lst), file="enhancer.underpredicted.txt", row.names=F,col.names=F,quote=F,sep="\t")
}


# enhancer
polygon(x=c(5,5,25,25), y=c(1,100,100,1), col="lemonchiffon")
polygon(x=c(5,5,25,25), y=c(1,ew,ew,1), col="darkolivegreen1")
polygon(x=c(5,5,25,25), y=c(1,eu,eu,1), col="darkslategray1")

text(adj=0.5, x=15, y=((100+ew)/2), 
     labels=paste(as.character(round((100-ew),1)), "%", sep=""), cex=text_cex)
text(adj=0.5, x=15, y=((ew+eu)/2), 
     labels=paste(as.character(round((ew-eu),1)), "%", sep=""), cex=text_cex)
text(adj=0.5, x=15, y=((eu)/2), 
     labels=paste(as.character(round((eu),1)), "%", sep=""), cex=text_cex)

#all
polygon(x=c(40,40,60,60), y=c(1,100,100,1), col="lemonchiffon")
polygon(x=c(40,40,60,60), y=c(1,ALLw,ALLw,1), col="darkolivegreen1")
polygon(x=c(40,40,60,60), y=c(1,ALLu,ALLu,1), col="darkslategray1")

text(adj=0.5, x=50, y=((100+ALLw)/2), 
     labels=paste(as.character(round((100-ALLw),1)), "%", sep=""), cex=text_cex)
text(adj=0.5, x=50, y=((ALLw+ALLu)/2), 
     labels=paste(as.character(round((ALLw-ALLu),1)), "%", sep=""), cex=text_cex)
text(adj=0.5, x=50, y=((ALLu)/2), 
     labels=paste(as.character(round((ALLu),1)), "%", sep=""), cex=text_cex)

#methylation
polygon(x=c(75,75,95,95), y=c(1,100,100,1), col="lemonchiffon")
polygon(x=c(75,75,95,95), y=c(1,mw,mw,1), col="darkolivegreen1")
polygon(x=c(75,75,95,95), y=c(1,mu,mu,1), col="darkslategray1")

text(adj=0.5, x=85, y=((100+mw)/2), 
     labels=paste(as.character(round((100-mw),1)), "%", sep=""), cex=text_cex)
text(adj=0.5, x=85, y=((mw+mu)/2), 
     labels=paste(as.character(round((mw-mu),1)), "%", sep=""), cex=text_cex)
text(adj=0.5, x=85, y=((mu)/2), 
     labels=paste(as.character(round((mu),1)), "%", sep=""), cex=text_cex)

text(adj=0.5, x=c(15, 50, 85), y=c(110, 110, 110), 
     labels=c("SE-Associated", "All\nGenes", "H3K27me3-\nAssociated"), 
     cex=text_cex)
text(adj=0.5, x=c(32.5, 67.5), y=c(124, 124), 
     labels=c("***", "***"), 
     cex=text_cex, col="grey40")
lines(x=c(17,48), y=c(120,120), col="grey40")
lines(x=c(83,52), y=c(120,120), col="grey40")
lines(x=c(17,17), y=c(117.5,120), col="grey40")
lines(x=c(48,48), y=c(117.5,120), col="grey40")
lines(x=c(83,83), y=c(117.5,120), col="grey40")
lines(x=c(52,52), y=c(117.5,120), col="grey40")
lines(x=c(32.5,32.5), y=c(120,122), col="grey40")
lines(x=c(67.5,67.5), y=c(120,122), col="grey40")

polygon(x=c(5,5,95,95), y=c(-25,-2,-2,-25))
text(x=c(16,16,54,54,54), y=c(-8,-18,-6, -13, -20), 
     labels=c("H3K27me3", "Enhancer", 
              "Overpredicted", "Well-Predicted", "Underpredicted"),
     adj=0, cex=0.9)
points(x=12, y=-8, pch=24, col=NULL, bg="magenta", cex=1.2)
points(x=12, y=-18, pch=25, col=NULL, bg="lightseagreen", cex=1.2)
points(x=50, y=-6, pch=22, col="black", bg="lemonchiffon", cex=1.5)
points(x=50, y=-13, pch=22, col="black", bg="darkolivegreen1", cex=1.5)
points(x=50, y=-20, pch=22, col="black", bg="darkslategray1", cex=1.5)





# B

init <- strsplit(org, split="", fixed=T)[[1]][1]
df <- read.table(paste0("../results/nemo/",org,"/masked_graphpart/nemo",init,"_preds/concat_preds.tsv"), header=T, sep="\t")
df$Diff <- df$Predicted - df$Actual

sub <- df[(df$Diff>=-2 & df$Diff<=2),]
sub <- sub[order(sub$Diff, decreasing=F),]
sub$Frac <- c(1:nrow(sub))/nrow(sub)


if (org == "Bnapus"){
  se_file <- "../data/plot_data/Superenhancers/cognategenes_SEs_chr_all_v3-v1.tr.tsv"
  me_file <- "../data/plot_data/H3K27me3/Marked.genes.Fig4.Ex.ids"
} else {
  se_file <- "../data/plot_data/H3K27me3/super.chr.pos.s.cognate.ids"
  me_file <- "../data/plot_data/H3K27me3/Marked.genes.Fig4"
}
eID <- read.table(se_file, header=F, sep="\t") 

if (org == "Bnapus"){
  eID <- eID[eID[,2] %in% sub$ID,2]
} else {
  eID <- eID[,1]
} 
mID <- read.table(me_file, header=F, sep="\t") 
mID <- gsub(pattern=".1", x=mID[,1], replacement="", fixed=T)


e.sub <- sub[sub$ID %in% eID,]
e.sub$Frac <- c(1:nrow(e.sub))/nrow(e.sub)
m.sub <- sub[sub$ID %in% mID,]
m.sub$Frac <- c(1:nrow(m.sub))/nrow(m.sub)

par(mar=c(3,3,2,0.1), mgp=c(1.4,0.4,0), tck=-0.02, xpd=F)
plot(sub$Diff, sub$Frac, type="n",
     ylab = "Cumulative Fraction",
     xlab = expression("Predicted - Observed ("~epsilon~")"),
     cex.axis=ax_cex, cex.lab=lab_cex)


polygon(x=c(-3,-3,3,3), y=c(-1,2,2,-1), col="darkolivegreen1", border=NA) #lightseagreen
polygon(x=c(-3,-3,-cutoff,-cutoff), y=c(-1,2,2,-1), col="darkslategray1", border=NA)
polygon(x=c(3,3,cutoff,cutoff), y=c(-1,2,2,-1), col="lemonchiffon", border=NA) #red4


lines(sub$Diff, sub$Frac, lwd=2, col="grey50")

############## enhancer ecdf
prev.y <- 0
prev.x <- -2
for(i in 1:nrow(e.sub)){
  lines(x=c(prev.x,e.sub$Diff[i]), y=c(prev.y,prev.y), 
        col="lightseagreen", lwd=2)
  prev.x <- e.sub$Diff[i]
  lines(x=c(e.sub$Diff[i], e.sub$Diff[i]), y=c(prev.y, e.sub$Frac[i]), 
        col="lightseagreen", lwd=2)
  prev.y <- e.sub$Frac[i]
}
lines(x=c(e.sub$Diff[nrow(e.sub)],2), y=c(1,1), col="lightseagreen", lwd=2)
############### methylation ecdf
prev.y <- 0
prev.x <- -2
for(i in 1:nrow(m.sub)){
  lines(x=c(prev.x,m.sub$Diff[i]), y=c(prev.y,prev.y), 
        col="magenta", lwd=2)
  prev.x <- m.sub$Diff[i]
  lines(x=c(m.sub$Diff[i], m.sub$Diff[i]), y=c(prev.y, m.sub$Frac[i]), 
        col="magenta", lwd=2)
  prev.y <- m.sub$Frac[i]
}
lines(x=c(m.sub$Diff[nrow(m.sub)],2), y=c(1,1), col="magenta", lwd=2)

legend("topleft", lty=c(1,1,1), seg.len=0.9, 
       legend=c("H3K27me3-Associated", "All Genes", "SE-Associated"), 
       lwd=2, cex=0.9,
       col=c("magenta", "grey50", "lightseagreen"), bty="n")



legend("bottomright", legend=expression(italic("B. napus")), pch=NA, bty="n", text.col="grey40")
fig_label(expression(bold("B")), cex=fig_lab_cex) 




# C

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

legend("topright", legend=expression(italic("B. napus")), text.col="grey40", bty="n")
fig_label(expression(bold("C")), cex=fig_lab_cex) 



dev.off()  

