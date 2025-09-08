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


# script to plot model
source("./plot_utils.R")

fig_lab_cex=1.5

png("fig_1.png", width=17, height=18, units="cm", res=1200)

par(mar=c(0,0,0,0), oma=c(0,0,0,0))


plot(x=c(7.5:392.5), y=c(-285:100), type="n", xlab="", ylab="", axes =F)


lc <- 100
rc <- 300
top <- 100
# promoter 1.5 * 62 units long
lines(x=c(lc-((62*1.5)/2), lc+((62*1.5)/2)), y=c(95, 95), lwd=2)
#terminator 1.5 * 62 units long
lines(x=c(rc-((62*1.5)/2), rc+((62*1.5)/2)), y=c(95, 95), lwd=2)

TSS<-lc-((62*1.5)/2)+(50*1.5)
lines(x=c(TSS, TSS), y=c(95, 99), lwd=2)
lines(x=c(TSS, TSS+3.5), y=c(99, 99), lwd=2)
polygon(x=c(TSS+3.5,TSS+3.5,TSS+6.5), y=c(97.5, 100.5, 99), col="black")
text(x=c(TSS-25,TSS,TSS+25), y=c(100,105,100), labels = c("-5 kb", "TSS", "+1.2 kb"), cex=0.65, adj=0.5, font=2)

# TTS at 1.5 * 12 units
TTS <- rc-((62*1.5)/2)+(12*1.5)
lines(x=c(TTS, TTS), y=c(95, 99), lwd=2)
lines(x=c(TTS-2.5, TTS+2.5), y=c(99, 99), lwd=2)
text(x=c(TTS-25,TTS,TTS+25), y=c(100,105,100), labels = c("-1.2 kb", "TTS", "+5 kb"), cex=0.65, adj=0.5, font=2)


polygon(x=c(190, 210, 210, 190), y=c(top-30, top-30, top-347, top-347), col=transparent("darkorange4", 0.85), border=NA)
polygon(x=c(170, 230, 200), y=c(top-347, top-347, top-380), col=transparent("darkorange4", 0.85), border=NA)


n <- 15
seq <- stri_rand_strings(4, n, '[ACGT]')
text(x=c(lc-40, lc+40, rc-40, rc+40), y=c(90,90,90,90), labels=seq, cex=0.5, font=2, adj=0.5)
text(x=c(lc, rc), y=c(90,90), labels=c("...", "..."), cex=0.5, font=2, adj=0.5)


for (i in seq(lc-65, lc-5, 4)){
  lines(x=c(i,i), y=c(69,85))
}
for (i in seq(lc+5, lc+65, 4)){
  lines(x=c(i,i), y=c(69,85))
}

for (i in seq(rc-65, rc-5, 4)){
  lines(x=c(i,i), y=c(69,85))
}
for (i in seq(rc+5, rc+65, 4)){
  lines(x=c(i,i), y=c(69,85))
}

for (i in seq(69, 85, 4)){
  lines(x=c(lc-65,lc-5), y=c(i,i))
  lines(x=c(lc+5,lc+65), y=c(i,i))
  lines(x=c(rc-65,rc-5), y=c(i,i))
  lines(x=c(rc+5,rc+65), y=c(i,i))
}


text(x=c(lc, rc), y=c(77,77), labels=c("...", "..."), cex=0.5, font=2, adj=0.5)


for (s in 1:4){
  if (s == 1){
    x <- lc-65 
  } else if (s == 2) {
    x <- lc+5
  } else if (s == 3) {
    x <- rc-65
  } else {
    x <- rc+5
  }
  for (i in strsplit(seq[s], split="", fixed=T)[[1]]){
    if (i == "T"){
      y <- 69
    } else if (i == "G"){
      y <- 73
    } else if (i == "C"){
      y <- 77
    } else {
      y <- 81
    }
    polygon(x=c(x,x,x+4,x+4), y=c(y,y+4,y+4,y), col="navy")
    x <- x+4
  }  
}

draw_array <- function(horiz=lc, vert=top-35, width=62, height=28.8, color="palegreen2", border="darkgreen"){
  polygon(x=c(horiz-(width/2),horiz-(width/2),horiz+(width/2),horiz+(width/2)), 
          y=c(vert, vert-height, vert-height, vert), col=color, border=border)
}

widths <- c(62,62,41.4,41.4,41.4,41.4,3.2,3.2,0.6)
text_vert <- c(top-47, top-87, top-120, top-135, top-165,
               top-195, top-225, top-240, top-270)


#promoter branch 
polygon(x=c(lc+(widths[1]/2),lc+(widths[1]/2),180), 
        y=c(top-35, top-35-28.8, text_vert[1]), col="darkseagreen1", border=NA)
polygon(x=c(lc+(widths[2]/2),lc+(widths[2]/2),180), 
        y=c(top-82, top-82-10.6, text_vert[2]), col="darkseagreen1", border=NA)
polygon(x=c(lc+(widths[3]/2),lc+(widths[3]/2),180), 
        y=c(top-130, top-130-10.6, text_vert[3]), col="lightcyan", border=NA)
polygon(x=c(lc+(widths[4]/2),lc+(widths[4]/2),180), 
        y=c(top-163, top-163-22.7, text_vert[4]), col="darkseagreen1", border=NA)
polygon(x=c(lc+(widths[5]/2),lc+(widths[5]/2),180), 
        y=c(top-193, top-193-18.7, text_vert[5]), col="darkseagreen1", border=NA)
polygon(x=c(lc+(widths[6]/2),lc+(widths[6]/2),180),
        y=c(top-216, top-216-6.4, text_vert[6]), col="darkseagreen1", border=NA)
polygon(x=c(lc+(widths[7]/2),lc+(widths[7]/2),180),
        y=c(top-228, top-228-6.4, text_vert[7]), col="lightcyan", border=NA)
polygon(x=c(lc+(widths[8]/2),lc+(widths[8]/2),180),
        y=c(top-243, top-243-15.2, text_vert[8]), col="darkseagreen1", border=NA)
polygon(x=c(lc+(widths[9]/2),lc+(widths[9]/2),180),
        y=c(top-263, top-263-15.2, text_vert[9]), col="lightcyan", border=NA)


draw_array(width=widths[1])
draw_array(vert=top-82, height=10.6, width=widths[2])
draw_array(vert=top-130, height=10.6, width=widths[3], color="lightblue", border="navy")
draw_array(vert=top-163, height=22.7, width=widths[4])
draw_array(vert=top-193, height=18.7, width=widths[5])
draw_array(vert=top-216, height=6.4, width=widths[6])
draw_array(vert=top-228, height=6.4, width=widths[7], color="lightblue", border="navy")
draw_array(vert=top-243, height=15.2, width=widths[8])
draw_array(vert=top-263, height=15.2, width=widths[9], color="lightblue", border="navy")


text(x=(lc-55), y=(top-49.4), labels="6200✕288", col="grey60", cex=0.65)
text(x=(lc-55), y=(top-82-(10.6/2)), labels="6200✕106", col="grey60", cex=0.65)
text(x=(lc-55), y=(top-130-(10.6/2)), labels="414✕106", col="grey60", cex=0.65)
text(x=(lc-55), y=(top-163-(22.7/2)), labels="414✕227", col="grey60", cex=0.65)
text(x=(lc-55), y=(top-193-(18.7/2)), labels="414✕187", col="grey60", cex=0.65)
text(x=(lc-55), y=(top-216-(6.4/2)), labels="414✕64", col="grey60", cex=0.65)
text(x=(lc-55), y=(top-228-(6.4/2)), labels="32✕64", col="grey60", cex=0.65)
text(x=(lc-55), y=(top-243-(15.2/2)), labels="32✕152", col="grey60", cex=0.65)
text(x=(lc-55), y=(top-263-(15.2/2)), labels="6✕152", col="grey60", cex=0.65)


#terminator branch 
polygon(x=c(rc-(widths[1]/2),rc-(widths[1]/2),220), 
        y=c(top-35, top-35-21.9, text_vert[1]), col="darkseagreen1", border=NA)
polygon(x=c(rc-(widths[2]/2),rc-(widths[2]/2),220), 
        y=c(top-63, top-63-44.9, text_vert[2]), col="darkseagreen1", border=NA)
polygon(x=c(rc-(widths[3]/2),rc-(widths[3]/2),220), 
        y=c(top-113, top-113-44.9, text_vert[3]), col="lightcyan", border=NA)
polygon(x=c(rc-(widths[4]/2),rc-(widths[4]/2),220), 
        y=c(top-163, top-163-25.9, text_vert[4]), col="darkseagreen1", border=NA)
polygon(x=c(rc-(widths[5]/2),rc-(widths[5]/2),220), 
        y=c(top-193, top-193-15.4, text_vert[5]), col="darkseagreen1", border=NA)
polygon(x=c(rc-(widths[6]/2),rc-(widths[6]/2),220),
        y=c(top-213, top-213-10.5, text_vert[6]), col="darkseagreen1", border=NA)
polygon(x=c(rc-(widths[7]/2),rc-(widths[7]/2),220),
        y=c(top-228, top-228-10.5, text_vert[7]), col="lightcyan", border=NA)
polygon(x=c(rc-(widths[8]/2),rc-(widths[8]/2),220),
        y=c(top-243, top-243-15.2, text_vert[8]), col="darkseagreen1", border=NA)
polygon(x=c(rc-(widths[9]/2),rc-(widths[9]/2),220),
        y=c(top-263, top-263-15.2, text_vert[9]), col="lightcyan", border=NA)



draw_array(horiz=rc, height=21.9, width=widths[1])
draw_array(horiz=rc, vert=top-63, height=44.9, width=widths[2])
draw_array(horiz=rc, vert=top-113, height=44.9, width=widths[3], color="lightblue", border="navy")
draw_array(horiz=rc, vert=top-163, height=25.9, width=widths[4])
draw_array(horiz=rc, vert=top-193, height=15.4, width=widths[5])
draw_array(horiz=rc, vert=top-213, height=10.5, width=widths[6])
draw_array(horiz=rc, vert=top-228, height=10.5, width=widths[7], color="lightblue", border="navy")
draw_array(horiz=rc, vert=top-243, height=15.2, width=widths[8])
draw_array(horiz=rc, vert=top-263, height=15.2, width=widths[9], color="lightblue", border="navy")

text(x=(rc+55), y=(top-35-(21.9/2)), labels="6200✕219", col="grey60", cex=0.65)
text(x=(rc+55), y=(top-63-(44.9/2)), labels="6200✕449", col="grey60", cex=0.65)
text(x=(rc+55), y=(top-113-(44.9/2)), labels="414✕449", col="grey60", cex=0.65)
text(x=(rc+55), y=(top-163-(25.9/2)), labels="414✕259", col="grey60", cex=0.65)
text(x=(rc+55), y=(top-193-(15.4/2)), labels="414✕154", col="grey60", cex=0.65)
text(x=(rc+55), y=(top-213-(10.5/2)), labels="414✕105", col="grey60", cex=0.65)
text(x=(rc+55), y=(top-228-(10.5/2)), labels="32✕105", col="grey60", cex=0.65)
text(x=(rc+55), y=(top-243-(15.2/2)), labels="32✕152", col="grey60", cex=0.65)
text(x=(rc+55), y=(top-263-(15.2/2)), labels="6✕152", col="grey60", cex=0.65)




#conv block1
text(x=c(150,200, 250), y=text_vert[1], adj=0.5, 
     labels=c("(5, 288)","1D-Convolution","(4, 219)"), cex=0.65, col="darkgreen")
text(x=200, y=text_vert[1]-10, adj=0.5, 
     labels="Batch Normalization (0.8167)", cex=0.65, col="grey20")
text(x=200, y=text_vert[1]-20, adj=0.5, 
     labels="GELU Activation", cex=0.65, col="grey20")

text(x=c(150,200, 250), y=text_vert[2], adj=0.5, labels=c("(8, 106)","1D-Convolution","(4, 449)"), cex=0.65, col="darkgreen")
text(x=200, y=text_vert[2]-10, adj=0.5, labels="Batch Normalization (0.8167)", cex=0.65, col="grey20")
text(x=200, y=text_vert[2]-20, adj=0.5, labels="GELU Activation", cex=0.65, col="grey20")


text(x=200, y=text_vert[3], adj=0.5, labels="Average Pooling (Size = 21, Stride = 15)", cex=0.65, col="navy")

text(x=c(150,200, 250), y=text_vert[4], adj=0.5, labels=c("(9, 227)", "1D-Convolution", "(16, 259)"), cex=0.65, col="darkgreen")
text(x=200, y=text_vert[4]-10, adj=0.5, labels="Batch Normalization (0.8167)", cex=0.65, col="grey20")
text(x=200, y=text_vert[4]-20, adj=0.5, labels="GELU Activation", cex=0.65, col="grey20")

text(x=c(150,200, 250), y=text_vert[5], adj=0.5, labels=c("(10, 187)", "1D-Convolution", "(15, 154)"), cex=0.65, col="darkgreen")
text(x=200, y=text_vert[5]-10, adj=0.5, labels="Batch Normalization (0.8167)", cex=0.65, col="grey20")
text(x=200, y=text_vert[5]-20, adj=0.5, labels="GELU Activation", cex=0.65, col="grey20")

text(x=c(150,200, 250), y=text_vert[6], adj=0.5, labels=c("(46, 64)", "1D-Convolution", "(20, 105)"), cex=0.65, col="darkgreen")
text(x=200, y=text_vert[6]-10, adj=0.5, labels="Batch Normalization (0.8167)", cex=0.65, col="grey20")
text(x=200, y=text_vert[6]-20, adj=0.5, labels="GELU Activation", cex=0.65, col="grey20")

text(x=200, y=text_vert[7], adj=0.5, labels="Average Pooling (Size = 18, Stride = 13)", cex=0.65, col="navy")

text(x=c(150,200, 250), y=text_vert[8], adj=0.5, labels=c("(62, 152)", "1D-Convolution", "(22, 152)"), cex=0.65, col="darkgreen")
text(x=200, y=text_vert[8]-10, adj=0.5, labels="Batch Normalization (0.8167)", cex=0.65, col="grey20")
text(x=200, y=text_vert[8]-20, adj=0.5, labels="GELU Activation", cex=0.65, col="grey20")

text(x=200, y=text_vert[9], adj=0.5, labels="Average Pooling (Size = 8, Stride = 6)", cex=0.65, col="navy")


polygon(x=c(200-(182.4/2), 200-(182.4/2), 200+(182.4/2), 200+(182.4/2)), y=c(text_vert[9]-12, text_vert[9]-12.1, text_vert[9]-12.1, text_vert[9]-12))

text(x=(lc-55), y=(text_vert[9]-12), labels="1824", col="grey60", cex=0.65)

for (i in seq(200-(182.4/2)+0.05, 200+(182.4/2)-0.05, 0.1)){
  for (j in seq(200-128, 200-8, 8)){
    lines(x=c(i,j), y=c(text_vert[9]-12.1, text_vert[9]-30), lwd=1, col=transparent("black", 0.996))
  }
  for (j in seq(200+128, 200+8, -8)){
    lines(x=c(i,j), y=c(text_vert[9]-12.1, text_vert[9]-30), lwd=1, col=transparent("black", 0.996))
  }
}

for (i in seq(200-128, 200-8, 8)){
  for (j in c(200-8, 200, 200+8)){
    lines(x=c(i,j), y=c(text_vert[9]-30, text_vert[9]-60), lwd=1, col=transparent("grey30", 0.5))
  }
}
for (i in seq(200+128, 200+8, -8)){
  for (j in c(200-8, 200, 200+8)){
    lines(x=c(i,j), y=c(text_vert[9]-30, text_vert[9]-60), lwd=1, col=transparent("grey30", 0.5))
  }
}

for (i in c(200-8, 200, 200+8)){
  lines(x=c(i,200), y=c(text_vert[9]-60, text_vert[9]-100), lwd=1, col="grey30")
}



text(x=200, y=text_vert[9]-16.5, labels="Flatten and Concatenate", col="ivory", cex=0.65)


points(x=seq(200-128, 200-8, 8), y=rep(text_vert[9]-30, 16), cex=1, pch=21, col="black", bg="aquamarine")
text(x=200, y=(text_vert[9]-30), adj=0.5, labels="...")
points(x=seq(200+128, 200+8, -8), y=rep(text_vert[9]-30, 16), cex=1, pch=21, col="black", bg="aquamarine")

text(x=(lc-55), y=(text_vert[9]-30), labels="750", col="grey60", cex=0.65)

text(x=rc+51.5, y=text_vert[9]-39, labels="Dense Layer (750)", col="aquamarine4", cex=0.65)
text(x=rc+51.5, y=text_vert[9]-46, labels="Batch Normalization (0.8167)", col="grey20", cex=0.65)
text(x=rc+51.5, y=text_vert[9]-53, labels="GELU Activation", col="grey20", cex=0.65)



points(x=c(200-8, 200, 200+8), y=rep(text_vert[9]-60, 3), cex=1, pch=21, col="black", bg="aquamarine")

text(x=(lc-55), y=(text_vert[9]-60), labels="3", col="grey60", cex=0.65)


text(x=rc+51.5, y=text_vert[9]-67, labels="Dense Layer (3)", col="aquamarine4", cex=0.65)
text(x=rc+51.5, y=text_vert[9]-74, labels="Batch Normalization (0.8167)", col="grey20", cex=0.65)
text(x=rc+51.5, y=text_vert[9]-81, labels="GELU Activation", col="grey20", cex=0.65)
text(x=rc+51.5, y=text_vert[9]-88, labels=expression("Dropout (1.243 ⋅ 10"^-3*")"), col="grey20", cex=0.65)


points(x=200, y=(text_vert[9]-100), cex=1, pch=21, col="black", bg="aquamarine")
text(x=(lc-55), y=(text_vert[9]-100), labels="1", col="grey60", cex=0.65)


text(x=200, y=text_vert[9]-114, labels=expression("Median Expression [log"[10]*"(TPM + 0.1)]"), col="black", cex=0.8)



dev.off()
