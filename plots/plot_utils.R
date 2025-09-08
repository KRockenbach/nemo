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


package_dir <- "~/R/nemo_plotting_packages/"
dependencies <- c("data.table", "pracma", "stringr", "viridis", "yarrr", "stringi", 
                  "MASS", "gsignal", "wesanderson", "scico", "hash")
for (lib in dependencies){ 
    if(!require(lib, character.only = TRUE)){
        dir.create(package_dir, recursive = TRUE, showWarnings = FALSE)  # create personal library
        .libPaths(package_dir)  # add to the path
        install.packages(lib, lib=package_dir)
    }
    library(lib, lib.loc=package_dir, character.only = TRUE)
}

get_density <- function(x, y, ...) {
  density <- kde2d(x, y, ...)
  i_x <- findInterval(x, density$x)
  i_y <- findInterval(y, density$y)
  i_xy <- cbind(i_x, i_y)
  return(density$z[i_xy])
}

#####################################################
#####################################################

get_tau <- function(x){
  x_hat <- x/max(x)
  tau <- sum(1-x_hat)/(length(x)-1)
  return(tau)
}

#####################################################
#####################################################

min_max.zero_center.scale <- function(x){
  to_return <- x/max(abs(x)) 
  return(to_return)
}

min_max.positive.scale <- function(x){
  to_return <- (x-min(x))/(max(x)-min(x)) 
  return(to_return)
}


# make importance sum to zero along sequence
mean_norm <- function(x){
  x <- x-mean(x)
  return(x)
}

# scale importances locally between groups
# amplify small differences in regions with low standard deviation
standard_scaling <- function(x){
  x <- (x-mean(x))/sd(x)
  return(x)
}



get_colors <- function(x, palette="viridis"){
  ## values are assumed to be between -1 and 1 for custom palette and 0 to 1 for any other
  if (palette == "custom"){
    colors <- c()
    for (i in 1:length(x)){
      if (x[i] >=0){
        r <- 1
        g <- 1 - abs(x[i])^(1.1)#(x[i])^2
        b <- 1 - abs(x[i])^(0.6) #0.88
      } else {
        r <- 1 - abs(x[i])^(1.7)
        g <- 1 - abs(x[i])^(0.88)#(x[i])^2
        b <- 1
      }
      colors[i] <- rgb(red=r, green=g, blue=b)
    }
  } else {  
    colors <- viridis(1000, option = palette)[1 + round((999 * x),0)]
    #colors <- scico(1000, palette = palette)[1 + round((999 * x),0)]
  }
  return(colors)
}


# function that takes vector of values from 0 to one and translates them into gray scale colors
grey_scale <- function(x){
  idx <- 1+(((1-x)**(3/2))*999)
  cols <- grey.colors(n=1000, start=0.0, end=0.95)[idx]
  return(cols)
}

#####################################################
#####################################################

fig_label <- function(text, region="figure", pos="topleft", cex=NULL, ...) {
  # this useful function was originally posted here:
  # https://logfc.wordpress.com/2017/03/15/adding-figure-labels-a-b-c-in-the-top-left-corner-of-the-plotting-region/
  
  region <- match.arg(region, c("figure", "plot", "device"))
  pos <- match.arg(pos, c("topleft", "top", "topright", 
                          "left", "center", "right", 
                          "bottomleft", "bottom", "bottomright"))
  
  if(region %in% c("figure", "device")) {
    ds <- dev.size("in")
    # xy coordinates of device corners in user coordinates
    x <- grconvertX(c(0, ds[1]), from="in", to="user")
    y <- grconvertY(c(0, ds[2]), from="in", to="user")
    
    # fragment of the device we use to plot
    if(region == "figure") {
      # account for the fragment of the device that 
      # the figure is using
      fig <- par("fig")
      dx <- (x[2] - x[1])
      dy <- (y[2] - y[1])
      x <- x[1] + dx * fig[1:2]
      y <- y[1] + dy * fig[3:4]
    } 
  }
  
  # much simpler if in plotting region
  if(region == "plot") {
    u <- par("usr")
    x <- u[1:2]
    y <- u[3:4]
  }
  
  sw <- strwidth(text, cex=cex) * 60/100
  sh <- strheight(text, cex=cex) * 60/100
  
  x1 <- switch(pos,
               topleft     =x[1] + sw, 
               left        =x[1] + sw,
               bottomleft  =x[1] + sw,
               top         =(x[1] + x[2])/2,
               center      =(x[1] + x[2])/2,
               bottom      =(x[1] + x[2])/2,
               topright    =x[2] - sw,
               right       =x[2] - sw,
               bottomright =x[2] - sw)
  
  y1 <- switch(pos,
               topleft     =y[2] - sh,
               top         =y[2] - sh,
               topright    =y[2] - sh,
               left        =(y[1] + y[2])/2,
               center      =(y[1] + y[2])/2,
               right       =(y[1] + y[2])/2,
               bottomleft  =y[1] + sh,
               bottom      =y[1] + sh,
               bottomright =y[1] + sh)
  
  old.par <- par(xpd=NA)
  on.exit(par(old.par))
  
  text(x1, y1, text, cex=cex, ...)
  return(invisible(c(x,y)))
}





