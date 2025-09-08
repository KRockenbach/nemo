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


organisms <- c("Athaliana", "Bnapus")


families <- c("MYBrelated", "Homeobox", "C2C2dof", "WRKY",
              "TCP", "bZIP", "MYB", "MADS", "bHLH", "Trihelix",
              "NAC", "G2like", "HSF", "AP2EREBP", "C2C2gata")




for (f in 1:length(families)){
    fam <- families[f]
    for (o in 1:2){
    
        organism <- organisms[o]
        TFs <- as.vector(read.table(paste0("../results/nemo/", organism, 
                             "/masked_graphpart/", fam, "_insertion/TF_names.lst"),
                      header=F, sep="\t")[,1])
    
        seqs <- c("promoter", "terminator")
    
        for (s in 1:2){
      
            for (t in 1:length(TFs)){
                TF <- TFs[t]
                print(TF)
                seq <- seqs[s]
                medium <- paste0("../results/nemo/", organism, "/masked_graphpart/", 
                             fam, "_insertion/", TF, "_medium_exp.tsv")
                medium_df <- read.table(medium, header=T, sep="\t")
                medium_df$prom_insert_idx <- (medium_df$prom_insert_idx - 5000)
                medium_df$term_insert_idx <- (medium_df$term_insert_idx - 1200)
            
                if (seq == "promoter"){
                    medium_resid <- (medium_df$prom_mut_preds - medium_df$baseline)
                } else {
                    medium_resid <- (medium_df$term_mut_preds - medium_df$baseline)
                }
                y <- c()
                x <- c()
                idx <- 1
        
                for (i in -1000:999){
                    x[idx] <- i
                    y[idx] <- mean(medium_resid[medium_df$prom_insert_idx==i], na.rm=T)
                    idx <- idx + 1
                }
        
                if (t == 1){
                    y_mat <- matrix(y)
                } else {
                    y_mat <- cbind(y_mat, y)
                }
                if (t == length(TFs)){
                    y <- apply(y_mat, MARGIN = 1, FUN = mean, na.rm=T)
                    outname <- paste0("../data/plot_data/", organism, "_", seq,"_TF_inserted.tsv")
                    if (file.exists(outname)){
                        df <- read.table(outname, sep="\t", header=T)
                        df[,f+1] <- y
                        colnames(df) <- c("x", families[1:f])
                    } else {
                        df <- data.frame(x,y)
                        colnames(df) <- c("x", families[1:f])
                    }
          
                    write.table(df, file=outname, col.names = T, row.names = F, 
                                sep="\t", quote = F)
                }
            
            }
      
        }
    
    }

}
