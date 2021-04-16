library(tidyverse)
library(ggplot2)

permute_axis <- function(dat,
                         var = 'resi_tumor',
                         plot.var = 'resi_tumor',
                         B=1000,
                         permute=T){
    
    ## canonical correlation
    cancor = cancor(dat %>% select(PC1,PC2),dat[,var])
    obs.cor = cancor$cor
    
    if(permute){
        
        new.idx = t(sapply(1:B,function(b){
            set.seed(b)
            vec = sample(1:nrow(dat),replace = F)
            if(identical(vec,1:nrow(dat))) vec = NULL
            
            return(vec)
            
        })) %>% unique
        print(paste0('#unique permuted rows: ',nrow(new.idx)))
        
        new.cor = sapply(1:nrow(new.idx),function(i){
                              
                              newdat = dat
                              newdat[,var] = dat[new.idx[i,],var]
                              cancor(newdat %>% select(PC1,PC2),newdat[,var])$cor
                              
                          })
        
        p.value = sum(new.cor>=obs.cor)/nrow(new.idx)
        print(paste0('=== Observed cancor: ',obs.cor))
        print(paste0('=== P-value: ',p.value))
        
        list(cancor = cancor, obs.cor = obs.cor,permute.cor = new.cor,p.value = p.value)
  
    }
}

##
file = 'dat_hvg_m2.rds'
dat = readRDS(paste0(output.dir,file)) # result file by running 2-xpca_xx.R (dat_hvg_m2.rds)
dat$response = ifelse(dat$response=='NR',0,1)
dat$tissue = ifelse(dat$tissue=='normal',0,1)

# Example: response
res = permute_axis(dat, var = 'response', plot.var = 'response', B = 10000, permute = T)
# Example: tissue
res = permute_axis(dat, var = 'tissue', plot.var = 'tissue', B = 10000, permute = T)
