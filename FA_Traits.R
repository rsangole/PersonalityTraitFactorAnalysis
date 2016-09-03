library('ggplot2')
library('dplyr')
library('tidyr')
library('corrplot')
library('nFactors')
library('NbClust')

setwd("~/Documents/NWU MSPA/410/Data_All_160811")

plot_factors <- function(df,f1,f2,nfact){
    factorloadings <- as.data.frame(df$loadings[,c(f1,f2)])
    g <- ggplot(data = factorloadings, aes(x = factorloadings[[1]], y = factorloadings[[2]]))+
        geom_label(aes(label=rownames(factorloadings)),cex=3.5,alpha=.5,fill='lightcyan')+
        theme_light()+
        labs(title=paste0('Factor ',f1,' vs Factor ',f2),
             x=paste0('Factor ',f1),
             y=paste0('Factor ',f2))+
        scale_x_continuous(breaks = seq(-1,1,.5),limits = c(-1,1))+
        scale_y_continuous(breaks = seq(-1,1,.5),limits = c(-1,1))
    ggsave(plot = g,filename = paste0('f',f1,'vsf',f2,'nf',nfact,'.png'),width = 150,height = 150,units = 'mm')
    return(g)
}

# Read the data
df <- tbl_df(read.csv('CSV/Sheet_1.csv'))
df <- df %>% dplyr::select(-1:-4)

# Correlation Plot
corrplot(corr = cor(df),
         method = 'square',
         type='full',
         order = 'hclust',
         hclust.method = 'ward.D2',
         tl.col = 'gray10',
         addrect = 6)

#How many factors are needed? - Method 1
plotuScree(x = cor(df))
ev <- as.data.frame(eigen(cor(df))$values) # get eigenvalues
colnames(ev) <- 'EigenVal'
ev <- ev %>%
    mutate(Prop=EigenVal/sum(EigenVal),
           CuSum=cumsum(Prop))
print(ev)

#Method 2
ev <- eigen(cor(df)) # get eigenvalues
ap <- parallel(subject=nrow(df),var=ncol(df),rep=100,cent=.05)
nS <- nScree(x=ev$values, cor = T,model = 'factors',aparallel=ap$eigen$qevpea)
plotnScree(nS)

nfact=9

# df.fa.unrotated  <-  factanal(df, factors = nfact, rotation = 'none')
df.fa.rotated <-  factanal(df, factors = nfact, rotation = 'varimax')
print(df.fa.rotated,digits = 3,cutoff=.2,sort=T)

plot_factors(df.fa.rotated,1,2,nfact)
plot_factors(df.fa.rotated,2,3,nfact)
plot_factors(df.fa.rotated,3,4,nfact)
plot_factors(df.fa.rotated,4,5,nfact)