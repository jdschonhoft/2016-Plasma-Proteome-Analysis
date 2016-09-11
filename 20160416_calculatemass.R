library(UniProt.ws)
library(seqinr)
library(xlsx)
library(mygene)



up <- UniProt.ws(taxId=9606) #loading homosapiens
columns(up) #list of uniprot 


#pulling from uniprot by entrez id.
PP_DB <- read.xlsx("20160417_plasma_MRM_orderedbyconc.xlsx", 1)
keys <- PP_DB$Gene.ID
column <- c("UNIGENE", "KEGG", "REVIEWED", "GENES", "SEQUENCE")
kt <- "ENTREZ_GENE"
res <- select(up, keys, column, kt)
res$GENES <- lapply(strsplit(as.character(res$GENES), " "), "[", 1) #splits column to first entry for gene name
res <- res[ which(res$REVIEWED=='reviewed'), ]  #only reviewed

seq_vector <- lapply(res$SEQUENCE, s2c)  #transforming sequences to character vector
res$MASS <- lapply(seq_vector, pmw) #calculating mass using seqinR
res$SEQUENCE <- NULL #removing column

#plotting
df <- res[c('GENES', 'MASS')]
df <- unique(df)
df['CONC'] = PP_DB$Protein.concentration.ug.ml
df <- df[1:45,]
df['MASS'] <- sapply(df$MASS, function(x) {x/1000.0})

jpeg("20160418_CvMW.jpeg", width = 8, height = 5, units = 'in', res = 1200)
plot(df$MASS, df$CONC, log = 'xy', main="Most Abundant Plasma Proteins", xlab="MW (kDa)", ylab="Concentration (ug/mL)", pch =19, col = 28)
text(df$MASS,df$CONC,labels=df$GENES, pos = 4, cex=0.5)
dev.off()


#useful line of code to map gene ids etc
queryMany(PP_DB$Gene.ID, scopes="symbol", fields=c("uniprot", 'entrezgene'), species="human")


