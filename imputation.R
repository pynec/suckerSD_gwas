#imputation script for running on a high performance computing cluster 
library(rhdf5)
library(stringi)

args <- commandArgs(trailingOnly = TRUE)

gprob1 = h5read(args[1], "gprob")
gprob2 = h5read(args[2], "gprob")
gprob3 = h5read(args[3], "gprob")

##following adjusted from hdf52rdata.R from beurklelab bitbucket
ploidy <- dim(gprob1)[1] - 1

gp_mean1 <- apply(gprob1, c(2,3), function(x){sum(x*(0:ploidy))})
gp_mean2 <- apply(gprob2, c(2,3), function(x){sum(x*(0:ploidy))})
gp_mean3 <- apply(gprob3, c(2,3), function(x){sum(x*(0:ploidy))})

gp.mean <- (gp_mean1+gp_mean2+gp_mean3)/3

nind = 89
nloci = 68186

##testing with all loci (68186)

##data frame to add genotype likelihoods to (GL for each loci and ind)

test_df = data.frame(matrix(ncol = nind, nrow = nloci))

for(i in 1:nloci){
      test_df[i, ] = gp.mean[, i]
 }
 
colnames(test_df) = NULL

variants = read.csv(args[4], nrows = nloci, skip = 1, header = F)

chr = variants$V1
pos = variants$V2 

##output to add to genotype file (with chr information and commas - for GEMMA)
output =  data.frame(matrix(ncol = nind + 1, nrow = nloci))

for(i in 1:nloci){
      output[i, 1] = paste0(chr[i], ", ", "x", ", ", "x")
      index = 0
      for(j in 2:(nind+2)){
      	    if(j == nind +2){
	    output[i,j] = paste0(test_df[i,index])
	    }
	    else{
	    output[i,j] = paste0(test_df[i,index], ", ")
	    index = index + 1
	    } } }
colnames(output) = NULL

##setting up genotype file
for(i in 1:nloci){
      write.table(output[i,], "geno_real_68186.txt", append = TRUE, row.names = FALSE, quote = FALSE)
}


##output for annotation file
anno = data.frame(matrix(ncol = 1, nrow = nloci))
for(i in 1:nloci){
anno[i,] = paste0(chr[i], sample(1:20,1), " ", pos[i], " ", chr[i])
}
colnames(anno) = NULL

##setting up annotation file
for(i in 1:nloci){
      write.table(anno[i,], "anno_real_68186.txt", append = TRUE, row.names = FALSE, quote = FALSE)
      }

##getting sex information on the individuals for phenotype file
details = read.csv(args[5], header = T)
sex = details$sex
sex = gsub("Male", "0", sex)
sex = gsub("Female", "1", sex)


