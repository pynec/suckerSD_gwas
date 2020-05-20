setwd("/Users/cassandrepyne/Documents")

##script that attempts to get genotype likelihoods from hdf5 files 
library(rhdf5)
library(stringi)

#h5ls("mm25maf05_k3_100k_rep1.hdf5")
##load in genotype likelihood files 
gprob1 = h5read("mm25maf05_k3_100k_rep1.hdf5", "gprob")  
gprob2 = h5read("mm25maf05_k3_100k_rep2.hdf5", "gprob")  
gprob3 = h5read("mm25maf05_k3_100k_rep3.hdf5", "gprob")  

##following from hdf52rdata.R from buerklelab bitbucket
ploidy <- dim(gprob1)[1] - 1
gp_mean1 <- apply(gprob1, c(2,3), function(x){sum(x*(0:ploidy))})
gp_mean2 <- apply(gprob2, c(2,3), function(x){sum(x*(0:ploidy))})
gp_mean3 <- apply(gprob3, c(2,3), function(x){sum(x*(0:ploidy))})

gp.mean <- (gp_mean1+gp_mean2+gp_mean3)/3

nind = 89
nloci = 68186

---------------------------------------------------------------------------
##Testing with loci = 1000

##data frame to add genotype likelihood to (GL for each loci & ind)
test_df <- data.frame(matrix(ncol = nind, nrow = 1000))
for(i in 1:1000){
	test_df[i, ] = gp.mean[, i]
}
colnames(test_df) = NULL

variants = read.csv("variants.txt", nrows = 1000, skip = 1, header = F)
chr = variants$V1
pos = variants$V2

##output to add to genotype file (with chromosome information and commas - for GEMMA)
output <- data.frame(matrix(ncol = nind + 1, nrow = 1000))
for(i in 1:1000){
	output[i, 1] = paste0(chr[i], ", ", "x", ", ", "x")
	index = 0 
	for(j in 2:(nind+2)){
		if(j == nind +2){
			output[i,j] = paste0(test_df[i,index])
		}
		else{
		output[i,j] = paste0(test_df[i,index], ", ")
		index = index + 1	
	}}	}
	colnames(output) = NULL  

##genotype file 
for(i in 1:1000){
	write.table(output[i,], "geno_real_1000.txt", append = TRUE, row.names = FALSE, quote = FALSE)
}

##output for annotation file
anno <- data.frame(matrix(ncol = 1, nrow = 1000))
for(i in 1:1000){
	anno[i, ] = paste0(chr[i], sample(1:20, 1), " ", pos[i], " ", chr[i])
}
colnames(anno) = NULL

##annotation file 
for(i in 1:1000){
	write.table(anno[i,], "anno_real_1000.txt", append = TRUE, row.names= FALSE, quote = FALSE)
}

##getting sex information on the individuals for pheontype file 
details <- read.csv("sucker_metadata_knownsex_89.csv", header=T)
sex <- details$sex 
sex <- gsub("Male", "0", sex)
sex <- gsub("Female", "1", sex)

##add to phenotype file 
##remember that it randomly adds an X to the beginning of the file 
write.table(sex, "pheno_real_1000.txt", append = TRUE, row.names = FALSE, quote = FALSE)



---------------------------------------------------------------------
##Testing with loci = 10000

##data frame to add genotype likelihood to (GL for each loci & ind)
test_df <- data.frame(matrix(ncol = nind, nrow = 10000))
for(i in 1:10000){
	test_df[i, ] = gp.mean[, i]
}
colnames(test_df) = NULL

variants = read.csv("variants.txt", nrows = 10000, skip = 1, header = F)

chr = variants$V1

pos = variants$V2

##output to add to genotype file (with chromosome information and commas - for GEMMA)
output <- data.frame(matrix(ncol = nind + 1, nrow = 10000))
for(i in 1:10000){
	output[i, 1] = paste0(chr[i], ", ", "x", ", ", "x")
	index = 0 
	for(j in 2:(nind+2)){
		if(j == nind +2){
			output[i,j] = paste0(test_df[i,index])
		}
		else{
		output[i,j] = paste0(test_df[i,index], ", ")
		index = index + 1	
	}}	}
	colnames(output) = NULL  

##genotype file 
for(i in 1:10000){
	write.table(output[i,], "geno_real_10000.txt", append = TRUE, row.names = FALSE, quote = FALSE)
}

##output for annotation file
anno <- data.frame(matrix(ncol = 1, nrow = 10000))
for(i in 1:10000){
	anno[i, ] = paste0(chr[i], sample(1:20, 1), " ", pos[i], " ", chr[i])
}
colnames(anno) = NULL

##annotation file 
for(i in 1:10000){
	write.table(anno[i,], "anno_real_10000.txt", append = TRUE, row.names= FALSE, quote = FALSE)
}

##getting sex information on the individuals for pheontype file 
details <- read.csv("sucker_metadata_knownsex_89.csv", header=T)
sex <- details$sex 
sex <- gsub("Male", "0", sex)
sex <- gsub("Female", "1", sex)

##add to phenotype file 
##remember that it randomly adds an X to the beginning of the file 
write.table(sex, "pheno_real_10000.txt", append = TRUE, row.names = FALSE, quote = FALSE)


---------------------------------------------------------------------
##Testing with loci = 20000

##data frame to add genotype likelihood to (GL for each loci & ind)
test_df <- data.frame(matrix(ncol = nind, nrow = 20000))
for(i in 1:20000){
	test_df[i, ] = gp.mean[, i]
}
colnames(test_df) = NULL

variants = read.csv("variants.txt", nrows = 20000, skip = 1, header = F)

chr = variants$V1

pos = variants$V2

##output to add to genotype file (with chromosome information and commas - for GEMMA)
output <- data.frame(matrix(ncol = nind + 1, nrow = 20000))
for(i in 1:20000){
	output[i, 1] = paste0(chr[i], ", ", "x", ", ", "x")
	index = 0 
	for(j in 2:(nind+2)){
		if(j == nind +2){
			output[i,j] = paste0(test_df[i,index])
		}
		else{
		output[i,j] = paste0(test_df[i,index], ", ")
		index = index + 1	
	}}	}
	colnames(output) = NULL  

##genotype file 
for(i in 1:20000){
	write.table(output[i,], "geno_real_20000.txt", append = TRUE, row.names = FALSE, quote = FALSE)
}

##output for annotation file
anno <- data.frame(matrix(ncol = 1, nrow = 20000))
for(i in 1:20000){
	anno[i, ] = paste0(chr[i], sample(1:20, 1), " ", pos[i], " ", chr[i])
}
colnames(anno) = NULL

##annotation file 
for(i in 1:34000){
	write.table(anno[i,], "anno_real_20000.txt", append = TRUE, row.names= FALSE, quote = FALSE)
}

details <- read.csv("sucker_metadata_knownsex_89.csv", header=T)
sex <- details$sex 
sex <- gsub("Male", "0", sex)
sex <- gsub("Female", "1", sex)

##getting sex information on the individuals for pheontype file 
##remember that it randomly adds an X to the beginning of the file 
write.table(sex, "pheno_real_34000.txt", append = TRUE, row.names = FALSE, quote = FALSE)

---------------------------------------------------------------------
##Testing with full data set

##data frame to add genotype likelihood to (GL for each loci & ind)

test_df <- data.frame(matrix(ncol = nind, nrow = 68186))
for(i in 1:68186){
	test_df[i, ] = gp.mean[, i]
}
colnames(test_df) = NULL

variants = read.csv("variants.txt", nrows = 20000, skip = 1, header = F)

chr = variants$V1

pos = variants$V2

##output to add to genotype file (with chromosome information and commas - for GEMMA)
output <- data.frame(matrix(ncol = nind + 1, nrow = 20000))
for(i in 1:20000){
	output[i, 1] = paste0(chr[i], ", ", "x", ", ", "x")
	index = 0 
	for(j in 2:(nind+2)){
		if(j == nind +2){
			output[i,j] = paste0(test_df[i,index])
		}
		else{
		output[i,j] = paste0(test_df[i,index], ", ")
		index = index + 1	
	}}	}
	colnames(output) = NULL  

##genotype file 
for(i in 1:20000){
	write.table(output[i,], "geno_real_20000.txt", append = TRUE, row.names = FALSE, quote = FALSE)
}

##output for annotation file
anno <- data.frame(matrix(ncol = 1, nrow = 20000))
for(i in 1:20000){
	anno[i, ] = paste0(chr[i], sample(1:20, 1), " ", pos[i], " ", chr[i])
}
colnames(anno) = NULL

##annotation file 
for(i in 1:34000){
	write.table(anno[i,], "anno_real_20000.txt", append = TRUE, row.names= FALSE, quote = FALSE)
}
##getting sex information on the individuals for pheontype file 
details <- read.csv("sucker_metadata_knownsex_89.csv", header=T)
sex <- details$sex 
sex <- gsub("Male", "0", sex)
sex <- gsub("Female", "1", sex)

##add to phenotype file 
##remember that it randomly adds an X to the beginning of the file 
write.table(sex, "pheno_real_34000.txt", append = TRUE, row.names = FALSE, quote = FALSE)

