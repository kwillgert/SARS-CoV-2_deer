#select sequences to be used for tree from a ML tree

#load sequences
seqs<- read.dna(file="./data/phylo_20211104/iqtree_all_vcf/all_vcf_alignment-2021-11-04_06-21-47.fasta", format = "fasta", as.character=T) %>% 
  as.data.frame() %>% 
  rownames_to_column("Sample_ID")
#seqs<-adegenet::fasta2DNAbin("./data/phylo_20211104/iqtree_all_vcf/all_vcf_alignment-2021-11-04_06-21-47.fasta")


#remove sequences without meta data
#create meta file for sequences
phylo_meta<-seqs %>% select(Sample_ID) %>% 
  left_join(combined_meta) %>% 
  replace_na(list(species="various"))

#write.xlsx(phylo_meta, "./data/phylo_meta.xlsx", sheetName="phylo_meta")

#identify tip samples without a match in meta data
missing_meta<-phylo_meta %>% 
  filter(!Sample_ID %in% combined_meta$Sample_ID) 

missing_meta$Sample_ID

#remove samples without meta data and append collection dates to sequence names
seqs<- seqs %>% 
  filter(!Sample_ID %in% missing_meta$Sample_ID) %>%
  left_join(combined_meta %>% select(Sample_ID, Date_Collected)) %>% 
  #mutate(Date_Collected=decimal_date(Date_Collected)) %>% #using decimal date
  mutate(Date_Collected=as.numeric(Date_Collected-min(Date_Collected, na.rm=TRUE))) %>%  #using difference in time between collection events (days)
  mutate(Sample_ID=str_c(Sample_ID, Date_Collected, sep="|")) %>% #add collection date to sequence name
  #mutate(Sample_ID=str_c(Sample_ID, Date_Collected, sep="-")) %>% #add collection date to sequence name
  arrange(Date_Collected) %>% 
  select(-Date_Collected) 

#update meta data for phylogeny
phylo_meta<-seqs %>% select(Sample_ID) %>% 
  left_join(combined_meta %>% mutate(Sample_ID=str_c(Sample_ID, Date_Collected, sep="|"))) %>% 
  replace_na(list(species="various"))

seqs<-seqs %>% 
  column_to_rownames(var="Sample_ID")#make sequence name row name


#convert sequences to matrix
seqs<-as.matrix(seqs)

#convert to DNAbin
seqs<-as.DNAbin(seqs)

#remove missing data (n) from sequence??? - see exercise 3


#save data as fasta file to use to build timed ML tree in IQ-Tree2
write.FASTA(seqs, "./data/covid19_seqs.fasta", header=FALSE)







#create fasta file with human samples only with date in decimal years
seqs2<- read.dna(file="./data/covid19_seqs.fasta", format = "fasta", as.character=T) %>% 
  as.data.frame() %>% 
  rownames_to_column("Sample_ID") %>% 
  mutate(Sample_ID=str_remove(Sample_ID, pattern="\\|.*$")) # remove sampling day at end of Sample ID to match with meta data

#remove samples without meta data and append collection dates to sequence names
seqs2<- seqs2 %>% 
  left_join(combined_meta %>% select(Sample_ID, Date_Collected, species)) %>% 
  filter(species=="human") %>% #only include human samples
  mutate(Date_Collected=decimal_date(Date_Collected)) %>% #using decimal date
  #mutate(Date_Collected=as.numeric(Date_Collected-min(Date_Collected, na.rm=TRUE))) %>%  #using difference in time between collection events (days)
  mutate(Sample_ID=str_c(Sample_ID, Date_Collected, sep="|")) %>% #add collection date to sequence name
  #mutate(Sample_ID=str_c(Sample_ID, Date_Collected, sep="-")) %>% #add collection date to sequence name
  arrange(Date_Collected) %>% 
  select(-Date_Collected, -species) #remove date collected and species

seqs2<-seqs2 %>% 
  column_to_rownames(var="Sample_ID")#make sequence name row name


#convert sequences to matrix
seqs2<-as.matrix(seqs2)

#convert to DNAbin
seqs2<-as.DNAbin(seqs2)



#save data as fasta file to use to run in BEAST
write.FASTA(seqs2, "./data/covid19_human_decDate.fasta", header=FALSE)

###NOTE: write.FASTA adds "FALSE" at top of text document, needs to be manually removed


seqs2
class(seqs2)

#compute distances
D<-dist.dna(seqs, model="TN93")
class(D)
length(D)

#Plot a heatmap of your pairwise distances
temp <- as.data.frame(as.matrix(D)) # turn D into a matrix re-shape it correctly
table.paint(temp, cleg = 0, clabel.row = 0.5, clabel.col = 0.5)
View(temp)
temp <- t(as.matrix(D))
temp <- temp[, ncol(temp):1]
par(mar = c(1, 5, 5, 1))
image(x = 1:nrow(seqs), y = 1:nrow(seqs), temp, col = rev(heat.colors(100)),
      xaxt = "n", yaxt = "n", xlab = "", ylab = "")
axis(side = 2, at = 1:nrow(seqs), lab = rownames(seqs), las = 2, cex.axis = 0.5)
axis(side = 3, at = 1:nrow(seqs), lab = rownames(seqs), las = 3, cex.axis = 0.5)

dna2 <- as.phyDat(seqs)
class(dna2)
dna2

