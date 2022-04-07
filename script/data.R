#read meta data

#read deer meta data
meta_deer <- read.xlsx("./data/meta_deer.xlsx", header=TRUE)
#deer data should contain the following columns: Sample_ID, Date_Collected, County, species, PCR (positive=1, negative=0), Longitude, Latitude)

#read human meta data
meta_human <- read.xlsx("./data/meta_human.xlsx", header=TRUE)

#create combined meta dataset
combined_meta <- meta_deer %>% 
  bind_rows(meta_human)
  
#combined data should contain the following columns: Sample_ID, Date_Collected, County, species
