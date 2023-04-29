# log in synapse ----------------------------------------------------------
install.packages("synapser", repos=c("http://ran.synapse.org", "http://cran.fhcrc.org"), type="win.binary")
library(synapser)
library(filesstrings)

# use your email and password on https://www.synapse.org/#!Synapse:syn4993293
# please apply for access to the data first
synLogin(email="xxx@xxx.xxx",password="XXXXXX")


# download the voice data -------------------------------------------------
for (offset in seq(0,65000,500)){
  results <- synTableQuery(paste0('SELECT * FROM syn5511444 LIMIT 500 OFFSET ',as.character(offset)))
  file_map <- synDownloadTableColumns(results,'audio_audio.m4a')
  for (j in 1:500){
    file.move(file_map[j],f = ,'./mpower_data/voice')
    name <- unlist(strsplit(file_map[j],'\\\\'))
    file.rename(name[[8]],paste0(name[[7]],'.m4a'))
  }
  # put the path to synapse cache in dir.remove function
  # for example, if you account name on Windows platform is Peter, then the path to synapse cache 
  # is 'C:/Users/Peter/Documents/.synapseCache'
  dir.remove('path to synapse Cache')
  print(offset)
}


# convert voice data to wave format ---------------------------------------
# convert m4a format to wav format
# here we use free software FormatFactory http://www.pcfreetime.com/formatfactory/index.php?language=en
# put the path to FormatFactory.exe to pathtoformatfactory variable
# for example, if you install the software under G:\FormatFactory\, then
# pathtoformatfactory <- 'G:\\FormatFactory\\FormatFactory.exe'
pathtoformatfactory <- 'path to Format Factory'

# put the path to voice folder in the last step to pathtovoicefolder variable
# for example, if you download all m4a files under G:\mpower\voice\, then
# pathtovoicefolder <- 'G:\\mpower\\voice\\'
pathtovoicefolder <- 'path to voice folder'

# put the path to converted voice files folder to pathtovoiceoutputfolder variable
# for example, if you save all wav files under G:\mpower\voicewave\, then
# pathtovoiceoutputfolder <- 'G:\\mpower\\voicewave'
pathtovoiceoutputfolder <- 'path to voice output folder'

namel <- list.files()
library(doParallel)
cl <- makePSOCKcluster(10)
registerDoParallel(cl)
foreach(i = 1:length(namel), .combine = c) %dopar% {
  input <- paste0(pathtoformatfactory, ' "-> WAV" "High quality" "', pathtovoicefolder, namel[i],'" "', pathtovoiceoutputfolder, '" /hide')
  system(input)
}
stopCluster(cl)


# extract voice features --------------------------------------------------
# go to './extract voice features' folder 
# continue in Matlab R2022a, open 'extract_voice_feature.m' and follow the instruction

# after running 'extract_voice_feature.m', obtain subid_tmp.csv and voice_tmp.csv under \mpower_data\extracted features folder
voice_F <- read.csv('./mpower_data/extracted features/voice_tmp.csv')
voice_F <- sapply(voice_F, as.numeric)

voice_n <- read.csv('./mpower_data/extracted features/subid_tmp.csv')

colnames(voice_n) <- 'audio_audio.m4a'
voice_F <- cbind(voice_n, voice_F)
voice_F$audio_audio.m4a <- as.numeric(strsplit(voice_F$audio_audio.m4a, '.wav'))

write.csv(voice_F, 'voiceFeatures.csv',row.names = F)


# download the walking data -----------------------------------------------
# use your email and password on https://www.synapse.org/#!Synapse:syn4993293/wiki/247860
# please apply for access to the data first
synLogin(email="xxx@xxx.xxx",password="XXXXXX")

library(filesstrings)
# download the walk outbound
for (offset in seq(0,35500,500)){
  results <- synTableQuery(paste0('SELECT * FROM syn5511449 LIMIT 500 OFFSET ',as.character(offset)))
  file_map <- synDownloadTableColumns(results,'accel_walking_outbound.json.items')
  for (j in 1:length(file_map)){
    file.move(file_map[j],f = ,'./mpower_data/walking')
    name <- unlist(strsplit(file_map[j],'\\\\'))
    file.rename(name[[8]],paste0(name[[7]],'.json'))
  }
  # put the path to synapse cache in dir.remove function
  # for example, if you account name on Windows platform is Peter, then the path to synapse cache 
  # is 'C:/Users/Peter/Documents/.synapseCache'
  dir.remove('path to synapse Cache')
  print(offset)
}

# download the walk rest
for (offset in seq(0,35500,500)){
  results <- synTableQuery(paste0('SELECT * FROM syn5511449 LIMIT 500 OFFSET ',as.character(offset)))
  file_map <- synDownloadTableColumns(results,'accel_walking_rest.json.items')
  for (j in 1:length(file_map)){
    file.move(file_map[j],f = ,'./mpower_data/rest')
    name <- unlist(strsplit(file_map[j],'\\\\'))
    file.rename(name[[8]],paste0(name[[7]],'.json'))
  }
  # put the path to synapse cache in dir.remove function
  # for example, if you account name on Windows platform is Peter, then the path to synapse cache 
  # is 'C:/Users/Peter/Documents/.synapseCache'
  dir.remove('path to synapse Cache')
  print(offset)
}


# extract walking features ------------------------------------------------
source('./extract walking features/walking.R')
setwd('./mpower_data/walking')
namel <- list.files()
walk_F <- NULL
for (i in 1:length(namel)){
  walk_F <- rbind(walk_F, 
                  c(as.numeric(unlist(strsplit(namel[i],split = '.json'))), getWalkFeatures(namel[i])))
}
colnames(walk_F)[1] <- "accel_walking_outbound.json.items"
setwd("../../")
write.csv(x = walk_F, './mpower_data/extracted features/walkFeatures.csv', row.names = F)

source('./extract walking features/resting.R')
setwd('./mpower_data/rest')
namel <- list.files()
walk_rest_F <- NULL
for (i in 1:length(namel)){
  walk_rest_F <- rbind(walk_rest_F, 
                       c(as.numeric(unlist(strsplit(namel[i],split = '.json'))), getRestFeatures(namel[i])))
}
colnames(walk_rest_F)[1] <- "accel_walking_rest.json.items"
setwd("../../")
write.csv(x = walk_rest_F, './mpower_data/extracted features/walk_rest_Features.csv', row.names = F)


# download tap features ---------------------------------------------------
# download and save 'tapFeatures.tsv' on https://www.synapse.org/#!Synapse:syn4993293/files/ to './mpower_data/extracted features/' as 'tapFeatures.csv'.


# other data --------------------------------------------------------------
voice <- as.data.frame(synTableQuery("select * from syn5511444"))
write.csv(voice, './mpower_data/other/Voice.csv')

updrs <- as.data.frame(synTableQuery("select * from syn5511432"))
write.csv(updrs, './mpower_data/other/UPDRS.csv')

demo <-  as.data.frame(synTableQuery("SELECT * FROM syn5511429"))
write.csv(demo, './mpower_data/other/demographics.csv')

tapping <- as.data.frame(synTableQuery("SELECT * FROM syn5511439"))
write.csv(tapping, './mpower_data/other/tapping.csv')

mem <- as.data.frame(synTableQuery("SELECT * FROM syn5511434"))
write.csv(mem, './mpower_data/other/Memory.csv')

walk <- as.data.frame(synTableQuery("SELECT * FROM syn5511449"))
write.csv(walk, './mpower_data/other/walking.csv')
