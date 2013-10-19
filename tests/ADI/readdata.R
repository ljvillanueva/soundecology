RMySQL_user=""
RMySQL_password=""
RMySQL_host="localhost"
RMySQL_dbname="tip2008"

#Load packages
require(RMySQL)
drv <- dbDriver("MySQL")
#Connect to MySQL
con <- dbConnect(drv, user=RMySQL_user, password = RMySQL_password, host=RMySQL_host, dbname = RMySQL_dbname)


csv_files <- dir(path = "results_tip2008", pattern="csv$", ignore.case=TRUE)

index_data <- data.frame()
index_data50 <- data.frame()

for (i in 1:length(csv_files)){
#for (i in 1:100){
  
	SoundID = unlist(strsplit(csv_files[i], split=".indices.csv"))
	index_data <- cbind(SoundID, read.table(paste("results_tip2008/", csv_files[i], sep=""), sep=",", 
				   col.names=c("adi", "shannon", "db"), 
				   fill=FALSE, 
				   strip.white=TRUE))

	dbWriteTable(con, "index_data", index_data, row.names = F, overwrite = F, append=T)
	
}

dbDisconnect(con)
