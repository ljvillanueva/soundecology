
RMySQL_user=""
RMySQL_password=""
RMySQL_host="localhost"
RMySQL_dbname="tip2008"

#Load packages
require(ggplot2)
require(RMySQL)
drv <- dbDriver("MySQL")
#Connect to MySQL
con <- dbConnect(drv, user=RMySQL_user, password = RMySQL_password, host=RMySQL_host, dbname = RMySQL_dbname)


#Execute query
res <- dbSendQuery(con, "SELECT index_data.adi AS adi, index_data.shannon AS shannon, Sites.SiteID AS SiteID FROM index_data, Sounds, Sites WHERE index_data.SoundID=Sounds.SoundID and Sounds.SiteID=Sites.SiteID")
#
#Assign results to data1
data1 <- fetch(res, n = -1)
#
#Get number of records
no_records<-dim(data1)[1]
#

cor.results <- data.frame(matrix(NA, nrow = 9, ncol = 4))

#ALL DATA
cat("All Sites:")
cat("\n")
a <- cor.test(data1$adi, data1$shannon)
cor.results[9,] <- c("All", a$estimate, a$statistic, a$p.value)

sink(file="results.txt", append=TRUE)
cat("All Sites:")
cat("\n")
print(cor.test(data1$adi, data1$shannon))
cat("\n\n")
sink()

png(file=paste("plots/all_", format(Sys.time(), "%d%b%Y_%H%M%S"), ".png", sep=""), height=800, width=1200)

p <- ggplot(data1, aes(adi, shannon))
this_plot <- p + geom_point(aes(colour = factor(SiteID)))
print(this_plot)
dev.off()


#By Site

#Execute query
res <- dbSendQuery(con, "SELECT * from Sites")
#
#Assign results to data1
sites <- fetch(res, n = -1)
#
#Get number of records
no_sites<-dim(sites)[1]
#
for (i in 1:no_sites) {
  this_siteID = sites[i,1]
  this_siteName = sites[i,2]
  cat(paste("Site:", this_siteName, "\n"))

  #Execute query
  res <- dbSendQuery(con, paste("SELECT index_data.adi, index_data.shannon from index_data, Sounds WHERE index_data.SoundID=Sounds.SoundID AND Sounds.SiteID='", i, "'", sep=""))
  #
  #Assign results to data1
  data1 <- fetch(res, n = -1)
  
  a <- cor.test(data1$adi, data1$shannon)
  cor.results[i,] <- c(this_siteName, a$estimate, a$statistic, a$p.value)
  
  print(a)
  
  sink(file="results.txt", append=TRUE)
  cat("Site:")
  cat(this_siteName)
  cat("\n")
  print(a)
  cat("\n\n")
  sink()
  
  png(file=paste("plots/", this_siteID, "_", format(Sys.time(), "%d%b%Y_%H%M%S"), ".png", sep=""), height=800, width=1200)
  
  p <- ggplot(data1, aes(adi, shannon))
  this_plot <- p + geom_point()
  print(this_plot)
  dev.off()
  
}

dbDisconnect(con)

names(cor.results) <- c("Site", "Estimate", "Statistic", "P Value")

write.csv(cor.results, file="results.table.csv", append=FALSE, row.names=FALSE)
