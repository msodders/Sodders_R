## R Homework Code 
Step 1: Get Data Into R

	library(tidyverse)
	SNPData <- read_tsv("https://raw.githubusercontent.com/EEOB-BioData/BCB546-Spring2021/main/assignments/UNIX_Assignment/snp_position.txt")
	FangData <- read_tsv("https://raw.githubusercontent.com/EEOB-BioData/BCB546-Spring2021/main/assignments/UNIX_Assignment/fang_et_al_genotypes.txt")

Step 2: Data Inspection

	ncol(SNPData)
	ncol(FangData)
	nrow(SNPData)
	nrow(FangData)
	typeof(SNPData)
	typeof(FangData)
	names(SNPData) 
	names(FangData)

Step 3: Data Processing 

*Separate Maize and Teosinte Data In Fang File and Head to check them*

	FangData_Maize <- FangData[which(FangData$Group=="ZMMIL" | FangData$Group =="ZMMLR" | FangData$Group == "ZMMMR"),]
	head(FangData_Maize)
	
	FangData_Teosinte <-FangData[which(FangData$Group=="ZMPBA" | FangData$Group =="ZMPIL" | FangData$Group == "ZMPJA"),]
	head(FangData_Teosinte)


*Prep SNP data*

	SNPData_SpecificC <- SNPData[,c(1,3,4)]

*Transpose Data*

	MaizeData_transposed <- as.data.frame(t(FangData_Maize))
	TeosinteData_transposed <- as.data.frame(t(FangData_Teosinte))

*Trim Rows and Name Column*

	colnames(MaizeData_transposed) <- as.character(unlist(MaizeData_transposed[1,])) 
	colnames(TeosinteData_transposed) <- as.character(unlist(TeosinteData_transposed[1,]))
	MaizeData_trim <- MaizeData_transposed[-c(1:3),]
	TeosinteData_trim <- TeosinteData_transposed[-c(1:3),]

*Join SNP and Fang Data*

	MaizeData_joined <- merge(x = SNPData_SpecificC, y = MaizeData_trim, by.x = "SNP_ID", by.y ="row.names", all.y = TRUE)

	TeosinteData_joined <- merge(x = SNPData_SpecificC, y = TeosinteData_trim, by.x = "SNP_ID", by.y ="row.names", all.y = TRUE)

*Make SNP Increasing*

	MaizeData_joined <- MaizeData_joined[order(MaizeData_joined$Position),] 
	TeosinteData_joined <- TeosinteData_joined[order(TeosinteData_joined$Position),]
	

*Put Each increase chromosome in its own .csv File*

	for (a in 1:10){
  		MaizeData_Chromosome_Ordered <- MaizeData_joined[MaizeData_joined$Chromosome == a, ]
  		write.csv(MaizeData_Chromosome_Ordered, file= paste("Maize_chromosome", a, ".csv", sep=""), row.names = F)
		}
	for (a in 1:10){
  		TeosinteData_Chromosome_Ordered <- TeosinteData_joined[TeosinteData_joined$Chromosome == a, ]
  		write.csv(TeosinteData_Chromosome_Ordered, file= paste("Teosinte_chromosome", a, ".csv", sep=""), row.names = F)
		}

*Make SNP Decreasing*

	MaizeData_joined_decrease <- MaizeData_joined[order(MaizeData_joined$Position, decreasing=T),]

	TeosinteData_joined_decrease <- TeosinteData_joined[order(TeosinteData_joined$Position, decreasing=T),]
	
*Put Each decrease chromosome in its own .csv File*

	for (a in 1:10){
 		MaizeData_DecreaseOrdered<-MaizeData_joined_decrease[MaizeData_joined_decrease$Chromosome == a, ]
  		write.csv(MaizeData_joined_decrease, file= paste("Maize_chromosome", a, "d.csv", sep=""), row.names = F)
		}
	for (a in 1:10){
 		TeosinteData_DecreaseOrdered<-TeosinteData_joined_decrease[TeosinteData_joined_decrease$Chromosome == a, ]
  		write.csv(TeosinteData_joined_decrease, file= paste("Teosinte_chromosome", a, "d.csv", sep=""), row.names = F)
		}


Step 4: Data Visualization 

*Package Installation/Setup*

	library(ggplot2)
	library(tidyverse)
	library(reshape2)
	library(dplyr)
	library(plyr)

	FangData_transposed2 <- as.data.frame(t(FangData))
	colnames(FangData_transposed2) <- as.character(unlist(FangData_transposed2[1,]))
	FangData_SNPData_joined <- merge(x = SNPData_SpecificC, y = FangData_transposed2, by.x = "SNP_ID", by.y ="row.names", all.y = TRUE)

SNP/Chromosome

	ggplot(FangData_SNPData_joined, aes((Chromosome))) + geom_bar() + ggtitle("SNPs Per Chromosome") + labs(x="Chromosome",y="SNP Number")


Missing data and amount of Heterozygosity

	headers <- colnames(FangData)[-c(1:3)]
	FangData_melted <- melt(FangData, measure.vars = headers)
	FangData_melted[ FangData_melted == "?/?" ] = NA
	FangData_melted$isHomozygous <- (FangData_melted$value=="A/A" | FangData_melted$value=="C/C" | FangData_melted$value=="G/G" | FangData_melted$value=="T/T")

*Ordered*

	FangData_melted_bysample <- FangData_melted[order(FangData_melted$Sample_ID),]
*Counted*

	FangData_bysample_count <- ddply(FangData_melted_bysample, c("Sample_ID"), summarise, counting_homozygous=sum(isHomozygous, na.rm=TRUE), counting_heterozygous=sum(!isHomozygous, na.rm=TRUE), isNA=sum(is.na(isHomozygous)))

	FangData_bysample_count_melted <- melt(FangData_bysample_count, measure.vars = c("counting_homozygous", "counting_heterozygous", "isNA"))

*Order Better for Plotting*

	FangData_grouped <- FangData_melted[order(FangData_melted$Group),]
	FangData_grouped_counted <- ddply(FangData_grouped, c("Group"), summarise, counting_homozygous=sum(isHomozygous, na.rm=TRUE), counting_heterozygous=sum(!isHomozygous, na.rm=TRUE), isNA=sum(is.na(isHomozygous)))
	FangData_grouped_counted_melted <- melt(FangData_grouped_counted, measure.vars = c("counting_homozygous", "counting_heterozygous", "isNA"))
*Plotting2*

	ggplot(FangData_grouped_counted_melted,aes(x = Group, y= value, fill=variable)) + geom_bar(stat = "identity", position = "stack")

Your own visualization- Plot of how many of each group of SNP is contained in FangData

	ggplot(data = FangData_grouped) + geom_density(mapping = aes(x=value), fill="blue")
	

