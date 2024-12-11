options(stringsAsFactors=FALSE)
options(scipen = 999)
setwd("/mnt/haoyj/Project/WubinMa_translocation/Translocation_20240311/processing/2.firstbatch/7.star")
chimout1 <- read.table("WM_SGWTS_293T_S_S1_L001.chimOut.bed")
chimreads <- strsplit(chimout1$V4, "\\/")
chimout1$readname <- sapply(chimreads, function(x) x[1])
chimout1$Rdirection <- sapply(chimreads, function(x) x[2])
readuniq <- unique(chimout1$readname)  ## 11211880
kept5pendexpandall <- c()
for(i in sample(1:1000000, 100)){
	chimoutu1 <- chimout1[chimout1$readname %in% readuniq[i], ]
	if(dim(chimoutu1)[1] > 2){
		keptu <- chimoutu1[chimoutu1$Rdirection == as.numeric(names(table(chimoutu1$Rdirection))[table(chimoutu1$Rdirection) == 1]), , drop=F]
		remianr <- chimoutu1[chimoutu1$Rdirection != as.numeric(names(table(chimoutu1$Rdirection))[table(chimoutu1$Rdirection) == 1]), , drop=F]
		argsforremain <- as.data.frame(t(apply(remianr, 1, function(x){
			arg1 <- x[1] == keptu[1, 1]		
			arg2 <- abs(as.numeric(x[2])-keptu[1, 2])
			c(arg1, arg2)
		})))
		if(argsforremain[1, 1] == argsforremain[2, 1]){
			remaink <- remianr[which.max(argsforremain[, 2]), ]
		} else {
			remaink <- remianr[argsforremain[, 1] == 0, ]
		}
		allkeptr <- rbind(keptu, remaink)
	} else {
		allkeptr <- chimoutu1
	}
	allkeptrs <- allkeptr[order(allkeptr$Rdirection, decreasing=F), ]
	kept5pend <- apply(allkeptrs, 1, function(x){
		if(x[6] == "+"){
			c(x[1], as.numeric(x[2]), as.numeric(x[2])+1, x[4], x[5], x[6], x[7])
		} else {
			c(x[1], as.numeric(x[3])-1, as.numeric(x[3]), x[4], x[5], x[6], x[7])
		}
	})
	kept5pendexpand <- c(kept5pend[, 1], kept5pend[, 2])
	kept5pendexpandall <- rbind(kept5pendexpandall, kept5pendexpand)
	print(i)
}
write.table(kept5pendexpandall, file="kept5pendexpandall.293T.sample100.txt", col.names=F, row.names=F, sep="\t", quote=F)
write.table(kept5pendexpandall[, 1:7], file="kept5pendexpandall.293T.sample100_bait.bed", col.names=F, row.names=F, sep="\t", quote=F)
write.table(kept5pendexpandall[, 8:14], file="kept5pendexpandall.293T.sample100_translocation.bed", col.names=F, row.names=F, sep="\t", quote=F)