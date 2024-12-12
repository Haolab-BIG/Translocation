options(stringsAsFactors=FALSE)
options(scipen = 999)
#
write.row <- function(x, file = filepath, append = TRUE, quote = F, sep="\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"), fileEncoding = ""){
    write.table(as.matrix(t(x)), file, append , quote, sep, eol, na, dec, row.names, col.names, qmethod, fileEncoding)
}
#
setwd("/mnt/haoyj/Project/WubinMa_translocation/Translocation_20240311/processing/2.firstbatch/7.star")
chimout1 <- read.table("WM_SGWTS_293T_S_S1_L001.chimOut.bed", header=F, sep="\t", quote="")
chimreads <- strsplit(chimout1$V4, "\\/")
chimout1$readname <- sapply(chimreads, function(x) x[1])
chimout1$Rdirection <- sapply(chimreads, function(x) x[2])
readuniq <- unique(chimout1$readname)  ## 11211880
junctiontxt <- read.table("WM_SGWTS_293T_S_S1_L0012Chimeric.out.junction.filter", header=F, sep="\t", quote="")
reads <- unique(junctiontxt$V10)
#
total_range <- length(reads)
num_parts <- 50
part_size <- floor(total_range / num_parts)
breaks <- seq(1, total_range, by = part_size)
if (length(breaks) < num_parts + 1) {
  breaks <- c(breaks, total_range + 1)
}
groups <- cut(1:total_range, breaks = breaks, include.lowest = TRUE, labels = FALSE)
ii <- 1
start <- breaks[ii]
end <- breaks[ii + 1] - 1
#
#selchimoreadsdfs <- c()
for(i in start:end){
	chimout1single <- chimout1[chimout1$readname %in% reads[i], ]
	junctionsingle <- junctiontxt[junctiontxt$V10 %in% reads[i], ]
	indexn <- unique(chimout1single$V5)
	chimoreads <- c()
	for(j in indexn){
		chimoutu1 <- chimout1single[chimout1single$V5 == indexn[j], ]
		chimoutu1[, 2] <- as.numeric(gsub(" ", "", chimoutu1[, 2]))
		chimoutu1[, 3] <- as.numeric(gsub(" ", "", chimoutu1[, 3]))
		if(dim(chimoutu1)[1] > 2){
			keptu <- chimoutu1[chimoutu1$Rdirection == as.numeric(names(table(chimoutu1$Rdirection))[table(chimoutu1$Rdirection) == 1]), , drop=F]
			remianr <- chimoutu1[chimoutu1$Rdirection != as.numeric(names(table(chimoutu1$Rdirection))[table(chimoutu1$Rdirection) == 1]), , drop=F]
			argsforremain <- as.data.frame(t(apply(remianr, 1, function(x){
				arg1 <- x[1] == keptu[1, 1]		
				if(x[9] == 1){
					arg2 <- abs(as.numeric(x[2])-junctionsingle[1, 2])
				}
				if(x[9] == 2){
					arg2 <- abs(as.numeric(x[2])-junctionsingle[1, 5])
				}
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
		allkeptrs <- c(t(allkeptr[order(allkeptr$Rdirection, decreasing=F), ]))
		chimoreads <- rbind(chimoreads, allkeptrs)
	}
	chimoreadsdf <- as.data.frame(chimoreads)
	chimoreadsdf[, 2] <- as.numeric(gsub(" ", "", chimoreadsdf[, 2]))
	chimoreadsdf[, 3] <- as.numeric(gsub(" ", "", chimoreadsdf[, 3]))
	chimoreadsdf[, 11] <- as.numeric(gsub(" ", "", chimoreadsdf[, 11]))
	chimoreadsdf[, 12] <- as.numeric(gsub(" ", "", chimoreadsdf[, 12]))
	selchimoreadsdf <- chimoreadsdf[apply(chimoreadsdf, 1, function(x){
		#x[7] <- gsub("H", "S", x[7])
		#x[16] <- gsub("H", "S", x[16])
		#junctionsingle[1, 12] <- gsub("H", "S", junctionsingle[1, 12])
		#junctionsingle[1, 16] <- gsub("H", "S", junctionsingle[1, 16])
		cigar_str1 <- x[7]
		cigar_str2 <- x[16]
		matches1 <- unlist(regmatches(cigar_str1, gregexpr("\\d+[A-Za-z]", cigar_str1)))
		matches11 <- matches1[grep("M", matches1)]
		#print(matches11)
		matches2 <- unlist(regmatches(cigar_str2, gregexpr("\\d+[A-Za-z]", cigar_str2)))
		matches22 <- matches2[grep("M", matches2)]
		#print(matches22)
		matchl1 <- length(grep(matches11, junctionsingle[1, 12], fixed=T)) == 1 | length(grep(matches11, junctionsingle[1, 14], fixed=T)) == 1
		matchl2 <- length(grep(matches22, junctionsingle[1, 14], fixed=T)) == 1 | length(grep(matches22, junctionsingle[1, 12], fixed=T)) == 1
		x[1] == junctionsingle[1, 1] & x[10] == junctionsingle[1, 4] & matchl1 & matchl2
	}), ]
	print(dim(selchimoreadsdf))
	print(i)
	kept5pend <- apply(selchimoreadsdf, 1, function(x){
		if(x[6] == "+"){
			x[3] <- as.numeric(x[2])+1
		}
		if(x[6] == "-"){
			x[2] <- as.numeric(x[3])-1
		}
		if(x[15] == "+"){
			x[12] <- as.numeric(x[11])+1
		}
		if(x[15] == "-"){
			x[11] <- as.numeric(x[12])-1
		}
		x
	})
	#selchimoreadsdfs <- rbind(selchimoreadsdfs, kept5pend)
	write.row(kept5pend, file=paste0("WM_SGWTS_293T_S_S1_L001.chimOut.best.jun.1bp.", ii, ".txt"))
}

#######
options(stringsAsFactors=FALSE)
setwd("/mnt/haoyj/Project/WubinMa_translocation/Translocation_20240311/processing/2.firstbatch/run")
rr <- read.table("2.starforchimeric.filter.1.r", header=F, sep="*", quote="")
for(ind in 3:50){
	rr[22, ] <- gsub(1, ind, rr[22, ])
	write.table(rr, file=paste0("2.starforchimeric.filter.", ind, ".r"), col.names=F, row.names=F, sep="*", quote=F)
	system(paste0("R CMD BATCH /mnt/haoyj/Project/WubinMa_translocation/Translocation_20240311/processing/2.firstbatch/run/2.starforchimeric.filter.", ind, ".r"), wait=F)
}