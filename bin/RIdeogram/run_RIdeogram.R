args <- commandArgs(TRUE)

if (length(args) < 2) {
	stop("Please input: xx.bed prefix")
}

require(RIdeogram)

bed=args[1]
prefix = args[2]

human_karyotype <- read.table("karyotype.txt", sep = "\t", header = T, stringsAsFactors = F)
gene_density <- read.table("genedensity.txt", sep = "\t", header = T, stringsAsFactors = F)
Random_RNAs_500 <- read.table(bed, sep = "\t", header = T, stringsAsFactors = F)

ideogram(karyotype = human_karyotype, overlaid = gene_density, label = Random_RNAs_500, label_type = "marker",output = paste(prefix,"offtarget.chromosome.svg",sep="."))
convertSVG(paste(prefix,"offtarget.chromosome.svg",sep="."), file = paste(prefix,"offtarget.chromosome",sep="."), device = "png")
#wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz
#less -S cytoBand.txt.hg19.gz | grep --color=auto acen | less -S
#https://cran.r-project.org/web/packages/RIdeogram/vignettes/RIdeogram.html
# ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/GRCh37_mapping/gencode.v32lift37.annotation.gff3.gz
#options(scipen=200)
#write.table(gene_density, file="human_genedensity.hg19.txt", quote=F, sep="\t", row.names=F)
#gene_density <- GFFex(input = "./gencode.v32lift37.annotation.gff3.gz", karyotype = "./human_karyotype.hg19.txt", feature = "gene", window = 1000000) 
