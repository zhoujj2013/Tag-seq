less human_karyotype.hg38.txt | cut -f 1,4,5 | sed '1d' | sed 's/^/chr/g' > human_cent.hg38.bed
CrossMap.py bed hg38ToHg19.over.chain.gz human_cent.hg38.bed human_cent.hg19.bed

paste hg19.len.txt acen.hg19.merged.txt | awk '{print $1"\t0\t"$2"\t"$4"\t"$5}' | cat human_karyotype.hg38.txt.header - > human_karyotype.hg19.txt


