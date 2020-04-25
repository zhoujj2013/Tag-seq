perl ../prepare_target_offtarget_for_Rideogram.pl ../CTLA4.CAS9.all.sites.merged.confirmed ../CTLA4.CAS9.parsing_water_for_visualization.offtarget.bed > ./offtarget.hg38.txt
Rscript ../run_RIdeogram.R offtarget.hg38.txt test
