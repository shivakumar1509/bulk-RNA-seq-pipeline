

#######################
### Run DESEQ2 in R ###
#######################

cd ~/Desktop/7_Deseq2_ready_2_30SR
mkdir ../8_Volcano_and_QC_Plots
mkdir ../8_Volcano_and_QC_Plots/QC_Plots
mv ./*allsamples.pdf ~/Desktop/8_Volcano_and_QC_Plots/QC_Plots
mkdir ../8_Volcano_and_QC_Plots/Volcano_Plots


#### Move Volcano Plots to proper folder
mv ./*.pdf ~/Desktop/8_Volcano_and_QC_Plots/Volcano_Plots

#### Move CSVs to proper folder
mkdir ../9_Deseq2_Results_CSV
mv ./*.csv ../9_Deseq2_Results_CSV

### Chanege directory to where working 
cd ~/Desktop/9_Deseq2_Results_CSV

###Convert CSV in Tab Delimited txt files
for file in ./*.csv; do echo ${file}; cat ${file} | tr "," "\\t" > ${file}.txt; done
mkdir ../9b_Deseq2_Results_txt
mv ./*.txt ../*Deseq2_Results_txt
cd ../*Deseq2_Results_txt
for filename in ./*.txt; do [ -f "$filename" ] | mv "$filename" "${filename//.csv/}"; done

### Make File with only normalized read counts
cat R2054Q_vs_Bl6.deseq2.txt| awk -v OFS='\t' '{print $2, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20, $21, $22, $23, $24, $25, $26, $27, $28, $29, $30, $31, $32, $33, $34, $35, $36, $37, $38, $39, $40, $41, $42, $43, $44, $45, $46, $47, $48}' > normalized_reads.txt
mkdir ./Normalized_Reads_File
mv normalized_reads.txt ./Normalized_Reads_File

### Trim normalized reads off of Deseq2 file to only contain stats
for file in ./*.txt; do echo ${file}; cat ${file} | awk -v OFS='\t' '{print $2,$3,$4,$5,$6,$7,$8}' | sed 's/"//g' > ${file}.trimmed.txt; done
mkdir ../9c_Deseq2_Results_txt_trimmed
mv ./*trimmed.txt ../*Deseq2_Results_txt_trimmed
cd ../*Deseq2_Results_txt_trimmed
for filename in ./*.txt; do [ -f "$filename" ] | mv "$filename" "${filename//.deseq2.txt/}"; done

### Make Trimmed Files with Gene Symnbol Column
mkdir ../6_HTSEQ_counts
cd ../6_HTSEQ_counts
touch HTSEQ.txt
echo -e 'Ensembl \t Gene'  >> HTSEQ.txt
cat Bl6_M01.sam.HTSEQ.txt| awk '{print $1 "\t" $2}' | grep -v '^__no_feature' | grep -v '^__ambiguous' | grep -v '^__too_low_aQual' | grep -v '^__not_aligned' | grep -v '^__alignment_not_unique' >> HTSEQ.txt
mkdir ../9d_Deseq2_Results_txt_trimmed_with_Genesymbol
cd ../*Deseq2_Results_txt_trimmed
for file in *.txt; do echo ${file} | paste -d' ' ${file} ../6_HTSEQ_counts/HTSEQ.txt > ${file}_temp.txt; done
rm ../6_HTSEQ_counts/HTSEQ.txt
for file in *temp.txt; do echo ${file}; cat ${file} | awk -v OFS='\t' '{print $1,$2,$3,$4,$5,$6,$7,$8,$9}' > ${file}.genesymbol.txt; done
mv ./*genesymbol.txt ../*Deseq2_Results_txt_trimmed_with_Genesymbol
rm *temp*
cd ../9d_Deseq2_Results_txt_trimmed_with_Genesymbol
for filename in ./*.txt; do [ -f "$filename" ] | mv "$filename" "${filename//txt_temp.txt./}"; done

### Filter out only significant genes with a Pvalue of 0.05 or less
for file in ./*genesymbol.txt; do echo ${file}; cat ${file} | awk 'NR==1{print}' > ${file}.SIG.txt; done
for file in ./*genesymbol.txt; do echo ${file}; cat ${file} | awk '{ if ($7 <= 0.05) print $0}' >> ${file}.SIG.txt; done
for filename in ./*SIG.txt; do [ -f "$filename" ] | mv "$filename" "${filename//trimmed.genesymbol.txt./}"; done
mkdir ../10a_Sig_Genes_Pvalue_0.05
mv ./*.SIG* ../*Sig_Genes_Pvalue*
cd ../*Sig_Genes_Pvalue*

### From list of significant genes, filter out further by either positively or negatively regulated
for file in ./*SIG.txt; do echo ${file}; cat ${file} | awk 'NR==1{print}' > ${file}.Neg.txt; done
for file in ./*SIG.txt; do echo ${file}; cat ${file} | awk '{ if ($3 < 0) print $0}' >> ${file}.Neg.txt; done
for filename in ./*Neg.txt; do [ -f "$filename" ] | mv "$filename" "${filename//SIG.txt./}"; done
mkdir ../10b_Neg_Sig_Genes_pvalue_0.05
mv *Neg* ../*Neg_Sig_Genes_pvalue*
#for file in ./*SIG.txt; do echo ${file}; cat ${file} | awk 'NR==1{print}' > ${file}.Pos.txt; done
for file in ./*SIG.txt; do echo ${file}; cat ${file} | awk '{ if ($3 > 0) print $0}' >> ${file}.Pos.txt; done
for filename in ./*Pos.txt; do [ -f "$filename" ] | mv "$filename" "${filename//SIG.txt./}"; done
mkdir ../10c_Pos_Sig_Genes_pvalue_0.05
mv *Pos* ../*Pos_Sig_Genes_pvalue*

### Make all SIG, Neg and Pos files have only the genes (without stats) in case want to splice lists apart/together later
for file in ./*SIG.txt; do echo ${file}; cat ${file} | awk '{print $9}' > ${file}._genes_only.txt; done
for filename in ./*genes*; do [ -f "$filename" ] | mv "$filename" "${filename//.txt./}"; done
mkdir ../11a_Sig_Genes_List_Only
mv ./*genes* ../*Sig_Genes_List_Only
cd ../*Neg_Sig_Genes_pvalue*
for file in ./*Neg.txt; do echo ${file}; cat ${file} | awk '{print $9}' > ${file}._genes_only.txt; done
for filename in ./*genes*; do [ -f "$filename" ] | mv "$filename" "${filename//.txt./}"; done
mkdir ../11b_Neg_Sig_Genes_List_Only
mv *genes* ../*Neg_Sig_Genes_List_Only
cd ../*Pos_Sig_Genes_pvalue*
for file in ./*Pos.txt; do echo ${file}; cat ${file} | awk '{print $9}' > ${file}._genes_only.txt; done
for filename in ./*genes*; do [ -f "$filename" ] | mv "$filename" "${filename//.txt./}"; done
mkdir ../11c_Pos_Sig_Genes_List_Only
mv *genes* ../*Pos_Sig_Genes_List_Only

cd ../11a*
find *txt -type f -print0 | xargs -0 sed -i '' /Row.names/d

cd ../11b*
find *txt -type f -print0 | xargs -0 sed -i '' /Row.names/d

cd ../11c*
find *txt -type f -print0 | xargs -0 sed -i '' /Row.names/d

