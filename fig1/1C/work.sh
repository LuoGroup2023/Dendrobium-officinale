cd /a/3_tiaozheng2/coverage/hifi
minimap2 -t 60 -ax map-hifi -L  --secondary=no  -K8G /a/2_gapfill/juicer_tgsgapcloser_hifi6/Allhifi_6.scaff_seqs /b/ALL_HiFi.fq > map.sam
cd /a/3_tiaozheng2/coverage/ont
minimap2 -t 60 -ax map-ont -L  --secondary=no  -K8G /a/2_gapfill/juicer_tgsgapcloser_hifi6/Allhifi_6.scaff_seqs /c/pass.all.fq > map.sam
cd /a/3_tiaozheng2/coverage/ul
minimap2 -t 60 -ax map-ont -L  --secondary=no  -K8G /a/2_gapfill/juicer_tgsgapcloser_hifi6/Allhifi_6.scaff_seqs "/c/pass.ul.fq.gz" > map.sam


## samtools view -@ 80 -b -S map.sam -o map.bam
## samtools sort -@ 80 -o  map_sorted.bam -T  hifi_tmp.ali  map.bam
samtools sort -@ 60 -O BAM map.sam > map_sorted.bam
samtools index -@ 60 map_sorted.bam

perl /a/3_coverage-anomalous/get_the_detailed_mapping_information_from_bwa_bam_file.pl map_sorted.bam 0 > map.bam.infor
perl /a/3_coverage-anomalous/get_mapped_bed.pl map.bam.infor mapped.bed maybemapped.bed 
cat mapped.bed| sort -k1,1 -k2,2n - > allCellMapped.bed
cat maybemapped.bed | sort -k1,1 -k2,2n - > allCellMaybeMapped.bed



grep "^chr" allCellMapped.bed > mapped.sort.bed
cat /a/2_gapfill/juicer_tgsgapcloser_hifi6/Allhifi_6.scaff_seqs.fai |cut -f 1,2 |head -n 19 > Dof_2.chr19s_len
bedtools genomecov -i mapped.sort.bed -g Dof_2.chr19s_len -bga -split > Do2.filtered.depth.bed
perl /a/3_tiaozheng2/coverage/make_1kb_depth_line.pl Dof_2.chr19s_len Do2.filtered.depth.bed Do2.filtered.1k.depth.bed
awk '{if ($4==0) print $0}' Do2.filtered.1k.depth.bed |sort -k1.4,1.5n -k1,1 -k2,2n  > coverage_0.bed


cd /a/3_tiaozheng2/coverage/ont
bedtools getfasta -fi /a/2_gapfill/juicer_tgsgapcloser_hifi6/Allhifi_6.scaff_seqs  \
-bed /a/3_tiaozheng2/coverage/ont/coverage_0.bed -fo coverage_0.fa -name

cd /a/3_tiaozheng2/coverage
cp /a/2_gapfill/juicer_tgsgapcloser_hifi6/  .
python /a/3_tiaozheng2/generate_bed.py
bedtools intersect -a  map_sorted.bam  -b /a/3_tiaozheng2/coverage/gap_100000.bed  > gap.region.bam




cd /a/3_tiaozheng/2_coverage2/hifi
samtools depth -@ 24 -a map_sorted.bam > depth.txt
awk '$3 == 0 {print $1"\t"$2-1"\t"$2}' depth.txt > uncovered.bed
bedtools merge -i uncovered.bed > uncovered_merged.bed
cat uncovered_merged.bed |awk 'BEGIN{OFS="\t"}{if ($3-$2 >1000) print $0,$3-$2}' |wc -l  #2458
cat uncovered_merged.bed |awk 'BEGIN{OFS="\t"}{if ($3-$2 >5000) print $0,$3-$2}' |wc -l  #801
cat uncovered_merged.bed |awk 'BEGIN{OFS="\t"}{if ($3-$2 >10000) print $0,$3-$2}' |wc -l  #347

cd /a/3_tiaozheng2/coverage/uncovered
ln -sf /a/3_tiaozheng2/coverage/hifi/uncovered_merged.bed hifi.uncovered_merged.bed
ln -sf /a/3_tiaozheng2/coverage/ul/uncovered_merged.bed ul.uncovered_merged.bed
ln -sf /a/3_tiaozheng2/coverage/ont/uncovered_merged.bed ont.uncovered_merged.bed
cat *bed > 3typs_Data.uncovered_merged.bed
cd /a/3_tiaozheng2/coverage/uncovered/10k
grep chr ../3typs_Data.uncovered_merged.merge.bed|awk 'BEGIN{OFS="\t"}{if ($3-$2 >10000) print $0,$3-$2}' > typs_Data.uncovered_merged.merge.10k.bed
