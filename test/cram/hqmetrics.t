

  $ BH_ROOT=$TESTDIR/../../

  $ bauhaus -o map -m -t ${BH_ROOT}test/data/hqexperiment.csv -w HQMetrics --noGrid generate
  Validation and input resolution succeeded.
  Runnable workflow written to directory "map"


  $ tree map
  map
  |-- R
  |   `-- hqrf_pbbamr.R
  |-- build.ninja
  |-- condition-table.csv
  |-- log
  `-- run.sh
  
  2 directories, 4 files


  $ cat map/condition-table.csv
  Condition,SubreadSet,Genome
  a,/pbi/collections/315/3150477/r54003_20161109_183934/2_B01/m54003_161110_020238.subreadset.xml,ecoliK12_pbi_March2013


  $ cat map/build.ninja
  # Variables
  ncpus = 8
  scratchDir = /scratch
  grid = 
  gridSMP = 
  
  # Rules
  rule hqrf_pbbamr
    command = $grid Rscript --verbose R/hqrf_pbbamr.R $origbam xscraps $in $out
  
  rule bam2bam_nohq
    command = $grid bam2bam -j $ncpus -b $ncpus --fullHQ --adapters $adapters $
        -o $outPrefix $in
  
  rule subreadset_create
    command = $grid dataset create $out $in
  
  rule map
    command = $gridSMP pbalign  --tmpDir=$scratchDir --nproc $ncpus $in $
        $reference $out
  
  rule splitByZmw
    command = $grid dataset split --zmws --targetSize 1 --chunks 8 --outdir $
        $outdir $in
  
  rule mergeDatasetsForCondition
    command = $grid dataset merge $out $in
  
  
  # Build targets
  build a/nohq/m54003_161110_020238.subreads.bam: bam2bam_nohq $
      /pbi/collections/315/3150477/r54003_20161109_183934/2_B01/m54003_161110_020238.subreadset.xml
    outPrefix = a/nohq/m54003_161110_020238
    adapters = $
        /pbi/collections/315/3150477/r54003_20161109_183934/2_B01/m54003_161110_020238.adapters.fasta
  
  build a/nohq/m54003_161110_020238.subreadset.xml: subreadset_create $
      a/nohq/m54003_161110_020238.subreads.bam
  
  build a/subreads_chunks/m54003_161110_020238.chunk0.subreadset.xml $
      a/subreads_chunks/m54003_161110_020238.chunk1.subreadset.xml $
      a/subreads_chunks/m54003_161110_020238.chunk2.subreadset.xml $
      a/subreads_chunks/m54003_161110_020238.chunk3.subreadset.xml $
      a/subreads_chunks/m54003_161110_020238.chunk4.subreadset.xml $
      a/subreads_chunks/m54003_161110_020238.chunk5.subreadset.xml $
      a/subreads_chunks/m54003_161110_020238.chunk6.subreadset.xml $
      a/subreads_chunks/m54003_161110_020238.chunk7.subreadset.xml: $
      splitByZmw a/nohq/m54003_161110_020238.subreadset.xml
    outdir = a/subreads_chunks
  
  build a/mapping_chunks/m54003_161110_020238.chunk0.alignmentset.xml: map $
      a/subreads_chunks/m54003_161110_020238.chunk0.subreadset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/ecoliK12_pbi_March2013/sequence/ecoliK12_pbi_March2013.fasta
  
  build a/mapping_chunks/m54003_161110_020238.chunk1.alignmentset.xml: map $
      a/subreads_chunks/m54003_161110_020238.chunk1.subreadset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/ecoliK12_pbi_March2013/sequence/ecoliK12_pbi_March2013.fasta
  
  build a/mapping_chunks/m54003_161110_020238.chunk2.alignmentset.xml: map $
      a/subreads_chunks/m54003_161110_020238.chunk2.subreadset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/ecoliK12_pbi_March2013/sequence/ecoliK12_pbi_March2013.fasta
  
  build a/mapping_chunks/m54003_161110_020238.chunk3.alignmentset.xml: map $
      a/subreads_chunks/m54003_161110_020238.chunk3.subreadset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/ecoliK12_pbi_March2013/sequence/ecoliK12_pbi_March2013.fasta
  
  build a/mapping_chunks/m54003_161110_020238.chunk4.alignmentset.xml: map $
      a/subreads_chunks/m54003_161110_020238.chunk4.subreadset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/ecoliK12_pbi_March2013/sequence/ecoliK12_pbi_March2013.fasta
  
  build a/mapping_chunks/m54003_161110_020238.chunk5.alignmentset.xml: map $
      a/subreads_chunks/m54003_161110_020238.chunk5.subreadset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/ecoliK12_pbi_March2013/sequence/ecoliK12_pbi_March2013.fasta
  
  build a/mapping_chunks/m54003_161110_020238.chunk6.alignmentset.xml: map $
      a/subreads_chunks/m54003_161110_020238.chunk6.subreadset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/ecoliK12_pbi_March2013/sequence/ecoliK12_pbi_March2013.fasta
  
  build a/mapping_chunks/m54003_161110_020238.chunk7.alignmentset.xml: map $
      a/subreads_chunks/m54003_161110_020238.chunk7.subreadset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/ecoliK12_pbi_March2013/sequence/ecoliK12_pbi_March2013.fasta
  
  build a/mapping/all_movies.alignmentset.xml: mergeDatasetsForCondition $
      a/mapping_chunks/m54003_161110_020238.chunk0.alignmentset.xml $
      a/mapping_chunks/m54003_161110_020238.chunk1.alignmentset.xml $
      a/mapping_chunks/m54003_161110_020238.chunk2.alignmentset.xml $
      a/mapping_chunks/m54003_161110_020238.chunk3.alignmentset.xml $
      a/mapping_chunks/m54003_161110_020238.chunk4.alignmentset.xml $
      a/mapping_chunks/m54003_161110_020238.chunk5.alignmentset.xml $
      a/mapping_chunks/m54003_161110_020238.chunk6.alignmentset.xml $
      a/mapping_chunks/m54003_161110_020238.chunk7.alignmentset.xml
  
  build a/m54003_161110_020238.hqrm.pdf $
      a/m54003_161110_020238.hqrm_metrics.csv: hqrf_pbbamr $
      a/mapping/all_movies.alignmentset.xml
    origbam = $
        /pbi/collections/315/3150477/r54003_20161109_183934/2_B01/m54003_161110_020238.subreads.bam
  

