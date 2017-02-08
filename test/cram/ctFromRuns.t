
Run a bigger workflow--coverage titration with reports.

  $ BH_ROOT=$TESTDIR/../../

  $ bauhaus -o ctFromRuns -m -t ${BH_ROOT}test/data/lambdaAndEcoli-with-p-var.csv -w CoverageTitrationReports generate
  Validation and input resolution succeeded.
  Runnable workflow written to directory "ctFromRuns"


  $ tree ctFromRuns
  ctFromRuns
  |-- build.ninja
  |-- condition-table.csv
  |-- log
  |-- run.sh
  `-- scripts
      `-- R
          `-- coverageTitrationPlots.R
  
  3 directories, 4 files



This ninja file is quite large, but here's hoping any bug will show up
as a small diff.


  $ cat ctFromRuns/build.ninja
  # Variables
  ncpus = 8
  scratchDir = /scratch
  grid = qsub -sync y -cwd -V -b y -e log -o log
  gridSMP = $grid -pe smp $ncpus
  
  # Rules
  rule copySubreadsDataset
    command = $grid dataset create $out $in
  
  rule map
    command = $gridSMP pbalign  --tmpDir=$scratchDir --nproc $ncpus $in $
        $reference $out
  
  rule splitByZmw
    command = $grid dataset split --zmws --targetSize 1 --chunks 8 --outdir $
        $outdir $in
  
  rule mergeDatasetsForCondition
    command = $grid dataset merge $out $in
  
  rule summarize_coverage
    command = $grid python -m $
        pbreports.report.summarize_coverage.summarize_coverage $
        --region_size=10000 $in $reference $out
  
  rule variantCalling
    command = $gridSMP variantCaller $modelPath $modelSpec --algorithm=arrow $
        $coverageLimitArgument -x0 -q0 -j $ncpus --reportEffectiveCoverage $
        $in -r $reference -o $out -o $consensusFasta -o $consensusFastq
  
  rule maskVariantsGff
    command = gffsubtract.pl $in $referenceMask > $out
  
  rule coverageTitrationSummaryAnalysis
    command = Rscript --vanilla scripts/R/coverageTitrationPlots.R .
  
  
  # Build targets
  build EcoliHigh/subreads/m54011_160306_050740.subreadset.xml: $
      copySubreadsDataset $
      /pbi/collections/315/3150122/r54011_20160305_235615/2_B01/m54011_160306_050740.subreadset.xml
  
  build EcoliHigh/subreads_chunks/m54011_160306_050740.chunk0.subreadset.xml $
      EcoliHigh/subreads_chunks/m54011_160306_050740.chunk1.subreadset.xml $
      EcoliHigh/subreads_chunks/m54011_160306_050740.chunk2.subreadset.xml $
      EcoliHigh/subreads_chunks/m54011_160306_050740.chunk3.subreadset.xml $
      EcoliHigh/subreads_chunks/m54011_160306_050740.chunk4.subreadset.xml $
      EcoliHigh/subreads_chunks/m54011_160306_050740.chunk5.subreadset.xml $
      EcoliHigh/subreads_chunks/m54011_160306_050740.chunk6.subreadset.xml $
      EcoliHigh/subreads_chunks/m54011_160306_050740.chunk7.subreadset.xml: $
      splitByZmw EcoliHigh/subreads/m54011_160306_050740.subreadset.xml
    outdir = EcoliHigh/subreads_chunks
  
  build $
      EcoliHigh/mapping_chunks/m54011_160306_050740.chunk0.alignmentset.xml: $
      map EcoliHigh/subreads_chunks/m54011_160306_050740.chunk0.subreadset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/ecoliK12_pbi_March2013/sequence/ecoliK12_pbi_March2013.fasta
  
  build $
      EcoliHigh/mapping_chunks/m54011_160306_050740.chunk1.alignmentset.xml: $
      map EcoliHigh/subreads_chunks/m54011_160306_050740.chunk1.subreadset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/ecoliK12_pbi_March2013/sequence/ecoliK12_pbi_March2013.fasta
  
  build $
      EcoliHigh/mapping_chunks/m54011_160306_050740.chunk2.alignmentset.xml: $
      map EcoliHigh/subreads_chunks/m54011_160306_050740.chunk2.subreadset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/ecoliK12_pbi_March2013/sequence/ecoliK12_pbi_March2013.fasta
  
  build $
      EcoliHigh/mapping_chunks/m54011_160306_050740.chunk3.alignmentset.xml: $
      map EcoliHigh/subreads_chunks/m54011_160306_050740.chunk3.subreadset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/ecoliK12_pbi_March2013/sequence/ecoliK12_pbi_March2013.fasta
  
  build $
      EcoliHigh/mapping_chunks/m54011_160306_050740.chunk4.alignmentset.xml: $
      map EcoliHigh/subreads_chunks/m54011_160306_050740.chunk4.subreadset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/ecoliK12_pbi_March2013/sequence/ecoliK12_pbi_March2013.fasta
  
  build $
      EcoliHigh/mapping_chunks/m54011_160306_050740.chunk5.alignmentset.xml: $
      map EcoliHigh/subreads_chunks/m54011_160306_050740.chunk5.subreadset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/ecoliK12_pbi_March2013/sequence/ecoliK12_pbi_March2013.fasta
  
  build $
      EcoliHigh/mapping_chunks/m54011_160306_050740.chunk6.alignmentset.xml: $
      map EcoliHigh/subreads_chunks/m54011_160306_050740.chunk6.subreadset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/ecoliK12_pbi_March2013/sequence/ecoliK12_pbi_March2013.fasta
  
  build $
      EcoliHigh/mapping_chunks/m54011_160306_050740.chunk7.alignmentset.xml: $
      map EcoliHigh/subreads_chunks/m54011_160306_050740.chunk7.subreadset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/ecoliK12_pbi_March2013/sequence/ecoliK12_pbi_March2013.fasta
  
  build EcoliHigh/mapping/all_movies.alignmentset.xml: $
      mergeDatasetsForCondition $
      EcoliHigh/mapping_chunks/m54011_160306_050740.chunk0.alignmentset.xml $
      EcoliHigh/mapping_chunks/m54011_160306_050740.chunk1.alignmentset.xml $
      EcoliHigh/mapping_chunks/m54011_160306_050740.chunk2.alignmentset.xml $
      EcoliHigh/mapping_chunks/m54011_160306_050740.chunk3.alignmentset.xml $
      EcoliHigh/mapping_chunks/m54011_160306_050740.chunk4.alignmentset.xml $
      EcoliHigh/mapping_chunks/m54011_160306_050740.chunk5.alignmentset.xml $
      EcoliHigh/mapping_chunks/m54011_160306_050740.chunk6.alignmentset.xml $
      EcoliHigh/mapping_chunks/m54011_160306_050740.chunk7.alignmentset.xml
  
  build EcoliLow/subreads/m54011_160305_235923.subreadset.xml: $
      copySubreadsDataset $
      /pbi/collections/315/3150122/r54011_20160305_235615/1_A01/m54011_160305_235923.subreadset.xml
  
  build EcoliLow/subreads_chunks/m54011_160305_235923.chunk0.subreadset.xml $
      EcoliLow/subreads_chunks/m54011_160305_235923.chunk1.subreadset.xml $
      EcoliLow/subreads_chunks/m54011_160305_235923.chunk2.subreadset.xml $
      EcoliLow/subreads_chunks/m54011_160305_235923.chunk3.subreadset.xml $
      EcoliLow/subreads_chunks/m54011_160305_235923.chunk4.subreadset.xml $
      EcoliLow/subreads_chunks/m54011_160305_235923.chunk5.subreadset.xml $
      EcoliLow/subreads_chunks/m54011_160305_235923.chunk6.subreadset.xml $
      EcoliLow/subreads_chunks/m54011_160305_235923.chunk7.subreadset.xml: $
      splitByZmw EcoliLow/subreads/m54011_160305_235923.subreadset.xml
    outdir = EcoliLow/subreads_chunks
  
  build EcoliLow/mapping_chunks/m54011_160305_235923.chunk0.alignmentset.xml: $
      map EcoliLow/subreads_chunks/m54011_160305_235923.chunk0.subreadset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/ecoliK12_pbi_March2013/sequence/ecoliK12_pbi_March2013.fasta
  
  build EcoliLow/mapping_chunks/m54011_160305_235923.chunk1.alignmentset.xml: $
      map EcoliLow/subreads_chunks/m54011_160305_235923.chunk1.subreadset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/ecoliK12_pbi_March2013/sequence/ecoliK12_pbi_March2013.fasta
  
  build EcoliLow/mapping_chunks/m54011_160305_235923.chunk2.alignmentset.xml: $
      map EcoliLow/subreads_chunks/m54011_160305_235923.chunk2.subreadset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/ecoliK12_pbi_March2013/sequence/ecoliK12_pbi_March2013.fasta
  
  build EcoliLow/mapping_chunks/m54011_160305_235923.chunk3.alignmentset.xml: $
      map EcoliLow/subreads_chunks/m54011_160305_235923.chunk3.subreadset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/ecoliK12_pbi_March2013/sequence/ecoliK12_pbi_March2013.fasta
  
  build EcoliLow/mapping_chunks/m54011_160305_235923.chunk4.alignmentset.xml: $
      map EcoliLow/subreads_chunks/m54011_160305_235923.chunk4.subreadset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/ecoliK12_pbi_March2013/sequence/ecoliK12_pbi_March2013.fasta
  
  build EcoliLow/mapping_chunks/m54011_160305_235923.chunk5.alignmentset.xml: $
      map EcoliLow/subreads_chunks/m54011_160305_235923.chunk5.subreadset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/ecoliK12_pbi_March2013/sequence/ecoliK12_pbi_March2013.fasta
  
  build EcoliLow/mapping_chunks/m54011_160305_235923.chunk6.alignmentset.xml: $
      map EcoliLow/subreads_chunks/m54011_160305_235923.chunk6.subreadset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/ecoliK12_pbi_March2013/sequence/ecoliK12_pbi_March2013.fasta
  
  build EcoliLow/mapping_chunks/m54011_160305_235923.chunk7.alignmentset.xml: $
      map EcoliLow/subreads_chunks/m54011_160305_235923.chunk7.subreadset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/ecoliK12_pbi_March2013/sequence/ecoliK12_pbi_March2013.fasta
  
  build EcoliLow/mapping/all_movies.alignmentset.xml: $
      mergeDatasetsForCondition $
      EcoliLow/mapping_chunks/m54011_160305_235923.chunk0.alignmentset.xml $
      EcoliLow/mapping_chunks/m54011_160305_235923.chunk1.alignmentset.xml $
      EcoliLow/mapping_chunks/m54011_160305_235923.chunk2.alignmentset.xml $
      EcoliLow/mapping_chunks/m54011_160305_235923.chunk3.alignmentset.xml $
      EcoliLow/mapping_chunks/m54011_160305_235923.chunk4.alignmentset.xml $
      EcoliLow/mapping_chunks/m54011_160305_235923.chunk5.alignmentset.xml $
      EcoliLow/mapping_chunks/m54011_160305_235923.chunk6.alignmentset.xml $
      EcoliLow/mapping_chunks/m54011_160305_235923.chunk7.alignmentset.xml
  
  build LambdaHigh/subreads/m54008_160308_053311.subreadset.xml: $
      copySubreadsDataset $
      /pbi/collections/315/3150128/r54008_20160308_001811/2_B01/m54008_160308_053311.subreadset.xml
  
  build LambdaHigh/subreads_chunks/m54008_160308_053311.chunk0.subreadset.xml $
      LambdaHigh/subreads_chunks/m54008_160308_053311.chunk1.subreadset.xml $
      LambdaHigh/subreads_chunks/m54008_160308_053311.chunk2.subreadset.xml $
      LambdaHigh/subreads_chunks/m54008_160308_053311.chunk3.subreadset.xml $
      LambdaHigh/subreads_chunks/m54008_160308_053311.chunk4.subreadset.xml $
      LambdaHigh/subreads_chunks/m54008_160308_053311.chunk5.subreadset.xml $
      LambdaHigh/subreads_chunks/m54008_160308_053311.chunk6.subreadset.xml $
      LambdaHigh/subreads_chunks/m54008_160308_053311.chunk7.subreadset.xml: $
      splitByZmw LambdaHigh/subreads/m54008_160308_053311.subreadset.xml
    outdir = LambdaHigh/subreads_chunks
  
  build $
      LambdaHigh/mapping_chunks/m54008_160308_053311.chunk0.alignmentset.xml: $
      map LambdaHigh/subreads_chunks/m54008_160308_053311.chunk0.subreadset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/lambdaNEB/sequence/lambdaNEB.fasta
  
  build $
      LambdaHigh/mapping_chunks/m54008_160308_053311.chunk1.alignmentset.xml: $
      map LambdaHigh/subreads_chunks/m54008_160308_053311.chunk1.subreadset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/lambdaNEB/sequence/lambdaNEB.fasta
  
  build $
      LambdaHigh/mapping_chunks/m54008_160308_053311.chunk2.alignmentset.xml: $
      map LambdaHigh/subreads_chunks/m54008_160308_053311.chunk2.subreadset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/lambdaNEB/sequence/lambdaNEB.fasta
  
  build $
      LambdaHigh/mapping_chunks/m54008_160308_053311.chunk3.alignmentset.xml: $
      map LambdaHigh/subreads_chunks/m54008_160308_053311.chunk3.subreadset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/lambdaNEB/sequence/lambdaNEB.fasta
  
  build $
      LambdaHigh/mapping_chunks/m54008_160308_053311.chunk4.alignmentset.xml: $
      map LambdaHigh/subreads_chunks/m54008_160308_053311.chunk4.subreadset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/lambdaNEB/sequence/lambdaNEB.fasta
  
  build $
      LambdaHigh/mapping_chunks/m54008_160308_053311.chunk5.alignmentset.xml: $
      map LambdaHigh/subreads_chunks/m54008_160308_053311.chunk5.subreadset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/lambdaNEB/sequence/lambdaNEB.fasta
  
  build $
      LambdaHigh/mapping_chunks/m54008_160308_053311.chunk6.alignmentset.xml: $
      map LambdaHigh/subreads_chunks/m54008_160308_053311.chunk6.subreadset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/lambdaNEB/sequence/lambdaNEB.fasta
  
  build $
      LambdaHigh/mapping_chunks/m54008_160308_053311.chunk7.alignmentset.xml: $
      map LambdaHigh/subreads_chunks/m54008_160308_053311.chunk7.subreadset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/lambdaNEB/sequence/lambdaNEB.fasta
  
  build LambdaHigh/mapping/all_movies.alignmentset.xml: $
      mergeDatasetsForCondition $
      LambdaHigh/mapping_chunks/m54008_160308_053311.chunk0.alignmentset.xml $
      LambdaHigh/mapping_chunks/m54008_160308_053311.chunk1.alignmentset.xml $
      LambdaHigh/mapping_chunks/m54008_160308_053311.chunk2.alignmentset.xml $
      LambdaHigh/mapping_chunks/m54008_160308_053311.chunk3.alignmentset.xml $
      LambdaHigh/mapping_chunks/m54008_160308_053311.chunk4.alignmentset.xml $
      LambdaHigh/mapping_chunks/m54008_160308_053311.chunk5.alignmentset.xml $
      LambdaHigh/mapping_chunks/m54008_160308_053311.chunk6.alignmentset.xml $
      LambdaHigh/mapping_chunks/m54008_160308_053311.chunk7.alignmentset.xml
  
  build LambdaLow/subreads/m54008_160308_002050.subreadset.xml: $
      copySubreadsDataset $
      /pbi/collections/315/3150128/r54008_20160308_001811/1_A01/m54008_160308_002050.subreadset.xml
  
  build LambdaLow/subreads_chunks/m54008_160308_002050.chunk0.subreadset.xml $
      LambdaLow/subreads_chunks/m54008_160308_002050.chunk1.subreadset.xml $
      LambdaLow/subreads_chunks/m54008_160308_002050.chunk2.subreadset.xml $
      LambdaLow/subreads_chunks/m54008_160308_002050.chunk3.subreadset.xml $
      LambdaLow/subreads_chunks/m54008_160308_002050.chunk4.subreadset.xml $
      LambdaLow/subreads_chunks/m54008_160308_002050.chunk5.subreadset.xml $
      LambdaLow/subreads_chunks/m54008_160308_002050.chunk6.subreadset.xml $
      LambdaLow/subreads_chunks/m54008_160308_002050.chunk7.subreadset.xml: $
      splitByZmw LambdaLow/subreads/m54008_160308_002050.subreadset.xml
    outdir = LambdaLow/subreads_chunks
  
  build $
      LambdaLow/mapping_chunks/m54008_160308_002050.chunk0.alignmentset.xml: $
      map LambdaLow/subreads_chunks/m54008_160308_002050.chunk0.subreadset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/lambdaNEB/sequence/lambdaNEB.fasta
  
  build $
      LambdaLow/mapping_chunks/m54008_160308_002050.chunk1.alignmentset.xml: $
      map LambdaLow/subreads_chunks/m54008_160308_002050.chunk1.subreadset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/lambdaNEB/sequence/lambdaNEB.fasta
  
  build $
      LambdaLow/mapping_chunks/m54008_160308_002050.chunk2.alignmentset.xml: $
      map LambdaLow/subreads_chunks/m54008_160308_002050.chunk2.subreadset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/lambdaNEB/sequence/lambdaNEB.fasta
  
  build $
      LambdaLow/mapping_chunks/m54008_160308_002050.chunk3.alignmentset.xml: $
      map LambdaLow/subreads_chunks/m54008_160308_002050.chunk3.subreadset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/lambdaNEB/sequence/lambdaNEB.fasta
  
  build $
      LambdaLow/mapping_chunks/m54008_160308_002050.chunk4.alignmentset.xml: $
      map LambdaLow/subreads_chunks/m54008_160308_002050.chunk4.subreadset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/lambdaNEB/sequence/lambdaNEB.fasta
  
  build $
      LambdaLow/mapping_chunks/m54008_160308_002050.chunk5.alignmentset.xml: $
      map LambdaLow/subreads_chunks/m54008_160308_002050.chunk5.subreadset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/lambdaNEB/sequence/lambdaNEB.fasta
  
  build $
      LambdaLow/mapping_chunks/m54008_160308_002050.chunk6.alignmentset.xml: $
      map LambdaLow/subreads_chunks/m54008_160308_002050.chunk6.subreadset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/lambdaNEB/sequence/lambdaNEB.fasta
  
  build $
      LambdaLow/mapping_chunks/m54008_160308_002050.chunk7.alignmentset.xml: $
      map LambdaLow/subreads_chunks/m54008_160308_002050.chunk7.subreadset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/lambdaNEB/sequence/lambdaNEB.fasta
  
  build LambdaLow/mapping/all_movies.alignmentset.xml: $
      mergeDatasetsForCondition $
      LambdaLow/mapping_chunks/m54008_160308_002050.chunk0.alignmentset.xml $
      LambdaLow/mapping_chunks/m54008_160308_002050.chunk1.alignmentset.xml $
      LambdaLow/mapping_chunks/m54008_160308_002050.chunk2.alignmentset.xml $
      LambdaLow/mapping_chunks/m54008_160308_002050.chunk3.alignmentset.xml $
      LambdaLow/mapping_chunks/m54008_160308_002050.chunk4.alignmentset.xml $
      LambdaLow/mapping_chunks/m54008_160308_002050.chunk5.alignmentset.xml $
      LambdaLow/mapping_chunks/m54008_160308_002050.chunk6.alignmentset.xml $
      LambdaLow/mapping_chunks/m54008_160308_002050.chunk7.alignmentset.xml
  
  build EcoliLow/variant_calling/alignments_summary.gff: summarize_coverage $
      EcoliLow/mapping/all_movies.alignmentset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/ecoliK12_pbi_March2013/sequence/ecoliK12_pbi_March2013.fasta
  
  build EcoliLow/variant_calling/arrow/variants-5.gff: variantCalling $
      EcoliLow/mapping/all_movies.alignmentset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/ecoliK12_pbi_March2013/sequence/ecoliK12_pbi_March2013.fasta
    consensusFasta = EcoliLow/variant_calling/arrow/consensus-5.fasta
    modelSpec = 
    modelPath = 
    consensusFastq = EcoliLow/variant_calling/arrow/consensus-5.fastq
    coverageLimitArgument = -X5
  
  build EcoliLow/variant_calling/arrow/masked-variants-5.gff: maskVariantsGff $
      EcoliLow/variant_calling/arrow/variants-5.gff
    referenceMask = $
        /pbi/dept/consensus/bauhaus/genome-masks/ecoliK12_pbi_March2013-mask.gff
  
  build EcoliLow/variant_calling/arrow/variants-10.gff: variantCalling $
      EcoliLow/mapping/all_movies.alignmentset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/ecoliK12_pbi_March2013/sequence/ecoliK12_pbi_March2013.fasta
    consensusFasta = EcoliLow/variant_calling/arrow/consensus-10.fasta
    modelSpec = 
    modelPath = 
    consensusFastq = EcoliLow/variant_calling/arrow/consensus-10.fastq
    coverageLimitArgument = -X10
  
  build EcoliLow/variant_calling/arrow/masked-variants-10.gff: $
      maskVariantsGff EcoliLow/variant_calling/arrow/variants-10.gff
    referenceMask = $
        /pbi/dept/consensus/bauhaus/genome-masks/ecoliK12_pbi_March2013-mask.gff
  
  build EcoliLow/variant_calling/arrow/variants-15.gff: variantCalling $
      EcoliLow/mapping/all_movies.alignmentset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/ecoliK12_pbi_March2013/sequence/ecoliK12_pbi_March2013.fasta
    consensusFasta = EcoliLow/variant_calling/arrow/consensus-15.fasta
    modelSpec = 
    modelPath = 
    consensusFastq = EcoliLow/variant_calling/arrow/consensus-15.fastq
    coverageLimitArgument = -X15
  
  build EcoliLow/variant_calling/arrow/masked-variants-15.gff: $
      maskVariantsGff EcoliLow/variant_calling/arrow/variants-15.gff
    referenceMask = $
        /pbi/dept/consensus/bauhaus/genome-masks/ecoliK12_pbi_March2013-mask.gff
  
  build EcoliLow/variant_calling/arrow/variants-20.gff: variantCalling $
      EcoliLow/mapping/all_movies.alignmentset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/ecoliK12_pbi_March2013/sequence/ecoliK12_pbi_March2013.fasta
    consensusFasta = EcoliLow/variant_calling/arrow/consensus-20.fasta
    modelSpec = 
    modelPath = 
    consensusFastq = EcoliLow/variant_calling/arrow/consensus-20.fastq
    coverageLimitArgument = -X20
  
  build EcoliLow/variant_calling/arrow/masked-variants-20.gff: $
      maskVariantsGff EcoliLow/variant_calling/arrow/variants-20.gff
    referenceMask = $
        /pbi/dept/consensus/bauhaus/genome-masks/ecoliK12_pbi_March2013-mask.gff
  
  build EcoliLow/variant_calling/arrow/variants-30.gff: variantCalling $
      EcoliLow/mapping/all_movies.alignmentset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/ecoliK12_pbi_March2013/sequence/ecoliK12_pbi_March2013.fasta
    consensusFasta = EcoliLow/variant_calling/arrow/consensus-30.fasta
    modelSpec = 
    modelPath = 
    consensusFastq = EcoliLow/variant_calling/arrow/consensus-30.fastq
    coverageLimitArgument = -X30
  
  build EcoliLow/variant_calling/arrow/masked-variants-30.gff: $
      maskVariantsGff EcoliLow/variant_calling/arrow/variants-30.gff
    referenceMask = $
        /pbi/dept/consensus/bauhaus/genome-masks/ecoliK12_pbi_March2013-mask.gff
  
  build EcoliLow/variant_calling/arrow/variants-40.gff: variantCalling $
      EcoliLow/mapping/all_movies.alignmentset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/ecoliK12_pbi_March2013/sequence/ecoliK12_pbi_March2013.fasta
    consensusFasta = EcoliLow/variant_calling/arrow/consensus-40.fasta
    modelSpec = 
    modelPath = 
    consensusFastq = EcoliLow/variant_calling/arrow/consensus-40.fastq
    coverageLimitArgument = -X40
  
  build EcoliLow/variant_calling/arrow/masked-variants-40.gff: $
      maskVariantsGff EcoliLow/variant_calling/arrow/variants-40.gff
    referenceMask = $
        /pbi/dept/consensus/bauhaus/genome-masks/ecoliK12_pbi_March2013-mask.gff
  
  build EcoliLow/variant_calling/arrow/variants-50.gff: variantCalling $
      EcoliLow/mapping/all_movies.alignmentset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/ecoliK12_pbi_March2013/sequence/ecoliK12_pbi_March2013.fasta
    consensusFasta = EcoliLow/variant_calling/arrow/consensus-50.fasta
    modelSpec = 
    modelPath = 
    consensusFastq = EcoliLow/variant_calling/arrow/consensus-50.fastq
    coverageLimitArgument = -X50
  
  build EcoliLow/variant_calling/arrow/masked-variants-50.gff: $
      maskVariantsGff EcoliLow/variant_calling/arrow/variants-50.gff
    referenceMask = $
        /pbi/dept/consensus/bauhaus/genome-masks/ecoliK12_pbi_March2013-mask.gff
  
  build EcoliLow/variant_calling/arrow/variants-60.gff: variantCalling $
      EcoliLow/mapping/all_movies.alignmentset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/ecoliK12_pbi_March2013/sequence/ecoliK12_pbi_March2013.fasta
    consensusFasta = EcoliLow/variant_calling/arrow/consensus-60.fasta
    modelSpec = 
    modelPath = 
    consensusFastq = EcoliLow/variant_calling/arrow/consensus-60.fastq
    coverageLimitArgument = -X60
  
  build EcoliLow/variant_calling/arrow/masked-variants-60.gff: $
      maskVariantsGff EcoliLow/variant_calling/arrow/variants-60.gff
    referenceMask = $
        /pbi/dept/consensus/bauhaus/genome-masks/ecoliK12_pbi_March2013-mask.gff
  
  build EcoliLow/variant_calling/arrow/variants-80.gff: variantCalling $
      EcoliLow/mapping/all_movies.alignmentset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/ecoliK12_pbi_March2013/sequence/ecoliK12_pbi_March2013.fasta
    consensusFasta = EcoliLow/variant_calling/arrow/consensus-80.fasta
    modelSpec = 
    modelPath = 
    consensusFastq = EcoliLow/variant_calling/arrow/consensus-80.fastq
    coverageLimitArgument = -X80
  
  build EcoliLow/variant_calling/arrow/masked-variants-80.gff: $
      maskVariantsGff EcoliLow/variant_calling/arrow/variants-80.gff
    referenceMask = $
        /pbi/dept/consensus/bauhaus/genome-masks/ecoliK12_pbi_March2013-mask.gff
  
  build EcoliLow/variant_calling/arrow/variants-100.gff: variantCalling $
      EcoliLow/mapping/all_movies.alignmentset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/ecoliK12_pbi_March2013/sequence/ecoliK12_pbi_March2013.fasta
    consensusFasta = EcoliLow/variant_calling/arrow/consensus-100.fasta
    modelSpec = 
    modelPath = 
    consensusFastq = EcoliLow/variant_calling/arrow/consensus-100.fastq
    coverageLimitArgument = -X100
  
  build EcoliLow/variant_calling/arrow/masked-variants-100.gff: $
      maskVariantsGff EcoliLow/variant_calling/arrow/variants-100.gff
    referenceMask = $
        /pbi/dept/consensus/bauhaus/genome-masks/ecoliK12_pbi_March2013-mask.gff
  
  build LambdaLow/variant_calling/alignments_summary.gff: summarize_coverage $
      LambdaLow/mapping/all_movies.alignmentset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/lambdaNEB/sequence/lambdaNEB.fasta
  
  build LambdaLow/variant_calling/arrow/variants-5.gff: variantCalling $
      LambdaLow/mapping/all_movies.alignmentset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/lambdaNEB/sequence/lambdaNEB.fasta
    consensusFasta = LambdaLow/variant_calling/arrow/consensus-5.fasta
    modelSpec = 
    modelPath = 
    consensusFastq = LambdaLow/variant_calling/arrow/consensus-5.fastq
    coverageLimitArgument = -X5
  
  build LambdaLow/variant_calling/arrow/masked-variants-5.gff: $
      maskVariantsGff LambdaLow/variant_calling/arrow/variants-5.gff
    referenceMask = /pbi/dept/consensus/bauhaus/genome-masks/lambdaNEB-mask.gff
  
  build LambdaLow/variant_calling/arrow/variants-10.gff: variantCalling $
      LambdaLow/mapping/all_movies.alignmentset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/lambdaNEB/sequence/lambdaNEB.fasta
    consensusFasta = LambdaLow/variant_calling/arrow/consensus-10.fasta
    modelSpec = 
    modelPath = 
    consensusFastq = LambdaLow/variant_calling/arrow/consensus-10.fastq
    coverageLimitArgument = -X10
  
  build LambdaLow/variant_calling/arrow/masked-variants-10.gff: $
      maskVariantsGff LambdaLow/variant_calling/arrow/variants-10.gff
    referenceMask = /pbi/dept/consensus/bauhaus/genome-masks/lambdaNEB-mask.gff
  
  build LambdaLow/variant_calling/arrow/variants-15.gff: variantCalling $
      LambdaLow/mapping/all_movies.alignmentset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/lambdaNEB/sequence/lambdaNEB.fasta
    consensusFasta = LambdaLow/variant_calling/arrow/consensus-15.fasta
    modelSpec = 
    modelPath = 
    consensusFastq = LambdaLow/variant_calling/arrow/consensus-15.fastq
    coverageLimitArgument = -X15
  
  build LambdaLow/variant_calling/arrow/masked-variants-15.gff: $
      maskVariantsGff LambdaLow/variant_calling/arrow/variants-15.gff
    referenceMask = /pbi/dept/consensus/bauhaus/genome-masks/lambdaNEB-mask.gff
  
  build LambdaLow/variant_calling/arrow/variants-20.gff: variantCalling $
      LambdaLow/mapping/all_movies.alignmentset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/lambdaNEB/sequence/lambdaNEB.fasta
    consensusFasta = LambdaLow/variant_calling/arrow/consensus-20.fasta
    modelSpec = 
    modelPath = 
    consensusFastq = LambdaLow/variant_calling/arrow/consensus-20.fastq
    coverageLimitArgument = -X20
  
  build LambdaLow/variant_calling/arrow/masked-variants-20.gff: $
      maskVariantsGff LambdaLow/variant_calling/arrow/variants-20.gff
    referenceMask = /pbi/dept/consensus/bauhaus/genome-masks/lambdaNEB-mask.gff
  
  build LambdaLow/variant_calling/arrow/variants-30.gff: variantCalling $
      LambdaLow/mapping/all_movies.alignmentset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/lambdaNEB/sequence/lambdaNEB.fasta
    consensusFasta = LambdaLow/variant_calling/arrow/consensus-30.fasta
    modelSpec = 
    modelPath = 
    consensusFastq = LambdaLow/variant_calling/arrow/consensus-30.fastq
    coverageLimitArgument = -X30
  
  build LambdaLow/variant_calling/arrow/masked-variants-30.gff: $
      maskVariantsGff LambdaLow/variant_calling/arrow/variants-30.gff
    referenceMask = /pbi/dept/consensus/bauhaus/genome-masks/lambdaNEB-mask.gff
  
  build LambdaLow/variant_calling/arrow/variants-40.gff: variantCalling $
      LambdaLow/mapping/all_movies.alignmentset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/lambdaNEB/sequence/lambdaNEB.fasta
    consensusFasta = LambdaLow/variant_calling/arrow/consensus-40.fasta
    modelSpec = 
    modelPath = 
    consensusFastq = LambdaLow/variant_calling/arrow/consensus-40.fastq
    coverageLimitArgument = -X40
  
  build LambdaLow/variant_calling/arrow/masked-variants-40.gff: $
      maskVariantsGff LambdaLow/variant_calling/arrow/variants-40.gff
    referenceMask = /pbi/dept/consensus/bauhaus/genome-masks/lambdaNEB-mask.gff
  
  build LambdaLow/variant_calling/arrow/variants-50.gff: variantCalling $
      LambdaLow/mapping/all_movies.alignmentset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/lambdaNEB/sequence/lambdaNEB.fasta
    consensusFasta = LambdaLow/variant_calling/arrow/consensus-50.fasta
    modelSpec = 
    modelPath = 
    consensusFastq = LambdaLow/variant_calling/arrow/consensus-50.fastq
    coverageLimitArgument = -X50
  
  build LambdaLow/variant_calling/arrow/masked-variants-50.gff: $
      maskVariantsGff LambdaLow/variant_calling/arrow/variants-50.gff
    referenceMask = /pbi/dept/consensus/bauhaus/genome-masks/lambdaNEB-mask.gff
  
  build LambdaLow/variant_calling/arrow/variants-60.gff: variantCalling $
      LambdaLow/mapping/all_movies.alignmentset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/lambdaNEB/sequence/lambdaNEB.fasta
    consensusFasta = LambdaLow/variant_calling/arrow/consensus-60.fasta
    modelSpec = 
    modelPath = 
    consensusFastq = LambdaLow/variant_calling/arrow/consensus-60.fastq
    coverageLimitArgument = -X60
  
  build LambdaLow/variant_calling/arrow/masked-variants-60.gff: $
      maskVariantsGff LambdaLow/variant_calling/arrow/variants-60.gff
    referenceMask = /pbi/dept/consensus/bauhaus/genome-masks/lambdaNEB-mask.gff
  
  build LambdaLow/variant_calling/arrow/variants-80.gff: variantCalling $
      LambdaLow/mapping/all_movies.alignmentset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/lambdaNEB/sequence/lambdaNEB.fasta
    consensusFasta = LambdaLow/variant_calling/arrow/consensus-80.fasta
    modelSpec = 
    modelPath = 
    consensusFastq = LambdaLow/variant_calling/arrow/consensus-80.fastq
    coverageLimitArgument = -X80
  
  build LambdaLow/variant_calling/arrow/masked-variants-80.gff: $
      maskVariantsGff LambdaLow/variant_calling/arrow/variants-80.gff
    referenceMask = /pbi/dept/consensus/bauhaus/genome-masks/lambdaNEB-mask.gff
  
  build LambdaLow/variant_calling/arrow/variants-100.gff: variantCalling $
      LambdaLow/mapping/all_movies.alignmentset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/lambdaNEB/sequence/lambdaNEB.fasta
    consensusFasta = LambdaLow/variant_calling/arrow/consensus-100.fasta
    modelSpec = 
    modelPath = 
    consensusFastq = LambdaLow/variant_calling/arrow/consensus-100.fastq
    coverageLimitArgument = -X100
  
  build LambdaLow/variant_calling/arrow/masked-variants-100.gff: $
      maskVariantsGff LambdaLow/variant_calling/arrow/variants-100.gff
    referenceMask = /pbi/dept/consensus/bauhaus/genome-masks/lambdaNEB-mask.gff
  
  build LambdaHigh/variant_calling/alignments_summary.gff: summarize_coverage $
      LambdaHigh/mapping/all_movies.alignmentset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/lambdaNEB/sequence/lambdaNEB.fasta
  
  build LambdaHigh/variant_calling/arrow/variants-5.gff: variantCalling $
      LambdaHigh/mapping/all_movies.alignmentset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/lambdaNEB/sequence/lambdaNEB.fasta
    consensusFasta = LambdaHigh/variant_calling/arrow/consensus-5.fasta
    modelSpec = 
    modelPath = 
    consensusFastq = LambdaHigh/variant_calling/arrow/consensus-5.fastq
    coverageLimitArgument = -X5
  
  build LambdaHigh/variant_calling/arrow/masked-variants-5.gff: $
      maskVariantsGff LambdaHigh/variant_calling/arrow/variants-5.gff
    referenceMask = /pbi/dept/consensus/bauhaus/genome-masks/lambdaNEB-mask.gff
  
  build LambdaHigh/variant_calling/arrow/variants-10.gff: variantCalling $
      LambdaHigh/mapping/all_movies.alignmentset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/lambdaNEB/sequence/lambdaNEB.fasta
    consensusFasta = LambdaHigh/variant_calling/arrow/consensus-10.fasta
    modelSpec = 
    modelPath = 
    consensusFastq = LambdaHigh/variant_calling/arrow/consensus-10.fastq
    coverageLimitArgument = -X10
  
  build LambdaHigh/variant_calling/arrow/masked-variants-10.gff: $
      maskVariantsGff LambdaHigh/variant_calling/arrow/variants-10.gff
    referenceMask = /pbi/dept/consensus/bauhaus/genome-masks/lambdaNEB-mask.gff
  
  build LambdaHigh/variant_calling/arrow/variants-15.gff: variantCalling $
      LambdaHigh/mapping/all_movies.alignmentset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/lambdaNEB/sequence/lambdaNEB.fasta
    consensusFasta = LambdaHigh/variant_calling/arrow/consensus-15.fasta
    modelSpec = 
    modelPath = 
    consensusFastq = LambdaHigh/variant_calling/arrow/consensus-15.fastq
    coverageLimitArgument = -X15
  
  build LambdaHigh/variant_calling/arrow/masked-variants-15.gff: $
      maskVariantsGff LambdaHigh/variant_calling/arrow/variants-15.gff
    referenceMask = /pbi/dept/consensus/bauhaus/genome-masks/lambdaNEB-mask.gff
  
  build LambdaHigh/variant_calling/arrow/variants-20.gff: variantCalling $
      LambdaHigh/mapping/all_movies.alignmentset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/lambdaNEB/sequence/lambdaNEB.fasta
    consensusFasta = LambdaHigh/variant_calling/arrow/consensus-20.fasta
    modelSpec = 
    modelPath = 
    consensusFastq = LambdaHigh/variant_calling/arrow/consensus-20.fastq
    coverageLimitArgument = -X20
  
  build LambdaHigh/variant_calling/arrow/masked-variants-20.gff: $
      maskVariantsGff LambdaHigh/variant_calling/arrow/variants-20.gff
    referenceMask = /pbi/dept/consensus/bauhaus/genome-masks/lambdaNEB-mask.gff
  
  build LambdaHigh/variant_calling/arrow/variants-30.gff: variantCalling $
      LambdaHigh/mapping/all_movies.alignmentset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/lambdaNEB/sequence/lambdaNEB.fasta
    consensusFasta = LambdaHigh/variant_calling/arrow/consensus-30.fasta
    modelSpec = 
    modelPath = 
    consensusFastq = LambdaHigh/variant_calling/arrow/consensus-30.fastq
    coverageLimitArgument = -X30
  
  build LambdaHigh/variant_calling/arrow/masked-variants-30.gff: $
      maskVariantsGff LambdaHigh/variant_calling/arrow/variants-30.gff
    referenceMask = /pbi/dept/consensus/bauhaus/genome-masks/lambdaNEB-mask.gff
  
  build LambdaHigh/variant_calling/arrow/variants-40.gff: variantCalling $
      LambdaHigh/mapping/all_movies.alignmentset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/lambdaNEB/sequence/lambdaNEB.fasta
    consensusFasta = LambdaHigh/variant_calling/arrow/consensus-40.fasta
    modelSpec = 
    modelPath = 
    consensusFastq = LambdaHigh/variant_calling/arrow/consensus-40.fastq
    coverageLimitArgument = -X40
  
  build LambdaHigh/variant_calling/arrow/masked-variants-40.gff: $
      maskVariantsGff LambdaHigh/variant_calling/arrow/variants-40.gff
    referenceMask = /pbi/dept/consensus/bauhaus/genome-masks/lambdaNEB-mask.gff
  
  build LambdaHigh/variant_calling/arrow/variants-50.gff: variantCalling $
      LambdaHigh/mapping/all_movies.alignmentset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/lambdaNEB/sequence/lambdaNEB.fasta
    consensusFasta = LambdaHigh/variant_calling/arrow/consensus-50.fasta
    modelSpec = 
    modelPath = 
    consensusFastq = LambdaHigh/variant_calling/arrow/consensus-50.fastq
    coverageLimitArgument = -X50
  
  build LambdaHigh/variant_calling/arrow/masked-variants-50.gff: $
      maskVariantsGff LambdaHigh/variant_calling/arrow/variants-50.gff
    referenceMask = /pbi/dept/consensus/bauhaus/genome-masks/lambdaNEB-mask.gff
  
  build LambdaHigh/variant_calling/arrow/variants-60.gff: variantCalling $
      LambdaHigh/mapping/all_movies.alignmentset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/lambdaNEB/sequence/lambdaNEB.fasta
    consensusFasta = LambdaHigh/variant_calling/arrow/consensus-60.fasta
    modelSpec = 
    modelPath = 
    consensusFastq = LambdaHigh/variant_calling/arrow/consensus-60.fastq
    coverageLimitArgument = -X60
  
  build LambdaHigh/variant_calling/arrow/masked-variants-60.gff: $
      maskVariantsGff LambdaHigh/variant_calling/arrow/variants-60.gff
    referenceMask = /pbi/dept/consensus/bauhaus/genome-masks/lambdaNEB-mask.gff
  
  build LambdaHigh/variant_calling/arrow/variants-80.gff: variantCalling $
      LambdaHigh/mapping/all_movies.alignmentset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/lambdaNEB/sequence/lambdaNEB.fasta
    consensusFasta = LambdaHigh/variant_calling/arrow/consensus-80.fasta
    modelSpec = 
    modelPath = 
    consensusFastq = LambdaHigh/variant_calling/arrow/consensus-80.fastq
    coverageLimitArgument = -X80
  
  build LambdaHigh/variant_calling/arrow/masked-variants-80.gff: $
      maskVariantsGff LambdaHigh/variant_calling/arrow/variants-80.gff
    referenceMask = /pbi/dept/consensus/bauhaus/genome-masks/lambdaNEB-mask.gff
  
  build LambdaHigh/variant_calling/arrow/variants-100.gff: variantCalling $
      LambdaHigh/mapping/all_movies.alignmentset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/lambdaNEB/sequence/lambdaNEB.fasta
    consensusFasta = LambdaHigh/variant_calling/arrow/consensus-100.fasta
    modelSpec = 
    modelPath = 
    consensusFastq = LambdaHigh/variant_calling/arrow/consensus-100.fastq
    coverageLimitArgument = -X100
  
  build LambdaHigh/variant_calling/arrow/masked-variants-100.gff: $
      maskVariantsGff LambdaHigh/variant_calling/arrow/variants-100.gff
    referenceMask = /pbi/dept/consensus/bauhaus/genome-masks/lambdaNEB-mask.gff
  
  build EcoliHigh/variant_calling/alignments_summary.gff: summarize_coverage $
      EcoliHigh/mapping/all_movies.alignmentset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/ecoliK12_pbi_March2013/sequence/ecoliK12_pbi_March2013.fasta
  
  build EcoliHigh/variant_calling/arrow/variants-5.gff: variantCalling $
      EcoliHigh/mapping/all_movies.alignmentset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/ecoliK12_pbi_March2013/sequence/ecoliK12_pbi_March2013.fasta
    consensusFasta = EcoliHigh/variant_calling/arrow/consensus-5.fasta
    modelSpec = 
    modelPath = 
    consensusFastq = EcoliHigh/variant_calling/arrow/consensus-5.fastq
    coverageLimitArgument = -X5
  
  build EcoliHigh/variant_calling/arrow/masked-variants-5.gff: $
      maskVariantsGff EcoliHigh/variant_calling/arrow/variants-5.gff
    referenceMask = $
        /pbi/dept/consensus/bauhaus/genome-masks/ecoliK12_pbi_March2013-mask.gff
  
  build EcoliHigh/variant_calling/arrow/variants-10.gff: variantCalling $
      EcoliHigh/mapping/all_movies.alignmentset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/ecoliK12_pbi_March2013/sequence/ecoliK12_pbi_March2013.fasta
    consensusFasta = EcoliHigh/variant_calling/arrow/consensus-10.fasta
    modelSpec = 
    modelPath = 
    consensusFastq = EcoliHigh/variant_calling/arrow/consensus-10.fastq
    coverageLimitArgument = -X10
  
  build EcoliHigh/variant_calling/arrow/masked-variants-10.gff: $
      maskVariantsGff EcoliHigh/variant_calling/arrow/variants-10.gff
    referenceMask = $
        /pbi/dept/consensus/bauhaus/genome-masks/ecoliK12_pbi_March2013-mask.gff
  
  build EcoliHigh/variant_calling/arrow/variants-15.gff: variantCalling $
      EcoliHigh/mapping/all_movies.alignmentset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/ecoliK12_pbi_March2013/sequence/ecoliK12_pbi_March2013.fasta
    consensusFasta = EcoliHigh/variant_calling/arrow/consensus-15.fasta
    modelSpec = 
    modelPath = 
    consensusFastq = EcoliHigh/variant_calling/arrow/consensus-15.fastq
    coverageLimitArgument = -X15
  
  build EcoliHigh/variant_calling/arrow/masked-variants-15.gff: $
      maskVariantsGff EcoliHigh/variant_calling/arrow/variants-15.gff
    referenceMask = $
        /pbi/dept/consensus/bauhaus/genome-masks/ecoliK12_pbi_March2013-mask.gff
  
  build EcoliHigh/variant_calling/arrow/variants-20.gff: variantCalling $
      EcoliHigh/mapping/all_movies.alignmentset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/ecoliK12_pbi_March2013/sequence/ecoliK12_pbi_March2013.fasta
    consensusFasta = EcoliHigh/variant_calling/arrow/consensus-20.fasta
    modelSpec = 
    modelPath = 
    consensusFastq = EcoliHigh/variant_calling/arrow/consensus-20.fastq
    coverageLimitArgument = -X20
  
  build EcoliHigh/variant_calling/arrow/masked-variants-20.gff: $
      maskVariantsGff EcoliHigh/variant_calling/arrow/variants-20.gff
    referenceMask = $
        /pbi/dept/consensus/bauhaus/genome-masks/ecoliK12_pbi_March2013-mask.gff
  
  build EcoliHigh/variant_calling/arrow/variants-30.gff: variantCalling $
      EcoliHigh/mapping/all_movies.alignmentset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/ecoliK12_pbi_March2013/sequence/ecoliK12_pbi_March2013.fasta
    consensusFasta = EcoliHigh/variant_calling/arrow/consensus-30.fasta
    modelSpec = 
    modelPath = 
    consensusFastq = EcoliHigh/variant_calling/arrow/consensus-30.fastq
    coverageLimitArgument = -X30
  
  build EcoliHigh/variant_calling/arrow/masked-variants-30.gff: $
      maskVariantsGff EcoliHigh/variant_calling/arrow/variants-30.gff
    referenceMask = $
        /pbi/dept/consensus/bauhaus/genome-masks/ecoliK12_pbi_March2013-mask.gff
  
  build EcoliHigh/variant_calling/arrow/variants-40.gff: variantCalling $
      EcoliHigh/mapping/all_movies.alignmentset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/ecoliK12_pbi_March2013/sequence/ecoliK12_pbi_March2013.fasta
    consensusFasta = EcoliHigh/variant_calling/arrow/consensus-40.fasta
    modelSpec = 
    modelPath = 
    consensusFastq = EcoliHigh/variant_calling/arrow/consensus-40.fastq
    coverageLimitArgument = -X40
  
  build EcoliHigh/variant_calling/arrow/masked-variants-40.gff: $
      maskVariantsGff EcoliHigh/variant_calling/arrow/variants-40.gff
    referenceMask = $
        /pbi/dept/consensus/bauhaus/genome-masks/ecoliK12_pbi_March2013-mask.gff
  
  build EcoliHigh/variant_calling/arrow/variants-50.gff: variantCalling $
      EcoliHigh/mapping/all_movies.alignmentset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/ecoliK12_pbi_March2013/sequence/ecoliK12_pbi_March2013.fasta
    consensusFasta = EcoliHigh/variant_calling/arrow/consensus-50.fasta
    modelSpec = 
    modelPath = 
    consensusFastq = EcoliHigh/variant_calling/arrow/consensus-50.fastq
    coverageLimitArgument = -X50
  
  build EcoliHigh/variant_calling/arrow/masked-variants-50.gff: $
      maskVariantsGff EcoliHigh/variant_calling/arrow/variants-50.gff
    referenceMask = $
        /pbi/dept/consensus/bauhaus/genome-masks/ecoliK12_pbi_March2013-mask.gff
  
  build EcoliHigh/variant_calling/arrow/variants-60.gff: variantCalling $
      EcoliHigh/mapping/all_movies.alignmentset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/ecoliK12_pbi_March2013/sequence/ecoliK12_pbi_March2013.fasta
    consensusFasta = EcoliHigh/variant_calling/arrow/consensus-60.fasta
    modelSpec = 
    modelPath = 
    consensusFastq = EcoliHigh/variant_calling/arrow/consensus-60.fastq
    coverageLimitArgument = -X60
  
  build EcoliHigh/variant_calling/arrow/masked-variants-60.gff: $
      maskVariantsGff EcoliHigh/variant_calling/arrow/variants-60.gff
    referenceMask = $
        /pbi/dept/consensus/bauhaus/genome-masks/ecoliK12_pbi_March2013-mask.gff
  
  build EcoliHigh/variant_calling/arrow/variants-80.gff: variantCalling $
      EcoliHigh/mapping/all_movies.alignmentset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/ecoliK12_pbi_March2013/sequence/ecoliK12_pbi_March2013.fasta
    consensusFasta = EcoliHigh/variant_calling/arrow/consensus-80.fasta
    modelSpec = 
    modelPath = 
    consensusFastq = EcoliHigh/variant_calling/arrow/consensus-80.fastq
    coverageLimitArgument = -X80
  
  build EcoliHigh/variant_calling/arrow/masked-variants-80.gff: $
      maskVariantsGff EcoliHigh/variant_calling/arrow/variants-80.gff
    referenceMask = $
        /pbi/dept/consensus/bauhaus/genome-masks/ecoliK12_pbi_March2013-mask.gff
  
  build EcoliHigh/variant_calling/arrow/variants-100.gff: variantCalling $
      EcoliHigh/mapping/all_movies.alignmentset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/ecoliK12_pbi_March2013/sequence/ecoliK12_pbi_March2013.fasta
    consensusFasta = EcoliHigh/variant_calling/arrow/consensus-100.fasta
    modelSpec = 
    modelPath = 
    consensusFastq = EcoliHigh/variant_calling/arrow/consensus-100.fastq
    coverageLimitArgument = -X100
  
  build EcoliHigh/variant_calling/arrow/masked-variants-100.gff: $
      maskVariantsGff EcoliHigh/variant_calling/arrow/variants-100.gff
    referenceMask = $
        /pbi/dept/consensus/bauhaus/genome-masks/ecoliK12_pbi_March2013-mask.gff
  
  build coverage-titration.csv coverage-titration.pdf: $
      coverageTitrationSummaryAnalysis $
      EcoliLow/variant_calling/alignments_summary.gff $
      EcoliLow/variant_calling/arrow/masked-variants-5.gff $
      EcoliLow/variant_calling/arrow/masked-variants-10.gff $
      EcoliLow/variant_calling/arrow/masked-variants-15.gff $
      EcoliLow/variant_calling/arrow/masked-variants-20.gff $
      EcoliLow/variant_calling/arrow/masked-variants-30.gff $
      EcoliLow/variant_calling/arrow/masked-variants-40.gff $
      EcoliLow/variant_calling/arrow/masked-variants-50.gff $
      EcoliLow/variant_calling/arrow/masked-variants-60.gff $
      EcoliLow/variant_calling/arrow/masked-variants-80.gff $
      EcoliLow/variant_calling/arrow/masked-variants-100.gff $
      LambdaLow/variant_calling/alignments_summary.gff $
      LambdaLow/variant_calling/arrow/masked-variants-5.gff $
      LambdaLow/variant_calling/arrow/masked-variants-10.gff $
      LambdaLow/variant_calling/arrow/masked-variants-15.gff $
      LambdaLow/variant_calling/arrow/masked-variants-20.gff $
      LambdaLow/variant_calling/arrow/masked-variants-30.gff $
      LambdaLow/variant_calling/arrow/masked-variants-40.gff $
      LambdaLow/variant_calling/arrow/masked-variants-50.gff $
      LambdaLow/variant_calling/arrow/masked-variants-60.gff $
      LambdaLow/variant_calling/arrow/masked-variants-80.gff $
      LambdaLow/variant_calling/arrow/masked-variants-100.gff $
      EcoliHigh/variant_calling/alignments_summary.gff $
      EcoliHigh/variant_calling/arrow/masked-variants-5.gff $
      EcoliHigh/variant_calling/arrow/masked-variants-10.gff $
      EcoliHigh/variant_calling/arrow/masked-variants-15.gff $
      EcoliHigh/variant_calling/arrow/masked-variants-20.gff $
      EcoliHigh/variant_calling/arrow/masked-variants-30.gff $
      EcoliHigh/variant_calling/arrow/masked-variants-40.gff $
      EcoliHigh/variant_calling/arrow/masked-variants-50.gff $
      EcoliHigh/variant_calling/arrow/masked-variants-60.gff $
      EcoliHigh/variant_calling/arrow/masked-variants-80.gff $
      EcoliHigh/variant_calling/arrow/masked-variants-100.gff $
      LambdaHigh/variant_calling/alignments_summary.gff $
      LambdaHigh/variant_calling/arrow/masked-variants-5.gff $
      LambdaHigh/variant_calling/arrow/masked-variants-10.gff $
      LambdaHigh/variant_calling/arrow/masked-variants-15.gff $
      LambdaHigh/variant_calling/arrow/masked-variants-20.gff $
      LambdaHigh/variant_calling/arrow/masked-variants-30.gff $
      LambdaHigh/variant_calling/arrow/masked-variants-40.gff $
      LambdaHigh/variant_calling/arrow/masked-variants-50.gff $
      LambdaHigh/variant_calling/arrow/masked-variants-60.gff $
      LambdaHigh/variant_calling/arrow/masked-variants-80.gff $
      LambdaHigh/variant_calling/arrow/masked-variants-100.gff
  


