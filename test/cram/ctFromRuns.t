
Run a bigger workflow--coverage titration with reports.

  $ BH_ROOT=$TESTDIR/../../

  $ bauhaus -o ctFromRuns -m -t ${BH_ROOT}test/data/lambdaAndEcoli.csv -w CoverageTitrationReports generate
  Validation and input resolution succeeded.
  Runnable workflow written to directory "ctFromRuns"


  $ tree ctFromRuns
  ctFromRuns
  |-- R
  |   `-- coverageTitrationPlots.R
  |-- build.ninja
  |-- condition-table.csv
  |-- log
  `-- run.sh
  
  2 directories, 4 files



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
    command = Rscript --vanilla R/coverageTitrationPlots.R .
  
  
  # Build targets
  build Ecoli/subreads/m54011_160305_235923.subreadset.xml: $
      copySubreadsDataset $
      /pbi/collections/315/3150122/r54011_20160305_235615/1_A01/m54011_160305_235923.subreadset.xml
  
  build Ecoli/subreads/m54011_160306_050740.subreadset.xml: $
      copySubreadsDataset $
      /pbi/collections/315/3150122/r54011_20160305_235615/2_B01/m54011_160306_050740.subreadset.xml
  
  build Ecoli/subreads_chunks/m54011_160305_235923.chunk0.subreadset.xml $
      Ecoli/subreads_chunks/m54011_160305_235923.chunk1.subreadset.xml $
      Ecoli/subreads_chunks/m54011_160305_235923.chunk2.subreadset.xml $
      Ecoli/subreads_chunks/m54011_160305_235923.chunk3.subreadset.xml $
      Ecoli/subreads_chunks/m54011_160305_235923.chunk4.subreadset.xml $
      Ecoli/subreads_chunks/m54011_160305_235923.chunk5.subreadset.xml $
      Ecoli/subreads_chunks/m54011_160305_235923.chunk6.subreadset.xml $
      Ecoli/subreads_chunks/m54011_160305_235923.chunk7.subreadset.xml: $
      splitByZmw Ecoli/subreads/m54011_160305_235923.subreadset.xml
    outdir = Ecoli/subreads_chunks
  
  build Ecoli/mapping_chunks/m54011_160305_235923.chunk0.alignmentset.xml: $
      map Ecoli/subreads_chunks/m54011_160305_235923.chunk0.subreadset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/ecoliK12_pbi_March2013/sequence/ecoliK12_pbi_March2013.fasta
  
  build Ecoli/mapping_chunks/m54011_160305_235923.chunk1.alignmentset.xml: $
      map Ecoli/subreads_chunks/m54011_160305_235923.chunk1.subreadset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/ecoliK12_pbi_March2013/sequence/ecoliK12_pbi_March2013.fasta
  
  build Ecoli/mapping_chunks/m54011_160305_235923.chunk2.alignmentset.xml: $
      map Ecoli/subreads_chunks/m54011_160305_235923.chunk2.subreadset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/ecoliK12_pbi_March2013/sequence/ecoliK12_pbi_March2013.fasta
  
  build Ecoli/mapping_chunks/m54011_160305_235923.chunk3.alignmentset.xml: $
      map Ecoli/subreads_chunks/m54011_160305_235923.chunk3.subreadset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/ecoliK12_pbi_March2013/sequence/ecoliK12_pbi_March2013.fasta
  
  build Ecoli/mapping_chunks/m54011_160305_235923.chunk4.alignmentset.xml: $
      map Ecoli/subreads_chunks/m54011_160305_235923.chunk4.subreadset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/ecoliK12_pbi_March2013/sequence/ecoliK12_pbi_March2013.fasta
  
  build Ecoli/mapping_chunks/m54011_160305_235923.chunk5.alignmentset.xml: $
      map Ecoli/subreads_chunks/m54011_160305_235923.chunk5.subreadset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/ecoliK12_pbi_March2013/sequence/ecoliK12_pbi_March2013.fasta
  
  build Ecoli/mapping_chunks/m54011_160305_235923.chunk6.alignmentset.xml: $
      map Ecoli/subreads_chunks/m54011_160305_235923.chunk6.subreadset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/ecoliK12_pbi_March2013/sequence/ecoliK12_pbi_March2013.fasta
  
  build Ecoli/mapping_chunks/m54011_160305_235923.chunk7.alignmentset.xml: $
      map Ecoli/subreads_chunks/m54011_160305_235923.chunk7.subreadset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/ecoliK12_pbi_March2013/sequence/ecoliK12_pbi_March2013.fasta
  
  build Ecoli/subreads_chunks/m54011_160306_050740.chunk0.subreadset.xml $
      Ecoli/subreads_chunks/m54011_160306_050740.chunk1.subreadset.xml $
      Ecoli/subreads_chunks/m54011_160306_050740.chunk2.subreadset.xml $
      Ecoli/subreads_chunks/m54011_160306_050740.chunk3.subreadset.xml $
      Ecoli/subreads_chunks/m54011_160306_050740.chunk4.subreadset.xml $
      Ecoli/subreads_chunks/m54011_160306_050740.chunk5.subreadset.xml $
      Ecoli/subreads_chunks/m54011_160306_050740.chunk6.subreadset.xml $
      Ecoli/subreads_chunks/m54011_160306_050740.chunk7.subreadset.xml: $
      splitByZmw Ecoli/subreads/m54011_160306_050740.subreadset.xml
    outdir = Ecoli/subreads_chunks
  
  build Ecoli/mapping_chunks/m54011_160306_050740.chunk0.alignmentset.xml: $
      map Ecoli/subreads_chunks/m54011_160306_050740.chunk0.subreadset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/ecoliK12_pbi_March2013/sequence/ecoliK12_pbi_March2013.fasta
  
  build Ecoli/mapping_chunks/m54011_160306_050740.chunk1.alignmentset.xml: $
      map Ecoli/subreads_chunks/m54011_160306_050740.chunk1.subreadset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/ecoliK12_pbi_March2013/sequence/ecoliK12_pbi_March2013.fasta
  
  build Ecoli/mapping_chunks/m54011_160306_050740.chunk2.alignmentset.xml: $
      map Ecoli/subreads_chunks/m54011_160306_050740.chunk2.subreadset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/ecoliK12_pbi_March2013/sequence/ecoliK12_pbi_March2013.fasta
  
  build Ecoli/mapping_chunks/m54011_160306_050740.chunk3.alignmentset.xml: $
      map Ecoli/subreads_chunks/m54011_160306_050740.chunk3.subreadset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/ecoliK12_pbi_March2013/sequence/ecoliK12_pbi_March2013.fasta
  
  build Ecoli/mapping_chunks/m54011_160306_050740.chunk4.alignmentset.xml: $
      map Ecoli/subreads_chunks/m54011_160306_050740.chunk4.subreadset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/ecoliK12_pbi_March2013/sequence/ecoliK12_pbi_March2013.fasta
  
  build Ecoli/mapping_chunks/m54011_160306_050740.chunk5.alignmentset.xml: $
      map Ecoli/subreads_chunks/m54011_160306_050740.chunk5.subreadset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/ecoliK12_pbi_March2013/sequence/ecoliK12_pbi_March2013.fasta
  
  build Ecoli/mapping_chunks/m54011_160306_050740.chunk6.alignmentset.xml: $
      map Ecoli/subreads_chunks/m54011_160306_050740.chunk6.subreadset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/ecoliK12_pbi_March2013/sequence/ecoliK12_pbi_March2013.fasta
  
  build Ecoli/mapping_chunks/m54011_160306_050740.chunk7.alignmentset.xml: $
      map Ecoli/subreads_chunks/m54011_160306_050740.chunk7.subreadset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/ecoliK12_pbi_March2013/sequence/ecoliK12_pbi_March2013.fasta
  
  build Ecoli/mapping/all_movies.alignmentset.xml: mergeDatasetsForCondition $
      Ecoli/mapping_chunks/m54011_160305_235923.chunk0.alignmentset.xml $
      Ecoli/mapping_chunks/m54011_160305_235923.chunk1.alignmentset.xml $
      Ecoli/mapping_chunks/m54011_160305_235923.chunk2.alignmentset.xml $
      Ecoli/mapping_chunks/m54011_160305_235923.chunk3.alignmentset.xml $
      Ecoli/mapping_chunks/m54011_160305_235923.chunk4.alignmentset.xml $
      Ecoli/mapping_chunks/m54011_160305_235923.chunk5.alignmentset.xml $
      Ecoli/mapping_chunks/m54011_160305_235923.chunk6.alignmentset.xml $
      Ecoli/mapping_chunks/m54011_160305_235923.chunk7.alignmentset.xml $
      Ecoli/mapping_chunks/m54011_160306_050740.chunk0.alignmentset.xml $
      Ecoli/mapping_chunks/m54011_160306_050740.chunk1.alignmentset.xml $
      Ecoli/mapping_chunks/m54011_160306_050740.chunk2.alignmentset.xml $
      Ecoli/mapping_chunks/m54011_160306_050740.chunk3.alignmentset.xml $
      Ecoli/mapping_chunks/m54011_160306_050740.chunk4.alignmentset.xml $
      Ecoli/mapping_chunks/m54011_160306_050740.chunk5.alignmentset.xml $
      Ecoli/mapping_chunks/m54011_160306_050740.chunk6.alignmentset.xml $
      Ecoli/mapping_chunks/m54011_160306_050740.chunk7.alignmentset.xml
  
  build Lambda/subreads/m54008_160308_002050.subreadset.xml: $
      copySubreadsDataset $
      /pbi/collections/315/3150128/r54008_20160308_001811/1_A01/m54008_160308_002050.subreadset.xml
  
  build Lambda/subreads/m54008_160308_053311.subreadset.xml: $
      copySubreadsDataset $
      /pbi/collections/315/3150128/r54008_20160308_001811/2_B01/m54008_160308_053311.subreadset.xml
  
  build Lambda/subreads_chunks/m54008_160308_002050.chunk0.subreadset.xml $
      Lambda/subreads_chunks/m54008_160308_002050.chunk1.subreadset.xml $
      Lambda/subreads_chunks/m54008_160308_002050.chunk2.subreadset.xml $
      Lambda/subreads_chunks/m54008_160308_002050.chunk3.subreadset.xml $
      Lambda/subreads_chunks/m54008_160308_002050.chunk4.subreadset.xml $
      Lambda/subreads_chunks/m54008_160308_002050.chunk5.subreadset.xml $
      Lambda/subreads_chunks/m54008_160308_002050.chunk6.subreadset.xml $
      Lambda/subreads_chunks/m54008_160308_002050.chunk7.subreadset.xml: $
      splitByZmw Lambda/subreads/m54008_160308_002050.subreadset.xml
    outdir = Lambda/subreads_chunks
  
  build Lambda/mapping_chunks/m54008_160308_002050.chunk0.alignmentset.xml: $
      map Lambda/subreads_chunks/m54008_160308_002050.chunk0.subreadset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/lambdaNEB/sequence/lambdaNEB.fasta
  
  build Lambda/mapping_chunks/m54008_160308_002050.chunk1.alignmentset.xml: $
      map Lambda/subreads_chunks/m54008_160308_002050.chunk1.subreadset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/lambdaNEB/sequence/lambdaNEB.fasta
  
  build Lambda/mapping_chunks/m54008_160308_002050.chunk2.alignmentset.xml: $
      map Lambda/subreads_chunks/m54008_160308_002050.chunk2.subreadset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/lambdaNEB/sequence/lambdaNEB.fasta
  
  build Lambda/mapping_chunks/m54008_160308_002050.chunk3.alignmentset.xml: $
      map Lambda/subreads_chunks/m54008_160308_002050.chunk3.subreadset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/lambdaNEB/sequence/lambdaNEB.fasta
  
  build Lambda/mapping_chunks/m54008_160308_002050.chunk4.alignmentset.xml: $
      map Lambda/subreads_chunks/m54008_160308_002050.chunk4.subreadset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/lambdaNEB/sequence/lambdaNEB.fasta
  
  build Lambda/mapping_chunks/m54008_160308_002050.chunk5.alignmentset.xml: $
      map Lambda/subreads_chunks/m54008_160308_002050.chunk5.subreadset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/lambdaNEB/sequence/lambdaNEB.fasta
  
  build Lambda/mapping_chunks/m54008_160308_002050.chunk6.alignmentset.xml: $
      map Lambda/subreads_chunks/m54008_160308_002050.chunk6.subreadset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/lambdaNEB/sequence/lambdaNEB.fasta
  
  build Lambda/mapping_chunks/m54008_160308_002050.chunk7.alignmentset.xml: $
      map Lambda/subreads_chunks/m54008_160308_002050.chunk7.subreadset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/lambdaNEB/sequence/lambdaNEB.fasta
  
  build Lambda/subreads_chunks/m54008_160308_053311.chunk0.subreadset.xml $
      Lambda/subreads_chunks/m54008_160308_053311.chunk1.subreadset.xml $
      Lambda/subreads_chunks/m54008_160308_053311.chunk2.subreadset.xml $
      Lambda/subreads_chunks/m54008_160308_053311.chunk3.subreadset.xml $
      Lambda/subreads_chunks/m54008_160308_053311.chunk4.subreadset.xml $
      Lambda/subreads_chunks/m54008_160308_053311.chunk5.subreadset.xml $
      Lambda/subreads_chunks/m54008_160308_053311.chunk6.subreadset.xml $
      Lambda/subreads_chunks/m54008_160308_053311.chunk7.subreadset.xml: $
      splitByZmw Lambda/subreads/m54008_160308_053311.subreadset.xml
    outdir = Lambda/subreads_chunks
  
  build Lambda/mapping_chunks/m54008_160308_053311.chunk0.alignmentset.xml: $
      map Lambda/subreads_chunks/m54008_160308_053311.chunk0.subreadset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/lambdaNEB/sequence/lambdaNEB.fasta
  
  build Lambda/mapping_chunks/m54008_160308_053311.chunk1.alignmentset.xml: $
      map Lambda/subreads_chunks/m54008_160308_053311.chunk1.subreadset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/lambdaNEB/sequence/lambdaNEB.fasta
  
  build Lambda/mapping_chunks/m54008_160308_053311.chunk2.alignmentset.xml: $
      map Lambda/subreads_chunks/m54008_160308_053311.chunk2.subreadset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/lambdaNEB/sequence/lambdaNEB.fasta
  
  build Lambda/mapping_chunks/m54008_160308_053311.chunk3.alignmentset.xml: $
      map Lambda/subreads_chunks/m54008_160308_053311.chunk3.subreadset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/lambdaNEB/sequence/lambdaNEB.fasta
  
  build Lambda/mapping_chunks/m54008_160308_053311.chunk4.alignmentset.xml: $
      map Lambda/subreads_chunks/m54008_160308_053311.chunk4.subreadset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/lambdaNEB/sequence/lambdaNEB.fasta
  
  build Lambda/mapping_chunks/m54008_160308_053311.chunk5.alignmentset.xml: $
      map Lambda/subreads_chunks/m54008_160308_053311.chunk5.subreadset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/lambdaNEB/sequence/lambdaNEB.fasta
  
  build Lambda/mapping_chunks/m54008_160308_053311.chunk6.alignmentset.xml: $
      map Lambda/subreads_chunks/m54008_160308_053311.chunk6.subreadset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/lambdaNEB/sequence/lambdaNEB.fasta
  
  build Lambda/mapping_chunks/m54008_160308_053311.chunk7.alignmentset.xml: $
      map Lambda/subreads_chunks/m54008_160308_053311.chunk7.subreadset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/lambdaNEB/sequence/lambdaNEB.fasta
  
  build Lambda/mapping/all_movies.alignmentset.xml: mergeDatasetsForCondition $
      Lambda/mapping_chunks/m54008_160308_002050.chunk0.alignmentset.xml $
      Lambda/mapping_chunks/m54008_160308_002050.chunk1.alignmentset.xml $
      Lambda/mapping_chunks/m54008_160308_002050.chunk2.alignmentset.xml $
      Lambda/mapping_chunks/m54008_160308_002050.chunk3.alignmentset.xml $
      Lambda/mapping_chunks/m54008_160308_002050.chunk4.alignmentset.xml $
      Lambda/mapping_chunks/m54008_160308_002050.chunk5.alignmentset.xml $
      Lambda/mapping_chunks/m54008_160308_002050.chunk6.alignmentset.xml $
      Lambda/mapping_chunks/m54008_160308_002050.chunk7.alignmentset.xml $
      Lambda/mapping_chunks/m54008_160308_053311.chunk0.alignmentset.xml $
      Lambda/mapping_chunks/m54008_160308_053311.chunk1.alignmentset.xml $
      Lambda/mapping_chunks/m54008_160308_053311.chunk2.alignmentset.xml $
      Lambda/mapping_chunks/m54008_160308_053311.chunk3.alignmentset.xml $
      Lambda/mapping_chunks/m54008_160308_053311.chunk4.alignmentset.xml $
      Lambda/mapping_chunks/m54008_160308_053311.chunk5.alignmentset.xml $
      Lambda/mapping_chunks/m54008_160308_053311.chunk6.alignmentset.xml $
      Lambda/mapping_chunks/m54008_160308_053311.chunk7.alignmentset.xml
  
  build Ecoli/variant_calling/alignments_summary.gff: summarize_coverage $
      Ecoli/mapping/all_movies.alignmentset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/ecoliK12_pbi_March2013/sequence/ecoliK12_pbi_March2013.fasta
  
  build Ecoli/variant_calling/arrow/variants-5.gff: variantCalling $
      Ecoli/mapping/all_movies.alignmentset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/ecoliK12_pbi_March2013/sequence/ecoliK12_pbi_March2013.fasta
    consensusFasta = Ecoli/variant_calling/arrow/consensus-5.fasta
    modelSpec = 
    modelPath = 
    consensusFastq = Ecoli/variant_calling/arrow/consensus-5.fastq
    coverageLimitArgument = -X5
  
  build Ecoli/variant_calling/arrow/masked-variants-5.gff: maskVariantsGff $
      Ecoli/variant_calling/arrow/variants-5.gff
    referenceMask = $
        /mnt/secondary/Share/VariantCalling/Quiver/GenomeMasks/ecoliK12_pbi_March2013-mask.gff
  
  build Ecoli/variant_calling/arrow/variants-10.gff: variantCalling $
      Ecoli/mapping/all_movies.alignmentset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/ecoliK12_pbi_March2013/sequence/ecoliK12_pbi_March2013.fasta
    consensusFasta = Ecoli/variant_calling/arrow/consensus-10.fasta
    modelSpec = 
    modelPath = 
    consensusFastq = Ecoli/variant_calling/arrow/consensus-10.fastq
    coverageLimitArgument = -X10
  
  build Ecoli/variant_calling/arrow/masked-variants-10.gff: maskVariantsGff $
      Ecoli/variant_calling/arrow/variants-10.gff
    referenceMask = $
        /mnt/secondary/Share/VariantCalling/Quiver/GenomeMasks/ecoliK12_pbi_March2013-mask.gff
  
  build Ecoli/variant_calling/arrow/variants-15.gff: variantCalling $
      Ecoli/mapping/all_movies.alignmentset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/ecoliK12_pbi_March2013/sequence/ecoliK12_pbi_March2013.fasta
    consensusFasta = Ecoli/variant_calling/arrow/consensus-15.fasta
    modelSpec = 
    modelPath = 
    consensusFastq = Ecoli/variant_calling/arrow/consensus-15.fastq
    coverageLimitArgument = -X15
  
  build Ecoli/variant_calling/arrow/masked-variants-15.gff: maskVariantsGff $
      Ecoli/variant_calling/arrow/variants-15.gff
    referenceMask = $
        /mnt/secondary/Share/VariantCalling/Quiver/GenomeMasks/ecoliK12_pbi_March2013-mask.gff
  
  build Ecoli/variant_calling/arrow/variants-20.gff: variantCalling $
      Ecoli/mapping/all_movies.alignmentset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/ecoliK12_pbi_March2013/sequence/ecoliK12_pbi_March2013.fasta
    consensusFasta = Ecoli/variant_calling/arrow/consensus-20.fasta
    modelSpec = 
    modelPath = 
    consensusFastq = Ecoli/variant_calling/arrow/consensus-20.fastq
    coverageLimitArgument = -X20
  
  build Ecoli/variant_calling/arrow/masked-variants-20.gff: maskVariantsGff $
      Ecoli/variant_calling/arrow/variants-20.gff
    referenceMask = $
        /mnt/secondary/Share/VariantCalling/Quiver/GenomeMasks/ecoliK12_pbi_March2013-mask.gff
  
  build Ecoli/variant_calling/arrow/variants-30.gff: variantCalling $
      Ecoli/mapping/all_movies.alignmentset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/ecoliK12_pbi_March2013/sequence/ecoliK12_pbi_March2013.fasta
    consensusFasta = Ecoli/variant_calling/arrow/consensus-30.fasta
    modelSpec = 
    modelPath = 
    consensusFastq = Ecoli/variant_calling/arrow/consensus-30.fastq
    coverageLimitArgument = -X30
  
  build Ecoli/variant_calling/arrow/masked-variants-30.gff: maskVariantsGff $
      Ecoli/variant_calling/arrow/variants-30.gff
    referenceMask = $
        /mnt/secondary/Share/VariantCalling/Quiver/GenomeMasks/ecoliK12_pbi_March2013-mask.gff
  
  build Ecoli/variant_calling/arrow/variants-40.gff: variantCalling $
      Ecoli/mapping/all_movies.alignmentset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/ecoliK12_pbi_March2013/sequence/ecoliK12_pbi_March2013.fasta
    consensusFasta = Ecoli/variant_calling/arrow/consensus-40.fasta
    modelSpec = 
    modelPath = 
    consensusFastq = Ecoli/variant_calling/arrow/consensus-40.fastq
    coverageLimitArgument = -X40
  
  build Ecoli/variant_calling/arrow/masked-variants-40.gff: maskVariantsGff $
      Ecoli/variant_calling/arrow/variants-40.gff
    referenceMask = $
        /mnt/secondary/Share/VariantCalling/Quiver/GenomeMasks/ecoliK12_pbi_March2013-mask.gff
  
  build Ecoli/variant_calling/arrow/variants-50.gff: variantCalling $
      Ecoli/mapping/all_movies.alignmentset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/ecoliK12_pbi_March2013/sequence/ecoliK12_pbi_March2013.fasta
    consensusFasta = Ecoli/variant_calling/arrow/consensus-50.fasta
    modelSpec = 
    modelPath = 
    consensusFastq = Ecoli/variant_calling/arrow/consensus-50.fastq
    coverageLimitArgument = -X50
  
  build Ecoli/variant_calling/arrow/masked-variants-50.gff: maskVariantsGff $
      Ecoli/variant_calling/arrow/variants-50.gff
    referenceMask = $
        /mnt/secondary/Share/VariantCalling/Quiver/GenomeMasks/ecoliK12_pbi_March2013-mask.gff
  
  build Ecoli/variant_calling/arrow/variants-60.gff: variantCalling $
      Ecoli/mapping/all_movies.alignmentset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/ecoliK12_pbi_March2013/sequence/ecoliK12_pbi_March2013.fasta
    consensusFasta = Ecoli/variant_calling/arrow/consensus-60.fasta
    modelSpec = 
    modelPath = 
    consensusFastq = Ecoli/variant_calling/arrow/consensus-60.fastq
    coverageLimitArgument = -X60
  
  build Ecoli/variant_calling/arrow/masked-variants-60.gff: maskVariantsGff $
      Ecoli/variant_calling/arrow/variants-60.gff
    referenceMask = $
        /mnt/secondary/Share/VariantCalling/Quiver/GenomeMasks/ecoliK12_pbi_March2013-mask.gff
  
  build Ecoli/variant_calling/arrow/variants-80.gff: variantCalling $
      Ecoli/mapping/all_movies.alignmentset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/ecoliK12_pbi_March2013/sequence/ecoliK12_pbi_March2013.fasta
    consensusFasta = Ecoli/variant_calling/arrow/consensus-80.fasta
    modelSpec = 
    modelPath = 
    consensusFastq = Ecoli/variant_calling/arrow/consensus-80.fastq
    coverageLimitArgument = -X80
  
  build Ecoli/variant_calling/arrow/masked-variants-80.gff: maskVariantsGff $
      Ecoli/variant_calling/arrow/variants-80.gff
    referenceMask = $
        /mnt/secondary/Share/VariantCalling/Quiver/GenomeMasks/ecoliK12_pbi_March2013-mask.gff
  
  build Ecoli/variant_calling/arrow/variants-100.gff: variantCalling $
      Ecoli/mapping/all_movies.alignmentset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/ecoliK12_pbi_March2013/sequence/ecoliK12_pbi_March2013.fasta
    consensusFasta = Ecoli/variant_calling/arrow/consensus-100.fasta
    modelSpec = 
    modelPath = 
    consensusFastq = Ecoli/variant_calling/arrow/consensus-100.fastq
    coverageLimitArgument = -X100
  
  build Ecoli/variant_calling/arrow/masked-variants-100.gff: maskVariantsGff $
      Ecoli/variant_calling/arrow/variants-100.gff
    referenceMask = $
        /mnt/secondary/Share/VariantCalling/Quiver/GenomeMasks/ecoliK12_pbi_March2013-mask.gff
  
  build Lambda/variant_calling/alignments_summary.gff: summarize_coverage $
      Lambda/mapping/all_movies.alignmentset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/lambdaNEB/sequence/lambdaNEB.fasta
  
  build Lambda/variant_calling/arrow/variants-5.gff: variantCalling $
      Lambda/mapping/all_movies.alignmentset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/lambdaNEB/sequence/lambdaNEB.fasta
    consensusFasta = Lambda/variant_calling/arrow/consensus-5.fasta
    modelSpec = 
    modelPath = 
    consensusFastq = Lambda/variant_calling/arrow/consensus-5.fastq
    coverageLimitArgument = -X5
  
  build Lambda/variant_calling/arrow/masked-variants-5.gff: maskVariantsGff $
      Lambda/variant_calling/arrow/variants-5.gff
    referenceMask = $
        /mnt/secondary/Share/VariantCalling/Quiver/GenomeMasks/lambdaNEB-mask.gff
  
  build Lambda/variant_calling/arrow/variants-10.gff: variantCalling $
      Lambda/mapping/all_movies.alignmentset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/lambdaNEB/sequence/lambdaNEB.fasta
    consensusFasta = Lambda/variant_calling/arrow/consensus-10.fasta
    modelSpec = 
    modelPath = 
    consensusFastq = Lambda/variant_calling/arrow/consensus-10.fastq
    coverageLimitArgument = -X10
  
  build Lambda/variant_calling/arrow/masked-variants-10.gff: maskVariantsGff $
      Lambda/variant_calling/arrow/variants-10.gff
    referenceMask = $
        /mnt/secondary/Share/VariantCalling/Quiver/GenomeMasks/lambdaNEB-mask.gff
  
  build Lambda/variant_calling/arrow/variants-15.gff: variantCalling $
      Lambda/mapping/all_movies.alignmentset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/lambdaNEB/sequence/lambdaNEB.fasta
    consensusFasta = Lambda/variant_calling/arrow/consensus-15.fasta
    modelSpec = 
    modelPath = 
    consensusFastq = Lambda/variant_calling/arrow/consensus-15.fastq
    coverageLimitArgument = -X15
  
  build Lambda/variant_calling/arrow/masked-variants-15.gff: maskVariantsGff $
      Lambda/variant_calling/arrow/variants-15.gff
    referenceMask = $
        /mnt/secondary/Share/VariantCalling/Quiver/GenomeMasks/lambdaNEB-mask.gff
  
  build Lambda/variant_calling/arrow/variants-20.gff: variantCalling $
      Lambda/mapping/all_movies.alignmentset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/lambdaNEB/sequence/lambdaNEB.fasta
    consensusFasta = Lambda/variant_calling/arrow/consensus-20.fasta
    modelSpec = 
    modelPath = 
    consensusFastq = Lambda/variant_calling/arrow/consensus-20.fastq
    coverageLimitArgument = -X20
  
  build Lambda/variant_calling/arrow/masked-variants-20.gff: maskVariantsGff $
      Lambda/variant_calling/arrow/variants-20.gff
    referenceMask = $
        /mnt/secondary/Share/VariantCalling/Quiver/GenomeMasks/lambdaNEB-mask.gff
  
  build Lambda/variant_calling/arrow/variants-30.gff: variantCalling $
      Lambda/mapping/all_movies.alignmentset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/lambdaNEB/sequence/lambdaNEB.fasta
    consensusFasta = Lambda/variant_calling/arrow/consensus-30.fasta
    modelSpec = 
    modelPath = 
    consensusFastq = Lambda/variant_calling/arrow/consensus-30.fastq
    coverageLimitArgument = -X30
  
  build Lambda/variant_calling/arrow/masked-variants-30.gff: maskVariantsGff $
      Lambda/variant_calling/arrow/variants-30.gff
    referenceMask = $
        /mnt/secondary/Share/VariantCalling/Quiver/GenomeMasks/lambdaNEB-mask.gff
  
  build Lambda/variant_calling/arrow/variants-40.gff: variantCalling $
      Lambda/mapping/all_movies.alignmentset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/lambdaNEB/sequence/lambdaNEB.fasta
    consensusFasta = Lambda/variant_calling/arrow/consensus-40.fasta
    modelSpec = 
    modelPath = 
    consensusFastq = Lambda/variant_calling/arrow/consensus-40.fastq
    coverageLimitArgument = -X40
  
  build Lambda/variant_calling/arrow/masked-variants-40.gff: maskVariantsGff $
      Lambda/variant_calling/arrow/variants-40.gff
    referenceMask = $
        /mnt/secondary/Share/VariantCalling/Quiver/GenomeMasks/lambdaNEB-mask.gff
  
  build Lambda/variant_calling/arrow/variants-50.gff: variantCalling $
      Lambda/mapping/all_movies.alignmentset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/lambdaNEB/sequence/lambdaNEB.fasta
    consensusFasta = Lambda/variant_calling/arrow/consensus-50.fasta
    modelSpec = 
    modelPath = 
    consensusFastq = Lambda/variant_calling/arrow/consensus-50.fastq
    coverageLimitArgument = -X50
  
  build Lambda/variant_calling/arrow/masked-variants-50.gff: maskVariantsGff $
      Lambda/variant_calling/arrow/variants-50.gff
    referenceMask = $
        /mnt/secondary/Share/VariantCalling/Quiver/GenomeMasks/lambdaNEB-mask.gff
  
  build Lambda/variant_calling/arrow/variants-60.gff: variantCalling $
      Lambda/mapping/all_movies.alignmentset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/lambdaNEB/sequence/lambdaNEB.fasta
    consensusFasta = Lambda/variant_calling/arrow/consensus-60.fasta
    modelSpec = 
    modelPath = 
    consensusFastq = Lambda/variant_calling/arrow/consensus-60.fastq
    coverageLimitArgument = -X60
  
  build Lambda/variant_calling/arrow/masked-variants-60.gff: maskVariantsGff $
      Lambda/variant_calling/arrow/variants-60.gff
    referenceMask = $
        /mnt/secondary/Share/VariantCalling/Quiver/GenomeMasks/lambdaNEB-mask.gff
  
  build Lambda/variant_calling/arrow/variants-80.gff: variantCalling $
      Lambda/mapping/all_movies.alignmentset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/lambdaNEB/sequence/lambdaNEB.fasta
    consensusFasta = Lambda/variant_calling/arrow/consensus-80.fasta
    modelSpec = 
    modelPath = 
    consensusFastq = Lambda/variant_calling/arrow/consensus-80.fastq
    coverageLimitArgument = -X80
  
  build Lambda/variant_calling/arrow/masked-variants-80.gff: maskVariantsGff $
      Lambda/variant_calling/arrow/variants-80.gff
    referenceMask = $
        /mnt/secondary/Share/VariantCalling/Quiver/GenomeMasks/lambdaNEB-mask.gff
  
  build Lambda/variant_calling/arrow/variants-100.gff: variantCalling $
      Lambda/mapping/all_movies.alignmentset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/lambdaNEB/sequence/lambdaNEB.fasta
    consensusFasta = Lambda/variant_calling/arrow/consensus-100.fasta
    modelSpec = 
    modelPath = 
    consensusFastq = Lambda/variant_calling/arrow/consensus-100.fastq
    coverageLimitArgument = -X100
  
  build Lambda/variant_calling/arrow/masked-variants-100.gff: maskVariantsGff $
      Lambda/variant_calling/arrow/variants-100.gff
    referenceMask = $
        /mnt/secondary/Share/VariantCalling/Quiver/GenomeMasks/lambdaNEB-mask.gff
  
  build coverage-titration.csv coverage-titration.pdf: $
      coverageTitrationSummaryAnalysis $
      Ecoli/variant_calling/alignments_summary.gff $
      Ecoli/variant_calling/arrow/masked-variants-5.gff $
      Ecoli/variant_calling/arrow/masked-variants-10.gff $
      Ecoli/variant_calling/arrow/masked-variants-15.gff $
      Ecoli/variant_calling/arrow/masked-variants-20.gff $
      Ecoli/variant_calling/arrow/masked-variants-30.gff $
      Ecoli/variant_calling/arrow/masked-variants-40.gff $
      Ecoli/variant_calling/arrow/masked-variants-50.gff $
      Ecoli/variant_calling/arrow/masked-variants-60.gff $
      Ecoli/variant_calling/arrow/masked-variants-80.gff $
      Ecoli/variant_calling/arrow/masked-variants-100.gff $
      Lambda/variant_calling/alignments_summary.gff $
      Lambda/variant_calling/arrow/masked-variants-5.gff $
      Lambda/variant_calling/arrow/masked-variants-10.gff $
      Lambda/variant_calling/arrow/masked-variants-15.gff $
      Lambda/variant_calling/arrow/masked-variants-20.gff $
      Lambda/variant_calling/arrow/masked-variants-30.gff $
      Lambda/variant_calling/arrow/masked-variants-40.gff $
      Lambda/variant_calling/arrow/masked-variants-50.gff $
      Lambda/variant_calling/arrow/masked-variants-60.gff $
      Lambda/variant_calling/arrow/masked-variants-80.gff $
      Lambda/variant_calling/arrow/masked-variants-100.gff
  


