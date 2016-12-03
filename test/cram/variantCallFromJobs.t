
Do variant calling from pre-existing mapping jobs.

  $ BH_ROOT=$TESTDIR/../../

  $ bauhaus -o vc -m -t ${BH_ROOT}test/data/lambdaAndEcoliJobs.csv -w VariantCalling generate
  Validation and input resolution succeeded.
  Runnable workflow written to directory "vc"

  $ tree vc
  vc
  |-- build.ninja
  |-- condition-table.csv
  |-- log
  `-- run.sh
  
  1 directory, 3 files



  $ cat vc/build.ninja
  # Variables
  ncpus = 8
  scratchDir = /scratch
  grid = qsub -sync y -cwd -V -b y -e log -o log
  gridSMP = $grid -pe smp $ncpus
  
  # Rules
  rule mergeDatasetsForCondition
    command = $grid dataset merge $out $in
  
  rule variantCalling
    command = $gridSMP variantCaller $modelPath $modelSpec --algorithm=arrow $
        $coverageLimitArgument -x0 -q0 -j $ncpus --reportEffectiveCoverage $
        $in -r $reference -o $out -o $consensusFasta -o $consensusFastq
  
  
  # Build targets
  build Ecoli/mapping/all_movies.alignmentset.xml: mergeDatasetsForCondition $
      /pbi/dept/secondary/siv/smrtlink/smrtlink-beta/smrtsuite/userdata/jobs_root/004/004110/tasks/pbalign.tasks.consolidate_bam-0/final.alignmentset.alignmentset.xml $
      /pbi/dept/secondary/siv/smrtlink/smrtlink-beta/smrtsuite/userdata/jobs_root/004/004111/tasks/pbalign.tasks.consolidate_bam-0/final.alignmentset.alignmentset.xml
  
  build Lambda/mapping/all_movies.alignmentset.xml: mergeDatasetsForCondition $
      /pbi/dept/secondary/siv/smrtlink/smrtlink-beta/smrtsuite/userdata/jobs_root/004/004183/tasks/pbalign.tasks.consolidate_bam-0/final.alignmentset.alignmentset.xml $
      /pbi/dept/secondary/siv/smrtlink/smrtlink-beta/smrtsuite/userdata/jobs_root/004/004206/tasks/pbalign.tasks.consolidate_bam-0/final.alignmentset.alignmentset.xml
  
  build Ecoli/variant_calling/arrow/variants.gff: variantCalling $
      Ecoli/mapping/all_movies.alignmentset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/ecoliK12_pbi_March2013/sequence/ecoliK12_pbi_March2013.fasta
    consensusFasta = Ecoli/variant_calling/arrow/consensus.fasta
    modelSpec = 
    modelPath = 
    consensusFastq = Ecoli/variant_calling/arrow/consensus.fastq
    coverageLimitArgument = 
  
  build Lambda/variant_calling/arrow/variants.gff: variantCalling $
      Lambda/mapping/all_movies.alignmentset.xml
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/lambdaNEB/sequence/lambdaNEB.fasta
    consensusFasta = Lambda/variant_calling/arrow/consensus.fasta
    modelSpec = 
    modelPath = 
    consensusFastq = Lambda/variant_calling/arrow/consensus.fastq
    coverageLimitArgument = 
  
