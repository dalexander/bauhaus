
  $ BH_ROOT=$TESTDIR/../../

  $ bauhaus -o ccs -m -t ${BH_ROOT}test//data/lambdaAndEcoli.csv -w CCS --chunks 0 generate
  Validation and input resolution succeeded.
  Runnable workflow written to directory "ccs"

  $ cat ccs/build.ninja
  # Variables
  ncpus = 8
  scratchDir = /scratch
  grid = qsub -sync y -cwd -V -b y -e log -o log
  gridSMP = $grid -pe smp $ncpus
  
  # Rules
  rule copySubreadsDataset
    command = $grid dataset create $out $in
  
  rule ccs
    command = $gridSMP ccs --force --numThreads=$ncpus $
        --reportFile=$ccsDiagnostics $modelPath $modelSpec $in $outBam && $
        dataset create --type ConsensusReadSet $out $outBam
  
  
  # Build targets
  build Ecoli/subreads/m54011_160305_235923.subreadset.xml: $
      copySubreadsDataset $
      /pbi/collections/315/3150122/r54011_20160305_235615/1_A01/m54011_160305_235923.subreadset.xml
  
  build Ecoli/subreads/m54011_160306_050740.subreadset.xml: $
      copySubreadsDataset $
      /pbi/collections/315/3150122/r54011_20160305_235615/2_B01/m54011_160306_050740.subreadset.xml
  
  build Ecoli/ccs/m54011_160305_235923.consensusreadset.xml: ccs $
      Ecoli/subreads/m54011_160305_235923.subreadset.xml
    ccsDiagnostics = Ecoli/m54011_160305_235923.ccs-report.txt
    outBam = Ecoli/ccs/m54011_160305_235923.ccs.bam
    modelPath = 
    modelSpec = 
  
  build Ecoli/ccs/m54011_160306_050740.consensusreadset.xml: ccs $
      Ecoli/subreads/m54011_160306_050740.subreadset.xml
    ccsDiagnostics = Ecoli/m54011_160306_050740.ccs-report.txt
    outBam = Ecoli/ccs/m54011_160306_050740.ccs.bam
    modelPath = 
    modelSpec = 
  
  build Lambda/subreads/m54008_160308_002050.subreadset.xml: $
      copySubreadsDataset $
      /pbi/collections/315/3150128/r54008_20160308_001811/1_A01/m54008_160308_002050.subreadset.xml
  
  build Lambda/subreads/m54008_160308_053311.subreadset.xml: $
      copySubreadsDataset $
      /pbi/collections/315/3150128/r54008_20160308_001811/2_B01/m54008_160308_053311.subreadset.xml
  
  build Lambda/ccs/m54008_160308_002050.consensusreadset.xml: ccs $
      Lambda/subreads/m54008_160308_002050.subreadset.xml
    ccsDiagnostics = Lambda/m54008_160308_002050.ccs-report.txt
    outBam = Lambda/ccs/m54008_160308_002050.ccs.bam
    modelPath = 
    modelSpec = 
  
  build Lambda/ccs/m54008_160308_053311.consensusreadset.xml: ccs $
      Lambda/subreads/m54008_160308_053311.subreadset.xml
    ccsDiagnostics = Lambda/m54008_160308_053311.ccs-report.txt
    outBam = Lambda/ccs/m54008_160308_053311.ccs.bam
    modelPath = 
    modelSpec = 
  


  $ bauhaus -o chunkedCCS -m -t ${BH_ROOT}test/data/lambdaAndEcoli.csv -w CCS --chunks 2 generate
  Validation and input resolution succeeded.
  Runnable workflow written to directory "chunkedCCS"

  $ cat chunkedCCS/build.ninja
  # Variables
  ncpus = 8
  scratchDir = /scratch
  grid = qsub -sync y -cwd -V -b y -e log -o log
  gridSMP = $grid -pe smp $ncpus
  
  # Rules
  rule copySubreadsDataset
    command = $grid dataset create $out $in
  
  rule ccs
    command = $gridSMP ccs --force --numThreads=$ncpus $
        --reportFile=$ccsDiagnostics $modelPath $modelSpec $in $outBam && $
        dataset create --type ConsensusReadSet $out $outBam
  
  rule splitByZmw
    command = $grid dataset split --zmws --targetSize 1 --chunks 2 --outdir $
        $outdir $in
  
  rule mergeDatasetsForMovie
    command = $grid dataset merge $out $in
  
  rule consolidateDatasetsForMovie
    command = $grid dataset consolidate $in $outBam $out
  
  rule mergeDatasetsForCondition
    command = $grid dataset merge $out $in
  
  
  # Build targets
  build Ecoli/subreads/m54011_160305_235923.subreadset.xml: $
      copySubreadsDataset $
      /pbi/collections/315/3150122/r54011_20160305_235615/1_A01/m54011_160305_235923.subreadset.xml
  
  build Ecoli/subreads/m54011_160306_050740.subreadset.xml: $
      copySubreadsDataset $
      /pbi/collections/315/3150122/r54011_20160305_235615/2_B01/m54011_160306_050740.subreadset.xml
  
  build Ecoli/subreads_chunks/m54011_160305_235923.chunk0.subreadset.xml $
      Ecoli/subreads_chunks/m54011_160305_235923.chunk1.subreadset.xml: $
      splitByZmw Ecoli/subreads/m54011_160305_235923.subreadset.xml
    outdir = Ecoli/subreads_chunks
  
  build Ecoli/ccs_chunks/m54011_160305_235923.chunk0.consensusreadset.xml: $
      ccs Ecoli/subreads_chunks/m54011_160305_235923.chunk0.subreadset.xml
    ccsDiagnostics = Ecoli/ccs_chunks/m54011_160305_235923.chunk0.ccs-report.txt
    outBam = Ecoli/ccs_chunks/m54011_160305_235923.chunk0.ccs.bam
    modelPath = 
    modelSpec = 
  
  build Ecoli/ccs_chunks/m54011_160305_235923.chunk1.consensusreadset.xml: $
      ccs Ecoli/subreads_chunks/m54011_160305_235923.chunk1.subreadset.xml
    ccsDiagnostics = Ecoli/ccs_chunks/m54011_160305_235923.chunk1.ccs-report.txt
    outBam = Ecoli/ccs_chunks/m54011_160305_235923.chunk1.ccs.bam
    modelPath = 
    modelSpec = 
  
  build Ecoli/ccs/m54011_160305_235923_preconsolidate.consensusreadset.xml: $
      mergeDatasetsForMovie $
      Ecoli/ccs_chunks/m54011_160305_235923.chunk0.consensusreadset.xml $
      Ecoli/ccs_chunks/m54011_160305_235923.chunk1.consensusreadset.xml
  
  build Ecoli/ccs/m54011_160305_235923.consensusreadset.xml: $
      consolidateDatasetsForMovie $
      Ecoli/ccs/m54011_160305_235923_preconsolidate.consensusreadset.xml
    outBam = Ecoli/ccs/m54011_160305_235923.ccs.bam
  
  build Ecoli/subreads_chunks/m54011_160306_050740.chunk0.subreadset.xml $
      Ecoli/subreads_chunks/m54011_160306_050740.chunk1.subreadset.xml: $
      splitByZmw Ecoli/subreads/m54011_160306_050740.subreadset.xml
    outdir = Ecoli/subreads_chunks
  
  build Ecoli/ccs_chunks/m54011_160306_050740.chunk0.consensusreadset.xml: $
      ccs Ecoli/subreads_chunks/m54011_160306_050740.chunk0.subreadset.xml
    ccsDiagnostics = Ecoli/ccs_chunks/m54011_160306_050740.chunk0.ccs-report.txt
    outBam = Ecoli/ccs_chunks/m54011_160306_050740.chunk0.ccs.bam
    modelPath = 
    modelSpec = 
  
  build Ecoli/ccs_chunks/m54011_160306_050740.chunk1.consensusreadset.xml: $
      ccs Ecoli/subreads_chunks/m54011_160306_050740.chunk1.subreadset.xml
    ccsDiagnostics = Ecoli/ccs_chunks/m54011_160306_050740.chunk1.ccs-report.txt
    outBam = Ecoli/ccs_chunks/m54011_160306_050740.chunk1.ccs.bam
    modelPath = 
    modelSpec = 
  
  build Ecoli/ccs/m54011_160306_050740_preconsolidate.consensusreadset.xml: $
      mergeDatasetsForMovie $
      Ecoli/ccs_chunks/m54011_160306_050740.chunk0.consensusreadset.xml $
      Ecoli/ccs_chunks/m54011_160306_050740.chunk1.consensusreadset.xml
  
  build Ecoli/ccs/m54011_160306_050740.consensusreadset.xml: $
      consolidateDatasetsForMovie $
      Ecoli/ccs/m54011_160306_050740_preconsolidate.consensusreadset.xml
    outBam = Ecoli/ccs/m54011_160306_050740.ccs.bam
  
  build Ecoli/ccs/all_movies.consensusreadset.xml: mergeDatasetsForCondition $
      Ecoli/ccs/m54011_160305_235923.consensusreadset.xml $
      Ecoli/ccs/m54011_160306_050740.consensusreadset.xml
  
  build Lambda/subreads/m54008_160308_002050.subreadset.xml: $
      copySubreadsDataset $
      /pbi/collections/315/3150128/r54008_20160308_001811/1_A01/m54008_160308_002050.subreadset.xml
  
  build Lambda/subreads/m54008_160308_053311.subreadset.xml: $
      copySubreadsDataset $
      /pbi/collections/315/3150128/r54008_20160308_001811/2_B01/m54008_160308_053311.subreadset.xml
  
  build Lambda/subreads_chunks/m54008_160308_002050.chunk0.subreadset.xml $
      Lambda/subreads_chunks/m54008_160308_002050.chunk1.subreadset.xml: $
      splitByZmw Lambda/subreads/m54008_160308_002050.subreadset.xml
    outdir = Lambda/subreads_chunks
  
  build Lambda/ccs_chunks/m54008_160308_002050.chunk0.consensusreadset.xml: $
      ccs Lambda/subreads_chunks/m54008_160308_002050.chunk0.subreadset.xml
    ccsDiagnostics = $
        Lambda/ccs_chunks/m54008_160308_002050.chunk0.ccs-report.txt
    outBam = Lambda/ccs_chunks/m54008_160308_002050.chunk0.ccs.bam
    modelPath = 
    modelSpec = 
  
  build Lambda/ccs_chunks/m54008_160308_002050.chunk1.consensusreadset.xml: $
      ccs Lambda/subreads_chunks/m54008_160308_002050.chunk1.subreadset.xml
    ccsDiagnostics = $
        Lambda/ccs_chunks/m54008_160308_002050.chunk1.ccs-report.txt
    outBam = Lambda/ccs_chunks/m54008_160308_002050.chunk1.ccs.bam
    modelPath = 
    modelSpec = 
  
  build Lambda/ccs/m54008_160308_002050_preconsolidate.consensusreadset.xml: $
      mergeDatasetsForMovie $
      Lambda/ccs_chunks/m54008_160308_002050.chunk0.consensusreadset.xml $
      Lambda/ccs_chunks/m54008_160308_002050.chunk1.consensusreadset.xml
  
  build Lambda/ccs/m54008_160308_002050.consensusreadset.xml: $
      consolidateDatasetsForMovie $
      Lambda/ccs/m54008_160308_002050_preconsolidate.consensusreadset.xml
    outBam = Lambda/ccs/m54008_160308_002050.ccs.bam
  
  build Lambda/subreads_chunks/m54008_160308_053311.chunk0.subreadset.xml $
      Lambda/subreads_chunks/m54008_160308_053311.chunk1.subreadset.xml: $
      splitByZmw Lambda/subreads/m54008_160308_053311.subreadset.xml
    outdir = Lambda/subreads_chunks
  
  build Lambda/ccs_chunks/m54008_160308_053311.chunk0.consensusreadset.xml: $
      ccs Lambda/subreads_chunks/m54008_160308_053311.chunk0.subreadset.xml
    ccsDiagnostics = $
        Lambda/ccs_chunks/m54008_160308_053311.chunk0.ccs-report.txt
    outBam = Lambda/ccs_chunks/m54008_160308_053311.chunk0.ccs.bam
    modelPath = 
    modelSpec = 
  
  build Lambda/ccs_chunks/m54008_160308_053311.chunk1.consensusreadset.xml: $
      ccs Lambda/subreads_chunks/m54008_160308_053311.chunk1.subreadset.xml
    ccsDiagnostics = $
        Lambda/ccs_chunks/m54008_160308_053311.chunk1.ccs-report.txt
    outBam = Lambda/ccs_chunks/m54008_160308_053311.chunk1.ccs.bam
    modelPath = 
    modelSpec = 
  
  build Lambda/ccs/m54008_160308_053311_preconsolidate.consensusreadset.xml: $
      mergeDatasetsForMovie $
      Lambda/ccs_chunks/m54008_160308_053311.chunk0.consensusreadset.xml $
      Lambda/ccs_chunks/m54008_160308_053311.chunk1.consensusreadset.xml
  
  build Lambda/ccs/m54008_160308_053311.consensusreadset.xml: $
      consolidateDatasetsForMovie $
      Lambda/ccs/m54008_160308_053311_preconsolidate.consensusreadset.xml
    outBam = Lambda/ccs/m54008_160308_053311.ccs.bam
  
  build Lambda/ccs/all_movies.consensusreadset.xml: mergeDatasetsForCondition $
      Lambda/ccs/m54008_160308_002050.consensusreadset.xml $
      Lambda/ccs/m54008_160308_053311.consensusreadset.xml
  



  $ bauhaus -o ccsMap -m -t ${BH_ROOT}test/data/lambdaAndEcoli.csv -w CCSMapping generate
  Validation and input resolution succeeded.
  Runnable workflow written to directory "ccsMap"
