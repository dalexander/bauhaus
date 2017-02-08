
  $ BH_ROOT=$TESTDIR/../../


First, try a couple cases where there's heterogeneity within a condition.

  $ bauhaus -m -w Mapping -t ${BH_ROOT}test/data/bad-cts/nonuniform-conditions-1.csv validate
  Condition table validation error: Conditions must be homogeneous--no variation allowed in variables/settings within a condition.  (Offending condition: "ControlChem"; offending column: "Genome")
  [1]


  $ bauhaus -m -w Mapping -t ${BH_ROOT}test/data/bad-cts/nonuniform-conditions-2.csv validate
  Condition table validation error: Conditions must be homogeneous--no variation allowed in variables/settings within a condition.  (Offending condition: "SparklyChem"; offending column: "p_Chemistry")
  [1]


Try a case where a genome can't be found:

  $ bauhaus -m -w Mapping -t ${BH_ROOT}test/data/bad-cts/unrecognized-genome.csv validate
  Condition table validation error: Reference not found: WillyWonka
  [1]

Try a case where an input can't be found:

  $ bauhaus -m -w Mapping -t ${BH_ROOT}test/data/bad-cts/nonexistent-data.csv validate
  Condition table validation error: Input data not found: 3150128-0002/WillyWonka
  [1]

Missing mandatory p_ variables for CoverageTitration workflow

  $ bauhaus -m -w CoverageTitration -t ${BH_ROOT}test/data/lambdaAndEcoli.csv validate
  Condition table validation error: For CoverageTitration, there must be at least one covariate ("p_" variable) in the condition table
  [1]

Missing genome mask (CoverageTitration workflow)

  $ bauhaus -m -w CoverageTitration -t ${BH_ROOT}test/data/bad-cts/unrecognized-genome-mask.csv validate
  Condition table validation error: Reference mask (required for CoverageTitration) not found for genome: plasmidbell_v1
  [1]
