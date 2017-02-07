  $ BH_ROOT=$TESTDIR/../../

Mapping reports

  $ bauhaus -o mapReports -m -t ${BH_ROOT}test/data/lambdaAndEcoli.csv -w MappingReports generate
  Validation and input resolution succeeded.
  Runnable workflow written to directory "mapReports"


  $ tree mapReports
  mapReports
  |-- build.ninja
  |-- condition-table.csv
  |-- conditions.json
  |-- log
  |-- reports
  |   |-- PbiPlots
  |   |   `-- rtc-PbiPlots.json
  |   `-- PbiSampledPlots
  |       `-- rtc-PbiSampledPlots.json
  `-- run.sh
  
  4 directories, 6 files

