
## `bauhaus` Tutorial

Let's suppose we want to run an experiment on a couple different
experimental conditions---the regular "Control" chemistry and a
special new "Sparkly" chemistry from the lab.  We need to shepherd
these runs through some secondary analysis and then do some
comparative analysis on the conditions.

First, we build a CSV file that describes the input sources and
associated variables. Here's the table we will encode in our CSV file
`inputs.csv`:

| Condition   |      RunCode | ReportsFolder | Genome    | p_Chemistry |
|-------------|--------------|---------------|-----------|-------------|
| ControlChem | 3150120-0001 |               | lambdaNEB | Control     |
| ControlChem | 3150120-0002 |               | lambdaNEB | Control     |
| SparklyChem | 3150128-0001 |               | lambdaNEB | Sparkly     |
| SparklyChem | 3150128-0002 |               | lambdaNEB | Sparkly     |


A couple things to note.  First, note that we are "symbolically"
referring to locations of input data using the PacBio-internal
"RunCode" idiom.  There are other ways to refer to input data---see
the [condition table specification][condition-table-spec].

Also note that we've specified multiple runs for each condition; this
means that those input data will be combined within an analysis
condition.

Now, we want to see some basic comparisons based on
mapping---alignment lengths, accuracy, etc.  We are going to feed this
table to the `bauhaus` command and it will generate a runnable mapping
workflow for us.  But first let's explicitly ask `bauhaus` to
*validate* the table for us, to make sure we've encoded things
correctly and the data can actually be found.

  ```sh
  $ bauhaus -t inputs.csv -w BasicMapping validate
  Validation and input resolution succeeded.
  ```

So far so good.  Note that we had to specify the "workflow" we want to
run, using `-w BasicMapping`.  In order to validate the table, we need to
know what workflow it will be used with---for example, the "Genome"
column is necessary for mapping-based workflows whereas it may not be
necessary for other workflows.

Next, we want to generate a runnable workflow for performing this mapping:

  ```sh
  $ bauhaus -t inputs.csv -w BasicMapping -o experiment1 generate
  Validation and input resolution succeeded.
  Runnable Workflow written to directory "experiment1"
  ```

Here we specified an output directory "experiment1" and called the
`generate` subcommand; the output directory was created and populated
as follows:

  ```sh
  $ tree experiment1
  experiment1
  ├── build.ninja
  ├── condition-table.csv
  ├── log
  └── run.sh

  1 directory, 3 files
  ```

The idea is that we can now just run the command `experiment1/run.sh`
and our workflow will execute.  Let's see how that works:

  ```sh
  $ cat experiment1/run.sh
  #!/bin/bash
  source /mnt/software/Modules/current/init/bash
  MODULEPATH=/pbi/dept/primary/modulefiles:$MODULEPATH

  module purge
  module load smrtanalysis/mainline
  module load gfftools/dalexander
  module load R/3.2.3-bauhaus

  THISDIR=$(cd "$(dirname "$0")" && pwd)
  cd $THISDIR
  export PATH=$THISDIR/scripts:$PATH

  ~dalexander/bin/ninja
  ```

Basically the idea is that the `run.sh` just sets up some paths to
software needed (SMRTanalysis, R, etc.) and then invokes the
`build.ninja` script using the [Ninja build tool][ninja].  Let's look
at that ninja script:


  ```sh
  $ cat experiment1/build.ninja
  # Variables
  ncpus = 8
  grid = qsub -sync y -cwd -V -b y -e log -o log
  gridSMP = $grid -pe smp $ncpus

  # Rules
  rule map
    command = $gridSMP pbalign --nproc $ncpus $in $reference $out

  rule mergeAlignmentSetsForCondition
    command = $grid dataset merge $out $in


  # Build targets
  build ControlChem/mapping/m54006_160304_234053.alignmentset.xml: map $
      /pbi/collections/315/3150120/r54006_20160304_233237/1_A01/m54006_160304_234053.subreadset.xml
    ncpus = 8
    reference = $
        /mnt/secondary/iSmrtanalysis/current/common/references/lambdaNEB/sequence/lambdaNEB.fasta
   ... snip
   ```

The ninja script specifies quite explicitly the commands to run in
order to make the workflow "go".  Commands are run on the cluster
where appropriate, but "synchronously" so that the foreground terminal
never loses control---if you do `Ctrl-C` after starting `run.sh`, the
workflow will shut down immediately.  And, naturally, if we restart
the ninja run after this, it will pick up where it left
off---resumption "just works" in ninja.

Another benefit of using an external build script format is that we
can leverage some of the tooling that others have developed.  So for
example, we can see a graph of what this workflow will do, using the
`ninja -t graph` command and graphviz:


  ```sh
  $ ninja -C experiment1 -t graph | dot -Tpng > build.png
  ```

Here's the `build.png`:

![simple mapping workflow graph](./img/simple-mapping.png)

So you can see that each run is mapped and then within each condition
the mapped outputs are combined (in a lightweight fashion, using the
["dataset" concept][dataset]---this just involves small operations to
an XML file).  Two observations about this: first, the mapping jobs
will be quite slow since they are using a single node to map subreads
from up to a million ZMWs. Second, the job isn't very interesting
because it doesn't result in any comparative analyses at all.

So let's try running a more realistic workflow, which will run more
efficiently: consensus coverate titration.  In consensus coverage
titration, we see how many residual errors remain from performing
consensus at different coverage levels.  Note that scientifically
speaking, lambda phage would not be a good genome for this analysis
because it's too small a genome.  But we're just showcasing `bauhaus`,
so let's just do it:


  ```sh
  $ bauhaus -t inputs.csv -w CoverageTitrationReports -o experiment2 generate
  ```

Here's a graph of the worfklow script:

![coverage titration reports workflow graph](./img/coverage-titration-reports.png)

What's going on here is that we are 1) (conceptuallly, again using the
dataset concept) splitting each input run into parts; 2) mapping each
part; 3) physically combining the mapped parts of each input run (an
unnecessary step, arguably); 4) *conceptually* (dataset...) combine
the mapped inputs for each condition; 5) perform unfiltered variant
calling at different coverage levels; 6) filter out variants in
genomic regions where the reference genome is suspect (this is an
unfortunate reality since some of our references generated a few years
ago have misassemblies); 7) finally, summarize the filtered variant
counts and types across conditions using an R script.  Quite a bit
more complicated.  But at the end of the day, upon running, the user
will get a PDF summarizing the analysis.


[ninja]: http://ninja-build.org/
[condition-table-spec]: ./ConditionTableSpec.org
[dataset]: http://pacbiofileformats.readthedocs.org/en/3.0/DataSet.html
