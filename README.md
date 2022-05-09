A project by Kelly Jones and Andrew Howe to explore the viability of
the APClustering algorithm in single-cell data.  

--------------------------------------------------------------------------------

All source files are in the directory structure in the following fashion:

local/. contains files too large to be pushed to git

Source/. contains compressed raw pancreas data, which was processed into
tissue_counts/pancreas_counts.rds using the script Rscripts/curatepancreas.R.
To repeat this processing, the file must be decompressed manually on a local
disk

tissue_analysis/pancreas/. contains all files that were output by the script
Rscripts/experiment.R.  All .rds files that have exceptionally long runtime are
included in this submission for the grader's sake of time

tissue_counts/. contains the processed pancreas data from Source/.

Rscripts/. contains all R scripts. Rscripts/experiment.R is the main working R
script, but we have also included test.R in the base directory for the grader's
convenience.

--------------------------------------------------------------------------------

To run sample code and examine outputs, open the R file test.R in the base
directory.  Make sure that the working directory is set to the base directory of
the submission on your local device. In order to run the code, the following
packages will have to be installed:

```
install.packages(c("cluster", "ggplot2", "devtools", "Seurat", 
                   "pheatmap", "BiocManager", "RColorBrewer", 
		   "stringr", "apcluster", "igraph", "Matrix",
                   "tidyverse", "mcclust", "philentropy"))
BiocManager::install("viper")
BiocManager::install("biomaRt")
devtools::install_github("JEFworks/MUDAN")
devtools::install_github(repo = "califano-lab/PISCES", force = TRUE, build_vignettes = TRUE)
devtools::install_github('https://github.com/arh2207/APTestPipeline')
devtools::install_github('https://github.com/jchiquet/aricode')
```

Disclosure: the APTestPipeline package was released by the author Andrew as part of 
a separate experiment on the APClustering algorithm for the class ML4FG.  From
the package, we call a few functions for conveniency in parts of this analysis.
No code from the ML4FG was duplicated in the submission for this class, and none
of the functions provided in the package are submitted as part of this
assignment.  Furthermore, the submitted papers do not reuse any material between
projects.

--------------------------------------------------------------------------------

Important notes:

Objects which take a long time to run have their outputs pre-computed and functions
commented out.  Functions may be un-commented to test.  

You must set the working directory at the top of R scripts.  Examples are shown from 
analyses done on local machines of the authors of this experiment

Package dependency installation and management tends to be difficult on various machines,
and troubleshooting is often peculiar to your operating system.  For example, users running
MacOS will typically have to install certain packages via homebrew and manage package
include paths.  Some of these modifications are performed automatically on Windows environments
and Linux environments.  We highly recommend using RStudio as opposed to running R scripts
using a barebones interpreter.

A grader who wants to see a simplified version of our analysis should look at test.R
in the base directory.  We have removed all long runtime functions, but there are still
some that may take a moderate time to run (~2 minutes).  It may help to know these functions
in advance:
SCTransform
RunUMAP
RunTSNE
