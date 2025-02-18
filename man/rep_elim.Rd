% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rep_elim.R
\name{rep_elim}
\alias{rep_elim}
\title{Elimination of replications within an experiment}
\usage{
rep_elim(processed_data, infoTables)
}
\arguments{
\item{processed_data}{An expression matrix with normalized expression values or a list of them.}

\item{infoTables}{A character matrix with information about experiment (or a list of them). The row names of it must be the same as column
names in expression matrix (first argument) and the column names must be 'Experiment', 'Cell line', 'Treatment', 'Time', 'Dose'
exactly in this order.}
}
\value{
Function returns list with two elements. The first one is an expression matrix without replications with changed column names.
The second one is a character matrix with information about UNIQUE samples and is needed to create stability ranking.
}
\description{
\code{rep_elim} function averages expression values for replications within an experiment and prepare expression matrices for being used
in \code{create_ranking} function. It also returns information about unique samples in each experiment what is also needed for
\code{create_ranking} function to run.

If there are no replications in the data \code{rep_elim} function will just prepared objects to use with \code{create_ranking} function.
}
\details{
Example table structure (second argument) should look like this:

\tabular{cccccc}{
                   \tab Experiment  \tab  Cell line  \tab  Treatment   \tab   Time   \tab   Dose   \cr
 array1.txt        \tab      Exp1    \tab   HCT116   \tab      C      \tab    0    \tab    0   \cr
 array2.txt        \tab      Exp1    \tab   HCT116   \tab      C      \tab     0    \tab     0    \cr
 array3.txt        \tab      Exp1    \tab   HCT116   \tab     IR      \tab   1    \tab    4    \cr
 array4.txt        \tab      Exp1    \tab    HCT116   \tab     IR      \tab    2    \tab     4    \cr
 array5.txt        \tab      Exp1    \tab    HCT116   \tab      IR     \tab     4    \tab    2
 }

 Note that in treatment column it shoul be C for control arrays and IR for treated ones (even if the samples were treated with
 for example chemicals and not with ionizing radiation).

 If you have downloaded data from ArrayExpress database you can find information needed for creating above in .sdrf files.
}
\examples{
\dontrun{
##### example with some experiments downloaded from ArrayExpress database

# download data from ArrayExpress database
to_download = c("E-GEOD-67309", "E-MTAB-966")
my_data = downloadAE(to_download, getwd())

# load data
platforms = c("Affymetrix", "Agilent")
loaded_data = load_multi_data(my_data, platforms)

# normalize and annotate
norm_data = multi_norm_and_annot(loaded_data$raw_expression_data, platforms)

# prepare tables for rep_elim function as shown in details
path_to_tables = system.file("inst/extdata", "tables_ex3.rds", package = "FindReference")
my_tables = readRDS(path_to_tables)

# eliminate replications and prepare object for create_ranking function
no_rep_data = rep_elim(norm_data, my_tables)


##### example with miRNA microarray data

# download data from ArrayExpress database
datamiRNA = downloadAE("E-MTAB-5197", "/home/emarek/")

# prepare table as shown in details load_miRNA help page
path_to_table = system.file("inst/extdata", "miRNA_ex1.rds", package = "FindReference")
my_table = readRDS(path_to_table)

# load data
loaded_data = load_miRNA(my_table, datamiRNA[[1]]$path)

# normalize and annotate data
norm_data = norm_and_annot_miRNA(loaded_data)

# eliminate replications and prepare object for create_ranking function
no_rep_data = rep_elim(norm_data, my_table)
}

}
\seealso{
\code{\link{create_ranking}}
}
