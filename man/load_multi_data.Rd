% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/load_multi_data.R
\name{load_multi_data}
\alias{load_multi_data}
\title{Loading microarray data from a list of experiments downloaded from ArrayExpress database}
\usage{
load_multi_data(dane, ExpInfoTable, sdrfFiles)
}
\arguments{
\item{dane}{A list of lists from downloadAE function.}

\item{platforma}{A character vector indicating platform for each experiment. Supported platforms are 'Affymetrix' and 'Agilent'.}

\item{donotread}{A list of character vectors indicating which microarrays should not be loaded for each experiment. Default is NA,
which means that all microarrays listed in dane$rawFiles will be loaded.
Note that list's length must be the same as the length of list 'dane'.}
}
\value{
Function returns a list with two elements. The first one is a list of objects specific for array platform with raw expression data.
The second element is a list of data frames with loaded .sdrf files for each experiment.
}
\description{
\code{load_multi_data} function loads data for a list of microarray experiments downloaded from ArrayExpress database.
}
\details{
This function is designed to read in microarray data from many experiments. For more detail see \code{load_data} function help.
}
\examples{
\dontrun{
### load all microarrays for two experiments
to_download = c("E-GEOD-21066", "E-MTAB-966")
platforms = c("Affymetrix", "Agilent")

dane = downloadAE(to_download, getwd())
loaded_experiments = load_multi_data(dane, platforms)

### do not load some microarrays from the first experiment
to_download = c("E-GEOD-21066", "E-MTAB-966")
platforms = c("Affymetrix", "Agilent")

dane = downloadAE(to_download, getwd())

# for second experiment all microarrays should be loaded
# -> the second element of the list is NA
unwanted_arrays = list(c("GSM526680.CEL", "GSM526756.CEL"), NA)

# or (which could be useful when there are a lot of experiments)
unwanted_arrays = rep(list(list(NA)), length(dane))
unwanted_arrays[[1]] = c("GSM526680.CEL", "GSM526756.CEL")

loaded_experiment = load_multi_data(dane, platforms, donotread = unwanted_arrays)

### when you downloaded data with downloadAE
### but didn't assign the output into variable

to_download = c("E-GEOD-21066", "E-MTAB-966")
platforms = c("Affymetrix", "Agilent")
downloadAE(to_download, getwd())

# read in data saved by downloadAE function
dane = readRDS('dataAE.rds')
loaded_experiment = load__multi_data(dane, platforms)

}

}
\seealso{
\code{\link{load_data}}
}
