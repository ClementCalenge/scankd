#' scankd: A package for implementing the approach of Kulldorff (2001)
#'
#' The package scankd provides all the functions required to implement
#' the scan statistics developed in Kulldorff (2001). The main help
#' page of interest is \code{ngkd}, which describes how to use this
#' algorithm.
#'
#'
#' @docType package
#' @name scankd
NULL

#' @title Raster map of France.
#'
#' @details
#' Raster map of France.
#'
#'
#' @format An object of class \code{SpatialPixelsDataFrame} (package
#' sp).
#'
"francekd"



#' @title Implements The Approach of Kulldorff (2001)
#' @export
#'
#' @details \code{ngkd} prepares the data for the scan (identification
#' of neighbourhood for each pixel. \code{Kulldorff2001} implements the
#' scan statistics, and calculates the most likely cluster (size,
#' duration and likelihood ratio) associated to each pixel of the map.
#' \code{identifyClusters} identifies and plots the most likely
#' non-overlapping clusters on the map. \code{circleskd} can be used
#' to plot the clusters to an already existing image of the area.
#'
#' @param pts an object inheriting the class SpatialPixels of the
#' package sp, containing the grid of pixels used in the scan (each
#' pixel being the possible center of a circular cluster).
#' @param radius a numeric vector of positive and increasing values
#' corresponding to the possible cluster radiuses (ordered
#' increasingly) tested in the scan.  The units should be the same as
#' those used in \code{pts}.
#' @param time a numeric vector of positive and increading values
#' corresponding to the possible cluster durations tested in the scan
#' (ordered increasingly). The time units are not important, as these
#' values are not used in the calculations (only for display
#' purposes).
#' @param x for \code{prind.ngkd}, an object of class "ngkd" returned
#' by the function \code{ngkd}.  For \code{circleskd}, an object of
#' class \code{SpatialPointsDataFrame} indicating the location of the
#' cluster centers, returned by \code{identifyClusters}.
#' @param observed a data.frame with \code{Nt} columns and \code{Np}
#' rows (where \code{Nt} is the number of time units considered in the
#' scan, and \code{Np} the number of pixels in \code{pts}), containing
#' the number of cases observed in each pixel during each time unit.
#' @param theoretical a data.frame with \code{Nt} columns and
#' \code{Np} rows (where \code{Nt} is the number of time units
#' considered in the scan, and \code{Np} the number of pixels in
#' \code{pts}), containing the number of cases expected in each pixel
#' during each time unit under "normal" conditions.
#' @param ngkd an object of class "ngkd" returned
#' by the function \code{ngkd}.
#' @param loglambda logical. Should the actual value of the likelihood
#' ratio be returned, or its logarithm (the latter often results in
#' clearer maps).
#' @param spixdf an object of class \code{SpatialPixelsDataFrame}
#' returned by the function \code{Kulldorff2001}.
#' @param nclustMax the maximum number of clusters that should be
#' returned by the function.
#' @param plot logical. If \code{TRUE}, the function plots the
#' location of the most likely clusters on a map of the
#' (log-)likelihood ratios for each pixel.
#' @param text Either a missing value, in which case the function
#' indicates the order of the clusters (first, second, etc.) with
#' their duration at the center of the cluster; or a character string
#' indicating the variable in \code{x} that should be used to plot the
#' text at the center of the cluster (e.g. "radius", "time", "index",
#' "loglambda"); or a character vector with the same number of
#' elements as \code{nrow(x)} indicating for each cluster what text to
#' write at the center of the cluster.
#' @param \dots additional arguments to be passed to or from other
#' functions.
#'
#' @return \code{ngkd} returns a list of class \code{"ngkd"} required
#' by all the functions of the packages (stores all the spatial and
#' temporal information about the scan). \code{Kulldorff2001}
#' implements the scan statistics of Kulldorff (2001) and return, for
#' each pixel of the map, the (log-)likelihood ratio associated to the
#' best cluster centered on this pixel, as well as the diameter and
#' duration of the pixel (these variables are returned as a
#' \code{SpatialPixelsDataFrame}). \code{identifyClusters} returns a
#' \code{SpatialPointsDataFrame} containing the most likely
#' nonoverlapping clusters. Each row corresponds to the center of a
#' cluster, and the data slot contains the pixel index in the map used
#' to create the \code{ngkd} object corresponding to this point, the
#' size and duration of the cluster, as well as the (log-) likelihood
#' ratio corresponding to this pixel.
#'
#' @author Clement Calenge, \email{clement.calenge@@oncfs.gouv.fr}
#' @references Kulldorff, M. (2001) Prospective time periodic
#' geographical disease surveillance using a scan
#' statistic. \emph{Journal of the Royal Statistical Society: Series
#' A}, \bold{164}, 61-72.
#' @importFrom plotrix draw.circle
#' @import graphics
#' @import stats
#' @import methods
#' @import sp
#' @useDynLib scankd
#'
#' @examples
#' ## We have a map of France
#' image(francekd)
#'
#'
#' ##############################################
#' ##
#' ## Simulation of a dataset. Beginning of the simulations
#' ##
#' ## We simulate the following situation:
#' ## (i) a "normal noise" of cases during 52 weeks (in average, 0.001
#' ##     cases expected per pixel):
#' set.seed(9809)
#' cases <- as.data.frame(matrix(as.double(rpois(nrow(francekd)*52,0.001)),
#'                               nrow=nrow(francekd)))
#'
#' ## (ii) we add a cluster of "unusual" cases in the North-east of
#' ## France for the last 8 weeks
#' df <- structure(list(gg = c(160L, 123L, 230L, 368L, 333L, 231L,
#' 271L, 307L, 161L, 159L, 126L, 157L, 270L, 235L, 336L, 333L, 193L,
#' 303L, 307L), gi = c(3L, 1L, 6L, 2L, 8L, 4L, 4L, 7L, 4L, 6L, 2L, 6L,
#' 2L, 2L, 2L, 1L, 6L, 8L, 8L)), .Names = c("gg", "gi"), row.names =
#' c(NA, -19L), class = "data.frame")
#' for (i in 1:nrow(df)) {
#'    cases[df[i,1], df[i,2]] <- cases[df[i,1], df[i,2]]+1
#' }
#'
#' ## End of the simulations
#' ##############################################
#'
#' ## Show these data (sum over the last year)
#' su <- apply(cases,1,sum)
#' francekd$totalCases <- su
#' image(francekd[,"totalCases"])
#' ## Note the cluster in North-eastern France
#'
#' ## The object cases contains the number of cases recorded in each pixel
#' ## (rows) for each week (column) with the first column corresponding to
#' ## the most recent week.
#'
#' ## Suppose, for example, that we expect in theory 0.001 cases per
#' ## pixel for each week (the model might be more complex, but is
#' ## good enough for this demonstration). Therefore, the expected
#' ## number of cases is:
#' expect <- as.data.frame(matrix(rep(0.001,nrow(francekd)*52),
#'                                nrow=nrow(francekd)))
#'
#' ## Now, we scan with the following radiuses (in meters)
#' (rad <- c(1:6)*50000)
#'
#' ## And every week, for one year, i.e.
#' (weekno <- 1:52)
#' ## each element of weekno corresponds to a column in cases or in expect
#' ## the most recent week is stored first in this data.frame
#'
#' ## Prepare the data:
#' ng <- ngkd(francekd, rad, weekno)
#'
#' ## Scan:
#' kd <- Kulldorff2001(cases, expect, ng)
#'
#' ## identify Clusters
#' clust <- identifyClusters(kd, ng)
#'
#' ## Note that the function identifies correctly both the size and
#' ## the duration on he cluster.
#'
#' ## or, to show these clusters on the map of all cases:
#' image(francekd[,"totalCases"])
#' circleskd(clust, c("Best", "Second", "Third"))
ngkd <- function(pts, radius, time)
{
    if (!inherits(pts, "SpatialPixels"))
        stop("pts should inherit the class \"SpatialPixels\"")

    if (any(diff(radius)<=0))
        stop("The radiuses should be ordered by increasing size")

    if (any(diff(time)<=0))
        stop("The time lags should be ordered by increasing size")

    ## Distance matrix
    coo <- coordinates(pts)
    di <- as.matrix(dist(coo))

    ## neigh
    ngkd <- lapply(1:nrow(pts), function(i) {
                       if (sum(di[i,]<radius[1])==0)
                           stop("The first radius is smaller than cellsize. Increase it.")
                       vois <- lapply(radius, function(j) {
                                          idv <- c(1:ncol(di))[di[i,]<j]
                                          idv <- idv[idv!=i]
                                      })

                       for (j in 2:length(vois)) {
                           deja <- unlist(vois[c(1:(j-1))])
                          if (is.null(vois[[j]][!(vois[[j]]%in%deja)]))
                               stop("two virtually identical radiuses (difference between two successive radiuses smaller than cellsize).")
                           vois[[j]] <- vois[[j]][!(vois[[j]]%in%deja)]
                       }
                       vec <- unlist(vois)
                       j1 <- 0
                       j2 <- 0
                       cur <- 1
                       for (i in 1:length(radius)) {
                           j1[i] <- (cur)-1 ## le -1 est pour l'indice en C
                           j2[i] <- (cur+length(vois[[i]])-1)-1 ## idem
                           vec[(j1[i]+1):(j2[i]+1)] <- vois[[i]]
                           cur <- cur+length(vois[[i]])
                       }
                       df <- data.frame(as.integer(j1), as.integer(j2))
                       return(list(as.integer(vec-1), df)) ## idem, le - 1
                   })
    class(ngkd) <- "ngkd"
    attr(ngkd, "map") <- pts
    attr(ngkd, "radius") <- radius
    attr(ngkd, "time") <- time
    return(ngkd)
}


#' @rdname ngkd
#' @export
print.ngkd <- function(x, ...)
{
    if (!inherits(x, "ngkd"))
        stop("x should inherit the class ngkd")
    cat("**************************\n",
        "* Object of class ngkd\n",
        "* Neighborhood of each pixel for each radius\n",
        "* required by the function Kulldorff2001\n\n",
        "Number of pixels in the map: ", length(x),"\n",
        "Number of radiuses considered in the analysis: ",
        nrow(x[[1]][[2]]),"\n",
        "Number of time units considered in the analysis: ",
        length(attr(x,"time")),"\n",
        sep="")
}




#' @rdname ngkd
#' @export
Kulldorff2001 <- function(observed, theoretical, ngkd, loglambda=TRUE)
{
    if (!inherits(observed, "data.frame"))
        stop("observed should be of class data.frame")
    if (!inherits(theoretical, "data.frame"))
        stop("theoretical should be of class data.frame")

    if (ncol(observed)!=length(attr(ngkd,"time")))
        stop("the number of columns in observed does not correspond to the number of time units")
    if (nrow(observed)!=length(ngkd))
        stop("the number of units in observed does not correspond to the number of pixels")
    if (nrow(observed)!=nrow(theoretical))
        stop("the numbers of units in observed and theoretical do not match")
    if (ncol(observed)!=ncol(theoretical))
        stop("the numbers of column in observed and theoretical do not match")

    NparT <- as.double(cumsum(unlist(apply(observed,2,sum))))
    if (NparT[1]<0.000001)
        stop("There are zero cases during the first period.\nPlease increase the duration of this period")
    resu <- .Call("trouveCluster", ngkd, observed, theoretical, NparT, PACKAGE="scankd")

    resu <- as.data.frame(resu)
    resu[,2] <- attr(ngkd, "radius")[resu[,2]]
    resu[,3] <- attr(ngkd, "time")[resu[,3]]
    if (loglambda) {
        resu[,1] <- log(resu[,1])
        names(resu) <- c("loglambda", "radius","time")
    } else {
        names(resu) <- c("lambda", "radius","time")
    }
    coordinates(resu) <- coordinates(attr(ngkd,"map"))
    gridded(resu) <- TRUE
    return(resu)
}


#' @rdname ngkd
#' @export
identifyClusters <- function(spixdf, ngkd, nclustMax=nrow(spixdf), plot=TRUE)
{
    if (!inherits(spixdf, "SpatialPixelsDataFrame"))
        stop("spixdf should have been generated by Kulldorff2001")
    if (ncol(spixdf)!=3)
        stop("spixdf should have been generated by Kulldorff2001")
    if (!names(spixdf)[1]%in%c("loglambda","lambda"))
        stop("spixdf should have been generated by Kulldorff2001")
    if (!all(names(spixdf)[2:3]==c("radius","time")))
        stop("spixdf should have been generated by Kulldorff2001")

    index <- 1:nrow(spixdf)
    ma <- which.max(spixdf[[1]])
    lam <- spixdf[[1]][ma]
    dia <- spixdf[[2]][ma]
    tim <- spixdf[[3]][ma]
    ind <- ma
    coo <- coordinates(spixdf)
    cons <- rep(0, length(index))

    if (nclustMax>1) {
        di <- as.matrix(dist(coo))
        for (i in 1:(nclustMax-1)) {
            ## which pixels are already in a cluster?
            cons <- cons + as.numeric(di[ind[i],]<=dia[i])
            cons <- as.numeric(cons>0)
            if (sum(cons) == length(cons)) {
                break
            }
            ## Identify all the clusters with none of these pixels:
            cons[cons==0] <- sapply(c(1:nrow(di))[cons==0], function(j) {
                                        ## for each focus pixel outside the present clusters,
                                        ## consider if the cluster associated with the focus
                                        ## includes pixels in already identified clusters
                                        bad <- any((di[j,]<spixdf[[2]][j])&(cons==1))
                                        ## return
                                        return(bad)
                                    })
            if (sum(cons) == length(cons)) {
                break
            }
            ##
            spixdf2 <- spixdf[cons==0,]
            index2 <- index[cons==0]
            ma <- which.max(spixdf2[[1]])
            lam[i+1] <- spixdf2[[1]][ma]
            dia[i+1] <- spixdf2[[2]][ma]
            tim[i+1] <- spixdf2[[3]][ma]
            ind[i+1] <- index2[ma]
        }
    }

    df <- data.frame(index=ind, radius=dia, time=tim, lam=lam)
    names(df)[4] <- names(spixdf)[1]
    coordinates(df) <- coo[ind,,drop=FALSE]

    if (plot) {
        co <- coordinates(df)
        image(spixdf)
        for (i in 1:nrow(df)) {
            text(co[i,1], co[i,2], paste0(i, " (",df[[3]][i],")"))
            draw.circle(co[i,1], co[i,2], df[[2]][i])
        }
    }
    attr(df, "type") <- "clusterskd"
    return(df)
}


#' @rdname ngkd
#' @export
circleskd <- function(x, text=NA)
{
    if (!inherits(x,"SpatialPointsDataFrame"))
        stop("x should inherits the class SpatialPixelsDataFrame")
    if (is.null(attr(x,"type")))
        stop("x should have been generated by identifyClusters")
    if (attr(x,"type")!="clusterskd")
        stop("x should have been generated by identifyClusters")
    if (length(text)>1) {
        if (length(text)!=nrow(x))
            stop("text should be either of length 1 or the same as nrow(x)")
    }
    co <- coordinates(x)

    for (i in 1:nrow(x)) {
        if (length(text)==1) {
            if (is.na(text)) {
                text2 <- paste0(i, " (",x[[3]][i],")")
            } else {
                if (any(names(x)==text)) {
                    text2 <- x[[text]][i]
                } else {
                    text2 <- text
                }
            }
        } else {
            text2 <- text[i]
        }
        text(co[i,1], co[i,2], text2)
        draw.circle(co[i,1], co[i,2], x[[2]][i])
    }
}
