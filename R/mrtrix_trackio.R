readTckHeader <- function(filename)
  {
    fp <- file(filename, open='rb')
    l <- readLines(fp, n=1)
    l <- gsub("[[:blank:]]+$", "", l)
    if (l != "mrtrix tracks") {
      print(l)
      stop("Not a track file")
    }
    hd <- NULL
    while (l != "END") {
      l <- readLines(fp, n=1)
      l <- gsub("[[:blank:]]+$", "", l)
      l <- gsub("^[[:blank:]]+", "", l)
      hd <- c(hd, l)
    }
    # ditch "END"
    hd <- head(hd, -1)
    hl <- strsplit(hd, ":")
    tags <- sapply(hl, function(x)x[1])
    vals <- sapply(hl, function(x)x[2])
    tags <- gsub("[[:blank:]]+$", "", tags)
    tags <- gsub("^[[:blank:]]+", "", tags)
    vals <- gsub("[[:blank:]]+$", "", vals)
    vals <- gsub("^[[:blank:]]+", "", vals)
    
    numvals <- as.numeric(vals)
    nn <- !is.na(numvals)
    ntags <- tags[nn]
    numvals <- numvals[nn]
    names(numvals) <- ntags
    ovals <- vals[!nn]
    names(ovals) <- tags[!nn]
    offset <- as.numeric(gsub("\\.", "", ovals["file"]))
    #browser()
    if (ovals["datatype"] != "Float32LE")
      stop("Only deals with floating point track files")
    close(fp)
    return(list(numeric=numvals, text=ovals, offset=offset))
  }

processTrackDefault <- function(tr)
  {
    ## sample track tool
    return(ncol(tr))
  }

trackLength <- function(trk)
  {
    x <- trk[1,]
    y <- trk[2,]
    z <- trk[3,]

    dx2 <- diff(x)^2
    dy2 <- diff(y)^2
    dz2 <- diff(z)^2

    length <- sum(sqrt(dx2 + dy2+ dz2))
    
  }
readTck <- function(filename, processTrack=processTrackDefault)
  {
    hd <- readTckHeader(filename)
    fp <- file(filename, open='rb')
    seek(fp, where=hd$offset)
    pbuff <- 1024 * 3
    prevdat <- NULL
    EOF <- FALSE
    res <- NULL
    repeat
      {
        incomplete <- TRUE
        while (incomplete) {
          ## keep reading until we have a 
          dat <- readBin(fp, n=pbuff, what="numeric", size=4)
          dim(dat) <- c(3, length(dat)/3)
          if (ncol(dat) < pbuff/3) {
            EOF <- TRUE
          }
            
          breaks <- which.max(is.na(dat[1,]))
          if (length(breaks) > 0) {
            incomplete <- FALSE
            dat <- cbind(prevdat, dat)
          }
        }
        ## now have a matrix with at least one complete track
        repeat
          {
            tend <- which(is.na(dat[1,]))
            if (length(tend) > 0) {
              tend <- tend[1]
              cc <- 1:(tend-1)
              track <- dat[, cc,drop=FALSE]
              dat <- dat[,-(1:tend),drop=FALSE]
              res <- c(res, processTrack(track))
            } else {
              prevdat <- dat
              break
            }
          }
        if (EOF) break
                       
      }
    close(fp)
    return(res)
  }
