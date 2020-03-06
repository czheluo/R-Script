read.plink <- function(bed, bim, fam, na.strings=c("0", "-9"), sep=".",
                       select.subjects=NULL, select.snps=NULL) {
  lb <- nchar(bed);
  ext <- substr(bed, lb-3, lb);
  if (ext==".bed") {
    stub <- substr(bed, 1, lb-4)
  } else {
    stub <- bed
    bed <- paste(bed, ".bed", sep="")
  }
  if (missing(bim))
    bim <- paste(stub, ".bim", sep="")
  if (missing(fam))
    fam <- paste(stub, ".fam", sep="")
  df.bim <- read.table(bim, comment.char="", as.is=TRUE, na.strings=na.strings)
  snps <- as.character(df.bim[,2])
  if (!is.null(select.snps)) {
    if (is.character(select.snps)) {
      select.snps <- match(select.snps, snps)
      if (any(is.na(select.snps)))
        stop("unrecognised snp selected")
    }
    else if (is.numeric(select.snps)) {
      if (min(select.snps)<1 || max(select.snps)>nrow(df.bim))
        stop("selected snp subscript out of range")
      select.snps <- as.integer(select.snps)
    }
    else
      stop("unrecognised snp selection")
    if (any(duplicated(select.snps)))
      stop("snp selection contains duplicates")
    select.snps <- sort(select.snps)
    df.bim <- df.bim[select.snps,]
    snps <- snps[select.snps]
  }     
  if (ncol(df.bim)==6) {
    names(df.bim) <- c("chromosome", "snp.name", "cM", "position",
                       "allele.1", "allele.2")
  }
  else {
    warning("non-standard .bim file")
  }
  rownames(df.bim) <- snps
  if (any(duplicated(snps)))
    stop("duplicated snp name(s)")
  df.fam <- read.table(fam, comment.char="", as.is=TRUE, na.strings=na.strings)
  if (!is.null(select.subjects)) {
    if (is.numeric(select.subjects)) {
      if (min(select.snps)<1 || max(select.snps)>nrow(df.fam))
        stop("selected subject subscript out of range")
      select.subjects <- as.integer(select.subjects)
    }
    else
      stop("unrecognised subject selection")
    if (any(duplicated(select.subjects)))
      stop("subject selection contains duplicates")
    select.subjects <- sort(select.subjects)
    df.fam <- df.bim[select.subjects,]
  }     
  ped <- as.character(df.fam[,1])
  mem <- as.character(df.fam[,2])
  if (any(duplicated(ped))) {
    if (any(duplicated(mem))) {
      id <- paste(ped, mem, sep=sep)
      if (any(duplicated(id)))
        stop("couldn't create unique subject identifiers")
    }
    else
      id <- mem
  }
  else
    id <- ped
  names(df.fam) <-  c("pedigree", "member", "father", "mother", "sex",
                      "affected")
  rownames(df.fam) <- id
  gt <- .Call("readbed", bed, id, snps, select.subjects, select.snps,
              PACKAGE="snpStats")
  list(genotypes=gt, fam=df.fam, map=df.bim)
}

write.plink <- function(file.base, snp.major=TRUE, snps,
                        subject.data, pedigree, id, father, mother, sex, phenotype,
                        snp.data, chromosome, genetic.distance, position, allele.1, allele.2,
                        na.code=0) {
  
  mcall <- match.call()
  if (!is(snps, "SnpMatrix"))
    stop("snps argument must be of class SnpMatrix")
  X <- is(snps, "XSnpMatrix")
  nr.snps = nrow(snps)
  nc.snps = ncol(snps)
  if (missing(subject.data)) { ## subject data are in calling environment
    
    if (missing(pedigree))
      pedigree <- rownames(snps)
    else {
      pedigree <- as.character(pedigree)
      len <- length(pedigree)
      if (len != nr.snps)
        stop("length of `pedigree' argument incompatible with `snps'")
      if (any(is.na(pedigree)))
        stop("NA illegal for pedigree identifier")
    }
    
    if (missing(id))
      id <- rep("1", nr.snps)
    else {
      id <- as.character(id)
      if (length(id)!=nr.snps)
        stop("length of `id' argument  incompatible with `snps'")
      if (any(is.na(id)))
        stop("NA illegal for subject identifier within pedigree")
    }
    
    if (missing(father))
      father <- rep(na.code, nr.snps)
    else {
      father <- as.character(father)
      if (length(father)!=nr.snps)
        stop("length of `father' argument incompatible with `snps'")
      father[is.na(father)] = as.character(na.code)
    }
    
    if (missing(mother))
      mother <- rep(na.code, nr.snps)
    else {
      mother <- as.character(mother)
      if (length(mother)!=nr.snps)
        stop("length of `mother' argument incompatible with `snps'")
      mother[is.na(mother)] = as.character(na.code)
    }
    
    if (missing(sex)) {
      if (X)
        sex <- ifelse(snps@diploid, 2, 1)
      else
        sex <- rep(na.code, nr.snps)
    }
    else {
      sex <- as.numeric(sex)
      sex[is.na(sex)] = as.numeric(na.code)
    }
    
    if (missing(phenotype))
      phenotype <- rep(na.code, nr.snps)
    else {
      phenotype <- as.numeric(phenotype)
      phenotype[is.na(phenotype)] <- as.numeric(na.code)
    }
    
  } else { ## ped data are in data dataframe
    
    subject.data <- as.data.frame(subject.data)
    sord <- match(rownames(snps), rownames(subject.data))
    if (any(is.na(sord)))
      stop("missing entries in `subject.data' frame")
    
    if (missing(pedigree))
      pedigree <- rownames(snps)
    else {
      pedigree <- as.character(eval(mcall$pedigree, envir=subject.data))[sord]
      if (any(is.na(pedigree)))
        stop("NA illegal for pedigree identifier")
    }
    
    if (missing(id))
      id <-  rownames(snps)
    else {
      id <- as.character(eval(mcall$id, envir=subject.data))[sord]
      if (any(is.na(id)))
        stop("NA illegal for subject identifier within pedigree")
    }
    
    if (missing(father))
      father <- rep(na.code, nr.snps)
    else {
      father <- as.character(eval(mcall$father, envir=subject.data))[sord]
      father[is.na(father)] = as.character(na.code)
    }
    if (missing(mother))
      mother <- rep(na.code, nr.snps)
    else {
      mother <- as.character(eval(mcall$mother, envir=subject.data))[sord]
      mother[is.na(mother)] = as.character(na.code)
    }
    if (missing(sex))
      sex <- rep(na.code, nr.snps)
    else {
      sex <-  as.numeric(eval(mcall$sex, envir=subject.data))[sord]
      sex[is.na(sex)] = as.numeric(na.code)
    }
    
    if (missing(phenotype))
      phenotype <- rep(na.code, nr.snps)
    else {
      phenotype <- as.numeric(eval(mcall$phenotype, envir=subject.data))[sord]
      phenotype[is.na(phenotype)] = as.character(na.code)
    }
  }
  
  ## Map file
  
  if (missing(snp.data)) { ## snp data are in calling environment
    if (missing(chromosome)) {
      if (X)
        chromosome <- rep(23, nc.snps)
      else
        chromosome <- rep(na.code, nc.snps)
    }
    else {
      len <- length(chromosome)
      if (is.character(chromosome)) {
        chromosome = match(chromosome,
                           c(as.character(1:22), "X", "Y", "XY", "MT"))
        if (any(is.na(chromosome)))
          stop("unrecognized chromosome name")
      }
      else
        chromosome <- as.numeric(chromosome)
      if (X && any(chromosome!=23))
        stop("chromosome argument conflicts with snp data type")
      if (len==1) 
        chromosome <- rep(as.numeric(chromosome), nc.snps)
      else if (len!=nc.snps)
        stop("length of `chromosome' argument incompatible with `snps'")
      chromosome[is.na(chromosome)] <- as.integer(na.code)
    }
    
    if (missing(genetic.distance))
      genetic.distance <- rep(na.code, nc.snps)
    else {
      len <- length(genetic.distance)
      if (len==nc.snps)
        genetic.distance <- as.numeric(genetic.distance)
      else
        stop("length of `genetic.distance' argument incompatible and `snps'")
      genetic.distance[is.na(genetic.distance)] <- as.numeric(na.code)
    }
    
    if (missing(position))
      position <- rep(na.code, nc.snps)
    else {
      len <- length(position)
      if (len==nc.snps)
        position <- as.numeric(position)
      else
        stop("length of `position' argument incompatibe and `snps'")
      position[is.na(position)] <- as.numeric(na.code)
    }
    
    if (missing(allele.1))
      allele.1 <- rep("A", nc.snps)
    else {
      len <- length(allele.1)
      if (len==nc.snps)
        allele.1 <- as.character(allele.1)
      else
        stop("length of `allele.1' argument incompatibe and `snps'")
      allele.1[is.na(allele.1)] <- "A"
    }
    
    if (missing(allele.2))
      allele.2 <- rep("B", nc.snps)
    else {
      len <- length(allele.2)
      if (len==nc.snps)
        allele.2 <- as.character(allele.2)
      else
        stop("length of `allele.2' argument incompatibe and `snps'")
      allele.2[is.na(allele.2)] <- "A"      
    }
  } else { ## snp data are in data dataframe
    
    snp.data <- as.data.frame(snp.data)
    sord <- match(colnames(snps), rownames(snp.data))
    if (any(is.na(sord)))
      stop("missing entries in `snp.data' frame")
    
    if (missing(chromosome))
      chromosome <- rep(na.code, nc.snps)
    else {
      chromosome <- as.character(eval(mcall$chromosome, envir=snp.data))[sord]
      chromosome[is.na(chromosome)] <- as.character(na.code)
    }
    
    if (missing(genetic.distance))
      genetic.distance <- rep(na.code, nc.snps)
    else {
      genetic.distance <-
        as.numeric(eval(mcall$genetic.distance, envir=snp.data))[sord]
      genetic.distance[is.na(genetic.distance)] <- as.numeric(na.code)
    }
    
    if (missing(position))
      position <- rep(na.code, nc.snps)
    else {
      position <- as.numeric(eval(mcall$position, envir=snp.data))[sord]
      position[is.na(position)] <- as.numeric(na.code)
    }
    
    if (missing(allele.1))
      allele.1 <- rep("A", nc.snps)
    else {
      allele.1 <- as.character(eval(mcall$allele.1, envir=snp.data))[sord]
      allele.1[is.na(allele.1)] <- "A"
    }
    
    if (missing(allele.2))
      allele.2 <- rep("B", nc.snps)
    else {
      allele.2 <- as.character(eval(mcall$allele.2, envir=snp.data))[sord]
      allele.2[is.na(allele.2)] <- "B"
    }
  }
  
  ## Write files
  
  famdf <- data.frame(pedigree, id, father, mother, sex, phenotype)
  famfn <- paste(file.base, "fam", sep=".")
  cat("Writing FAM file to", famfn, "\n")  
  write.table(famdf, file=famfn, row.names=FALSE, col.names=FALSE,
              quote=FALSE, sep="\t")
  
  
  mapdf <- data.frame(chromosome, colnames(snps), genetic.distance,
                      position, allele.1, allele.2)
  mapfn <- paste(file.base, "bim", sep=".")
  cat("Writing extended MAP file to", mapfn, "\n")
  write.table(mapdf, file=mapfn, row.names=FALSE, col.names=FALSE,
              quote=FALSE, sep="\t")
  
  bedfn <- paste(file.base, "bed", sep=".")
  cat("Writing BED file to ", bedfn, " (", sep="")
  if (snp.major)
    cat("SNP-major mode)\n")
  else
    cat("Subject-major mode)\n")
  .Call("writebed", snps, bedfn, snp.major, PACKAGE="snpStats")
}
