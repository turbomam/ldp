library(reshape)

pdbXchain.names <- names(atom.dat)
prep4xyz <- lapply(pdbXchain.names, function(one.pdbXchain) {
  print(one.pdbXchain)
  temp <- atom.dat[[one.pdbXchain]]
  rownames(temp) <- temp$resno
  temp <- temp[, c('x', 'y', 'z')]
})
names(prep4xyz) <- names(atom.dat)


max.hood.size <- 9

pair.dists <- lapply(prep4xyz, function(current.pdb) {
  CA.dists <- dist.xyz(current.pdb)
  CAd.melt <- melt(CA.dists)
  CAd.melt$hops <- CAd.melt$X1 - CAd.melt$X2
  within.max <-
    CAd.melt[CAd.melt$hops > 0 &
               CAd.melt$hops <= max.hood.size , ]
  return(within.max)
})

actual.hood.size <- 4

actual.ldps <- lapply(pair.dists, function(current.pdb.strand) {
  # current.pdb.strand <- pair.dists[['1DF8-A']]
  
  current.pdb.strand <-
    current.pdb.strand[current.pdb.strand$hops <= actual.hood.size , ]
  
  if (nrow(current.pdb.strand) > 0) {
    hood.sum <- aggregate(
      current.pdb.strand$value,
      by = list(current.pdb.strand$X2),
      FUN = sum,
      na.rm = TRUE
    )
    names(hood.sum) <- c('resno', 'distsum')
    
    hood.contributors <- aggregate(
      current.pdb.strand$value,
      by = list(current.pdb.strand$X2),
      FUN = length
    )
    names(hood.contributors) <- c('resno', 'contributors')
    hood.res <- merge(x = hood.sum,
                      y = hood.contributors,
                      by = 'resno')
    
    return(hood.res)
  }
  
})

actual.min.max.res <-
  lapply(actual.ldps, function(current.pdb.strand) {
    return(current.pdb.strand$resno[current.pdb.strand$contributors == actual.hood.size])
  })

actual.min.max.res <- unlist(actual.min.max.res)
summary(actual.min.max.res)

# some pdb entries use "inflated" residue number for chains > A ?
hist(actual.min.max.res, breaks = 99)

which.max(actual.min.max.res)
actual.min.max.res[[which.max(actual.min.max.res)]]
