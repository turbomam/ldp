library(reshape)
library(bio3d)

Sys.setenv(JAVA_HOME = 'C:\\Program Files\\Java\\jre1.8.0_111')

# for daisy/gower clustering of LDP sums, even with NAs
# distances didn't seem that sueful for finding interesing paris
# coudl probably take out of production
library(cluster)

# for fast df -> list conversion of ldp diff results
library(purrr)

# only precompute linear distances wiht a span of max.hood.size or less
max.hood.size <- 9

distfunc <- function(x)
  daisy(x, metric = "gower")


load("C:/Users/mark/Desktop/ldp/atom_dat.Rdata")
load("C:/Users/mark/Desktop/ldp/")

pdbXchain.names <- names(atom.dat)
prep4xyz <- lapply(pdbXchain.names, function(one.pdbXchain) {
  print(one.pdbXchain)
  temp <- atom.dat[[one.pdbXchain]]
  rownames(temp) <- temp$resno
  temp <- temp[, c('x', 'y', 'z')]
})
names(prep4xyz) <- names(atom.dat)


# max.hood.size <- 9

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

if (interactive()) {
  # show which residue numebers are commonly used in the loaded entries
  print(summary(actual.min.max.res))
  # some pdb entries use "inflated" residue number for chains > A ?
  hist(actual.min.max.res, breaks = 99)
  
  print(which.max(actual.min.max.res))
  print(actual.min.max.res[[which.max(actual.min.max.res)]])
}

### omitting some residue-numer accounting that was in asym_cleaner.R

range2check <- min(actual.min.max.res):max(actual.min.max.res)


all.sorted <- names(actual.ldps)
locateds.overlap <- lapply(all.sorted, function(current.first) {
  # current.first <- '1DF8-A'
  # print(current.first)
  temp <- actual.ldps[[current.first]]
  temp <- temp[temp$contributors == actual.hood.size ,]
  first.adequate.sums <- temp$resno
  # print(first.located)
  eligible <-
    all.sorted[all.sorted > current.first]
  # print(sort(eligible))
  inner.temp <- lapply(eligible, function(current.last) {
    # current.last <- '1DF8-B'
    
    temp <- actual.ldps[[current.last]]
    temp <- temp[temp$contributors == actual.hood.size ,]
    last.adequate.sums <- temp$resno
    
    inter.size <-
      length(intersect(first.adequate.sums, last.adequate.sums))
    union.size <-
      length(union(first.adequate.sums, last.adequate.sums))
    inter.rat.union <- inter.size / union.size
    return(list(current.first, current.last, inter.rat.union))
  })
  if (length(inner.temp) > 0) {
    inner.temp <- do.call(rbind.data.frame, inner.temp)
    names(inner.temp) <- c('first', 'last', 'ratio')
    return(inner.temp)
  }
})
locateds.overlap <- do.call(rbind.data.frame, locateds.overlap)

if (interactive()) {
  # show the degree to which different pairs of entries
  # share the same residue numbers
  hist(locateds.overlap$ratio, breaks = 99)
}


recon.frame <- cbind.data.frame(range2check)
names(recon.frame) <- "resno"
recon.frame$distsum <- NA
recon.frame$contributors <- 0
recon.frame$priority <- 2

### is this really useful... aren't the distances recalculated in real time
### in other places?

recon.ldps <- lapply(actual.ldps, function(one.ldp) {
  # one.ldp <- actual.ldps[[1]]
  one.ldp <- one.ldp[one.ldp$contributors == actual.hood.size ,]
  one.ldp$priority <- 1
  temp <- rbind.data.frame(one.ldp, recon.frame)
  
  aggdata <-
    aggregate(
      temp$priority,
      by = list(temp$resno),
      FUN = min,
      na.rm = TRUE
    )
  names(aggdata) <- c('resno', 'priority')
  
  priority.frame <- merge(x = aggdata,
                          y = temp,
                          by = c('resno', 'priority'))
  priority.frame <- priority.frame[order(priority.frame$resno) , ]
  # return(priority.frame)
  return(priority.frame$distsum)
})

recon.ldps <- do.call(rbind.data.frame, recon.ldps)
# names(recon.ldps) <- range2check
names(recon.ldps) <- range2check

ldp.dist <- distfunc(recon.ldps)

ldp.d.mat <- as.matrix(ldp.dist)

ldp.d.table <- melt(data = ldp.d.mat)

# plot(hclust(ldp.dist))
# may need some cleanup
# Error in hclust(ldp.dist) : NA/NaN/Inf in foreign function call (arg 11)

ldp.d.table$ordered <-
  as.character(ldp.d.table$X2) > as.character(ldp.d.table$X1)

ldp.d.table <- ldp.d.table[ldp.d.table$ordered, ]
ldp.d.table <- ldp.d.table[, 1:3]

# names(locateds.overlap)
# names(ldp.d.table)

names(ldp.d.table) <- c("first" , "last" , "distance")

ldp.d.table <- melt(data = ldp.d.mat)

# plot(hclust(ldp.dist))
# may need some cleanup
# Error in hclust(ldp.dist) : NA/NaN/Inf in foreign function call (arg 11)

ldp.d.table$ordered <-
  as.character(ldp.d.table$X2) > as.character(ldp.d.table$X1)

ldp.d.table <- ldp.d.table[ldp.d.table$ordered,]
ldp.d.table <- ldp.d.table[, 1:3]

# names(locateds.overlap)
# names(ldp.d.table)

names(ldp.d.table) <- c("first" , "last" , "distance")

reasonalbe.interesting <- merge(x = locateds.overlap,
                                y = ldp.d.table,
                                by = c("first" , "last"))
reasonalbe.interesting$combo.name <-
  paste0(reasonalbe.interesting$first,
         ' - ',
         reasonalbe.interesting$last)

if (interactive()) {
  plot(reasonalbe.interesting$ratio,
       reasonalbe.interesting$distance)
}

# it doesn't look like the distance metric in reasonable.interesting
# is a good indicator of two pdbs with different folds
# maybe it ahs soemthign to do with lots of NAs
# because a very small proportion of the entries use residue numbers
# greater than the protein length,
# so there are lots of mostly-NA columns?

# NEEDED performance improvement
# had been using two nested loops
# use multiple cores?  just think of somehting more effienct?

all.sorted <- names(actual.ldps)
system.time(all.diffs <-
              lapply(all.sorted, function(current.first) {
                # current.first <- '1DF8-A'
                print(current.first)
                first.ldp <-
                  as.numeric(recon.ldps[current.first, ])
                eligible <-
                  all.sorted[all.sorted > current.first]
                if (length(eligible) > 0) {
                  eligible <- recon.ldps[eligible,]
                  
                  diff.frame <-
                    as.data.frame(t(apply(eligible, 1, '-', first.ldp)))
                  
                  diff.list <-
                    flatten(by_row(
                      diff.frame,
                      ..f = function(x)
                        flatten_dbl(x),
                      .labels = FALSE
                    ))
                  names(diff.list) <- rownames(diff.frame)
                  
                  return(diff.list)
                  
                }
                
              }))
names(all.diffs) <- all.sorted

# calculate # NAs, mean & sd
diff.stats <- lapply(names(all.diffs), function(current.diff) {
  # current.diff <- '4GDA-B'
  print(current.diff)
  temp <- all.diffs[[current.diff]]
  temp.names <- names(temp)
  # print(temp.names)
  inner.temp <- lapply(temp.names, function(current.inner) {
    print(current.inner)
    diff.values <- all.diffs[[current.diff]][[current.inner]]
    # print(diff.values)
    diff.mean <- mean(diff.values, na.rm = TRUE)
    diff.sd <-   sd(diff.values, na.rm = TRUE)
    diff.nas <- sum(is.na(diff.values))
    diff.max <- max(diff.values, na.rm = TRUE)
    return(list(
      current.diff,
      current.inner,
      diff.mean,
      diff.sd,
      diff.nas,
      diff.max
    ))
  })
  if (length(inner.temp) > 0) {
    inner.temp <- do.call(rbind.data.frame, inner.temp)
    names(inner.temp) <-
      c('first', 'last', 'mean', 'sd', 'nas', 'max')
    return(inner.temp)
  }
  
})
diff.stats <- do.call(rbind.data.frame, diff.stats)

diff.stats$pct.res.located <-
  1 - (diff.stats$nas / max(diff.stats$nas))


### select a pair of ldps to plot as difference
### here, based on number of nas resulting from subtraction
### (including those outside of the uniprot expected length)
### and the maximum difference signal (more sensitive to single, extreme differences?)
### could also use sd (more sensitive to man smaller sfferences?)
### or use soem signal processing

### should also include amino acid percent id, caluclated below (pid.res)
###  pid.res uses unique sequence serial numbers

if (interactive()) {
  plot(diff.stats$nas, diff.stats$max)
}

# code follows in asym_cleaner.R
# for identifying interesting (high max LDP diff)
# but reaasonalbe (low # or NAs in terms of incompatible residue numbers)

combo.brief <-
  unique(combo.res[, c("pdbid", "cgph", "nonpolents", "polents")])

chains4pdbid <- unique(combo.res[, c("pdbid", "strand")])

# chains4pdbid <- stack(unstack(temp))
names(chains4pdbid) <- c("pdbid", "chains")

combo.brief <- merge(x = chains4pdbid,
                     y = combo.brief,
                     by = 'pdbid')

small.piv <- strsplit(x = combo.res$nonpolents, split = ', ')
names(small.piv) <- combo.res$pdbid
small.piv <- stack(small.piv)
names(small.piv) <- c("nonpolent", "pdbid")
small.piv$pdbid <- as.character(small.piv$pdbid)

# i feel like i used to know how to do this with cast
small.piv <-
  as.data.frame.matrix(table(small.piv$pdbid, small.piv$nonpolent))
# some crystals claim no water present!!!
small.piv <- small.piv[, !names(small.piv) %in% c("HOH")]

# smalls.per.entry  <- rowSums(small.piv)
# entries.per.small <- colSums(small.piv)

nonpol.dist <- as.matrix(dist(small.piv, method = "binary"))

save(combo.brief, nonpol.dist, small.piv, file = "nonpol_dist.Rdata")
