library(ggvis)
library(dplyr)

load("atom_dat.Rdata")
load("nonpol_dist.Rdata")


get.fav.small <- function(ref.entry, fav.small) {
  parta <- combo.res[, c('pdbid', 'cgph')]
  partb <-
    cbind.data.frame(rownames(nonpol.dist), nonpol.dist[, ref.entry])
  names(partb) <- c('pdbid', 'dist')
  mergedparts <- merge(
    x = parta,
    y = partb,
    by.x = 'pdbid',
    by.y = 'pdbid'
  )
  mergedparts$self <- FALSE
  mergedparts$self <- mergedparts$pdbid == ref.entry
  
  partc <-
    cbind.data.frame(rownames(small.piv), as.logical(small.piv[, fav.small]))
  names(partc) <- c('pdbid', fav.small)
  
  mergedparts <- merge(
    x = mergedparts,
    y = partc,
    by.x = 'pdbid',
    by.y = 'pdbid'
  )
  
  return(mergedparts)
}

temp <- get.fav.small(ref.entry = '3L9P', fav.small = '4LO')
names.temp <- names(temp)
names.temp[[5]] <- "fav.pres"
names(temp) <- names.temp

temp$id <- 1:nrow(temp)  # Add an id column to use ask the key

temp <- merge(x = temp,
              y = combo.brief,
              by.x = "pdbid",
              by.y = "pdbid")

# dput(names(temp))

all_values <- function(x) {
  if (is.null(x))
    return(NULL)
  # x$pdbid
  row <- temp[temp$id == x$id, ]
  # c("pdbid", "cgph.x", "dist", "self", "fav.pres", "id", "chains",
  #   "cgph.y", "nonpolents", "polents")
  row <- row[c("pdbid", "chains", "nonpolents", "polents")]
  
  paste0(names(row), ": ", format(row), collapse = "<br />")
  # return(x)
}

temp %>%
  ggvis( ~ cgph.x, ~ dist, key := ~ id) %>%
  layer_points(fill = ~ factor(fav.pres),
               shape = ~ factor(self)) %>%
  add_axis("x", title = "Crystal Growth pH") %>%
  add_axis("y", title = "Distance from 'self' by Ligand Usage")  %>%
  add_tooltip(all_values, on = "hover")

# add titles and axis labels
# legend looks funny/overlapping?
# hover tips