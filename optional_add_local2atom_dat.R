library(bio3d)

load("C:/Users/mark/Desktop/ldp/atom_dat.Rdata")

# setwd("c:\\Users\\mark\\Documents")
# list.files()
# setwd("wang4alks")

local.pdbs <-
  list.files(path = "c:\\Users\\mark\\Documents\\wang4alks",
             pattern = "*.pdb",
             full.names = TRUE)

# temp <-
#   read.pdb("c:\\Users\\mark\\Documents\\wang4alks\\cluster.out.rep.c0.pdb")
#
# paste0(aa321(temp$atom$resid[temp$atom$elety == "CA"]), collapse = "")

name.bits <-
  strsplit(x = local.pdbs, split = ".", fixed = TRUE)

last.bit <- lapply(name.bits, function(one.name) {
  temp <- length(one.name)
  return(one.name[temp - 1])
})


temp <- lapply(local.pdbs, function(current.local) {
  # current.local <- local.pdbs[[1]]
  print(current.local)
  name.bits <-
    strsplit(x = current.local,
             split = ".",
             fixed = TRUE)
  name.bits <- name.bits[[1]]
  # print(name.bits)
  temp.len <- length(name.bits)
  # print(temp.len)
  basename.last <- name.bits[temp.len - 1]
  as.was <- read.pdb(current.local)
  new.fn <- paste0(basename.last, ".pdb")
  print(new.fn)
  
  as.was$atom$resno <- as.was$atom$resno + 1092
  
  # write.pdb(pdb = as.was, file = new.fn)
  # writes same thing... not necessary
  # just amke nice names
  #
  
  # paste0(aa321(as.was$atom$resid[as.was$atom$elety == "CA"]), collapse = "")
  
  awat <- as.was$atom
  awat <- awat[awat$elety == "CA" & awat$type == "ATOM" ,]
  
})

names(temp) <- last.bit

# local.pdbs[[1]]

atom.dat <- append(temp, atom.dat)

save(atom.dat, file = "atom_dat.Rdata")
