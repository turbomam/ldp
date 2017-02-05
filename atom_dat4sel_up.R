library(rrdf)
library(bio3d)

# ALK_HUMAN
selected.up <- '<http://purl.uniprot.org/uniprot/Q9UM73>'

# CBPA1_BOVIN
# selected.up <- '<http://purl.uniprot.org/uniprot/P00730>'

# selected.up <- readline(prompt = 'enter UniProt IRI: ')

###   ###   ###

###   ###   ###

# searching v40 and v42 namespaces seperately
# and then unioning
# i have a funny feeling about this

pdb2strand.template <- 'PREFIX PDBo: %s
# PDB.val goes above
PREFIX PDBr: <http://rdf.wwpdb.org/pdb/>
PREFIX rdf:<http://www.w3.org/1999/02/22-rdf-syntax-ns#>
PREFIX dcterms:<http://purl.org/dc/terms/>
PREFIX xsd:<http://www.w3.org/2001/XMLSchema#>
SELECT distinct ?pdbid ?strand ?seq ?epseq
FROM <http://rdf.integbio.jp/dataset/pdbj>
WHERE {
?structref a PDBo:struct_ref .
?structref PDBo:of_datablock ?db .
?db dcterms:identifier ?pdbid .
?structref PDBo:reference_to_entity ?ref2ent .
# ?ref2ent PDBo:entity.id ?entid .
# ?ref2ent PDBo:entity.pdbx_description ?desc .
# ?ref2ent PDBo:entity.pdbx_number_of_molecules ?molct .
# ?ref2ent PDBo:entity.src_method ?srcmeth .
# ?ref2ent PDBo:entity.type ?enttype .
?ref2ent PDBo:referenced_by_entity_poly ?rbyentpoly .
# ?rbyentpoly PDBo:entity_poly.entity_id ?enpid .
# nonstandard monomers and linkages
# one letter code
# ?rbyentpoly PDBo:entity_poly.pdbx_seq_one_letter_code_can ?canseq .
?rbyentpoly PDBo:entity_poly.pdbx_strand_id ?strand .
optional { ?rbyentpoly PDBo:entity_poly.pdbx_seq_one_letter_code_can ?epseq } .
# ?rbyentpoly PDBo:entity_poly.type ?polytype .
# ?ref2ent PDBo:referenced_by_entity_src_gen ?rbysource .
# ?rbysource PDBo:entity_src_gen.pdbx_gene_src_scientific_name ?sciname .
# entity id, genus, taxids
# ?ref2ent PDBo:referenced_by_struct_asym ?rbyasym .
# ?rbyasym PDBo:struct_asym.id ?asymid .
# asym ent id?, baknk chain id, modified
# ?ref2ent PDBo:referenced_by_struct_ref ?rbysr .
# ?structref PDBo:struct_ref.id ?srid .
# ?structref PDBo:struct_ref.entity_id ?entid .
optional { ?structref PDBo:struct_ref.pdbx_seq_one_letter_code ?seq } .
# selected.up goes here
?structref PDBo:link_to_uniprot %s
}'


# ACTUALLY SHOUD REPORT ALL ENTITIES, SO AS TO INCLUDE PEPTIDE LIGANDS

pdb.ph.het.template <- 'PREFIX PDBo: %s
# PDB.val goes here
PREFIX PDBr: <http://rdf.wwpdb.org/pdb/>
PREFIX rdf:<http://www.w3.org/1999/02/22-rdf-syntax-ns#>
PREFIX dcterms:<http://purl.org/dc/terms/>
PREFIX xsd:<http://www.w3.org/2001/XMLSchema#>
SELECT ?pdbid ?cgph
(group_concat(distinct ?enplab ; separator = ", ") AS ?nonpolents)
(group_concat(distinct ?eplab ; separator = ", ") AS ?polents)
FROM <http://rdf.integbio.jp/dataset/pdbj>
WHERE {
?structref a PDBo:struct_ref .
?structref PDBo:of_datablock ?db .
?db dcterms:identifier ?pdbid .
optional { ?cg PDBo:of_datablock ?db .
?cg rdf:type PDBo:exptl_crystal_grow .
?cg PDBo:exptl_crystal_grow.pH ?cgph .
# bind(xsd:double(?cgph) as ?phval) .
} .
optional { ?enp rdf:type PDBo:pdbx_entity_nonpoly  .
?enp PDBo:of_datablock ?db  .
?enp PDBo:pdbx_entity_nonpoly.comp_id ?enplab .
# names may be hard to parse because they can contaion just about any character
# ?enp PDBo:pdbx_entity_nonpoly.name ?enplab .
} .
optional { ?ep rdf:type PDBo:entity_poly .
?ep PDBo:of_datablock ?db  .
?ep PDBo:reference_to_entity ?ent .
?ent PDBo:entity.pdbx_description ?eplab .
# names may be hard to parse because they can contaion just about any character
} .
# selected.up goes here
?structref PDBo:link_to_uniprot %s
}
group by ?pdbid ?cgph'

PDBo.val <- '<http://rdf.wwpdb.org/schema/pdbx-v40.owl#>'

my.query <-
  sprintf(pdb2strand.template, PDBo.val, selected.up)

# cat(my.query)

my.endpoint <- "https://integbio.jp/rdf/sparql"

system.time(up.40.frame <-
              sparql.remote(
                endpoint = my.endpoint,
                sparql = my.query,
                jena = TRUE
              ))


my.query <-
  sprintf(pdb.ph.het.template, PDBo.val, selected.up)

# cat(my.query)

my.endpoint <- "https://integbio.jp/rdf/sparql"

system.time(
  features.40.frame <-
    sparql.remote(
      endpoint = my.endpoint,
      sparql = my.query,
      jena = TRUE
    )
)

merged.40 <- merge(x = up.40.frame,
                   y = features.40.frame,
                   by = 'pdbid',
                   all = TRUE)

merged.40[] <- lapply(merged.40[], as.character)
merged.40$cgph <- as.numeric(merged.40$cgph)


PDBo.val <- '<http://rdf.wwpdb.org/schema/pdbx-v42.owl#>'

my.query <-
  sprintf(pdb2strand.template, PDBo.val, selected.up)

my.endpoint <- "https://integbio.jp/rdf/sparql"

system.time(up.42.frame <-
              sparql.remote(
                endpoint = my.endpoint,
                sparql = my.query,
                jena = TRUE
              ))

my.query <-
  sprintf(pdb.ph.het.template, PDBo.val, selected.up)

# cat(my.query)

my.endpoint <- "https://integbio.jp/rdf/sparql"

system.time(
  features.42.frame <-
    sparql.remote(
      endpoint = my.endpoint,
      sparql = my.query,
      jena = TRUE
    )
)

merged.42 <- merge(x = up.42.frame,
                   y = features.42.frame,
                   by = 'pdbid',
                   all = TRUE)

merged.42[] <- lapply(merged.42[], as.character)
merged.42$cgph <- as.numeric(merged.42$cgph)



combo.res <-
  unique(rbind.data.frame(merged.40, merged.42, stringsAsFactors = FALSE))

desired.pdbs <- sort(unique(combo.res$pdbid))

atom.dat <- lapply(desired.pdbs, function(current.pdb) {
  # current.pdb <- "1KFF"
  print(current.pdb)
  my.strands <- combo.res$strand[combo.res$pdbid == current.pdb]
  my.strands <- unlist(strsplit(my.strands, ','))
  # print(my.strands)
  my.structure <-
    read.pdb(current.pdb, rm.alt = FALSE, rm.insert = FALSE)
  my.atoms <- my.structure$atom
  # print(unique(my.atoms$chain))
  my.alphas <- my.atoms[my.atoms$type == 'ATOM' &
                          my.atoms$elety == 'CA' &
                          my.atoms$chain %in% my.strands &
                          is.na(my.atoms$alt) &
                          is.na(my.atoms$insert), ]
  # print(unique(my.alphas$chain))
  pdbXstrand <- split(my.alphas, my.alphas$chain)
  temp <- names(pdbXstrand)
  # print(temp)
  temp <- paste0(current.pdb, "-", temp)
  # print(temp)
  names(pdbXstrand) <- temp
  # print(str(pdbXstrand))
  return(pdbXstrand)
  # # my.atoms$pdb <- current.pdb
  # # return(my.atoms)
})

atom.dat <- unlist(atom.dat, recursive = FALSE)