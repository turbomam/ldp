Sys.setenv(JAVA_HOME = 'C:\\Program Files\\Java\\jre1.8.0_111')

library(rrdf)
library(reshape)

my.query <-
  'PREFIX PDBo: <http://rdf.wwpdb.org/schema/pdbx-v40.owl#>
PREFIX dcterms:<http://purl.org/dc/terms/>
PREFIX rdf:<http://www.w3.org/1999/02/22-rdf-syntax-ns#>
PREFIX xsd:<http://www.w3.org/2001/XMLSchema#>
SELECT ?up ?dbcode (count( distinct ?pdbid) as ?pdbct) (min(?phval) as ?phmin) (max(?phval) as ?phmax)
(count( distinct ?compid) as ?nonpolyct)
FROM <http://rdf.integbio.jp/dataset/pdbj>
WHERE {
?sr PDBo:link_to_uniprot ?up .
?sr PDBo:struct_ref.db_code ?dbcode .
?sr PDBo:of_datablock ?pdb .
?pdb dcterms:identifier ?pdbid .
optional { ?cg rdf:type PDBo:exptl_crystal_grow .
?cg PDBo:of_datablock ?pdb .
?cg PDBo:exptl_crystal_grow.pH ?cgph .
bind(xsd:double(?cgph) as ?phval) .
} .
optional {?enp rdf:type	PDBo:pdbx_entity_nonpoly .
?enp PDBo:of_datablock ?pdb .
?enp PDBo:pdbx_entity_nonpoly.comp_id ?compid
}
}
group by ?up ?dbcode
having (count( distinct ?pdbid) > 1 )
# limit 99
'
my.endpoint <- "https://integbio.jp/rdf/sparql"

system.time(frequent.ups <-
              sparql.remote(
                endpoint = my.endpoint,
                sparql = my.query,
                jena = TRUE
              ))

frequent.ups <- as.data.frame(frequent.ups)
frequent.ups[] <- lapply(frequent.ups[], as.character)
frequent.ups[, 3:6] <- lapply(frequent.ups[, 3:6], as.numeric)
frequent.ups$phrange <- frequent.ups$phmax - frequent.ups$phmin

###

# the descriptions for pdb entries may not be consistent,
# even if they have the same exclusive protein, by uniprot identifier
# so get the uniprot description form uniprot
# for proteins that are know to ahve a pdb entry
# would be nice to federate, but has been timing out for me

my.endpoint <- "http://sparql.uniprot.org/sparql"

my.query <- "PREFIX up:<http://purl.uniprot.org/core/>
PREFIX rdfs:<http://www.w3.org/2000/01/rdf-schema#>
PREFIX skos:<http://www.w3.org/2004/02/skos/core#>
SELECT distinct ?protein ?name
WHERE
{
  ?protein a up:Protein .
  ?protein rdfs:seeAlso ?db .
  ?db up:database <http://purl.uniprot.org/database/PDB> .
  ?protein up:recommendedName ?recommended .
  ?recommended up:fullName ?name .
  # ?protein up:encodedBy ?gene .
  # ?gene skos:prefLabel ?text .
}"

  system.time(up2pdb <-
                sparql.remote(
                  endpoint = my.endpoint,
                  sparql = my.query,
                  jena = TRUE
                ))
  
  up2pdb <- as.data.frame(up2pdb)
  
  # merge, since we're not federating (yet)
  
  pdb.ct.up.desc <- merge(
    x = frequent.ups,
    y = up2pdb,
    by.x = "up",
    by.y = "protein",
    all.x = TRUE
  )
  
  save(pdb.ct.up.desc, file = "pdb_ct_up_desc.Rdata")
  