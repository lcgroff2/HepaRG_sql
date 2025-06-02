library(tidyverse)
library(httr)
library(jsonlite)
library(png)
library(magick)
library(digest)
library(RMariaDB)

setwd("~/Data/VSSTox")

#Need Python version of rdkit because R version is rubbish.
library(reticulate)
# install_miniconda(update = T)
use_condaenv(condaenv = "r-reticulate")
# py_install("rdkit")

rdkit <- import("rdkit.Chem")
rdMolDescriptors <- import("rdkit.Chem.rdMolDescriptors")
Descriptors <- import("rdkit.Chem.Descriptors")

#Set the structures to be enumerated by pulling directly from DSSTox
#Connect to database
dbDisconnect(DSSTox)


#VSSTox Temp

dbDisconnect(VSSTox)

dbListTables(VSSTox)

checksumFx <- function(index) {
  split <- str_split(index,"")[[1]]
  checksum <- 0
  for(i in 1:length(split)) {
    checksum <- checksum + as.numeric(split[i])*i
    
  }
  checksum <- checksum%%10
  return(checksum)
}

sid_from_inchikey <- function(inchikey) {
  #try query on the inchikey derived from SMILES result of in silico tool:
  new_sid <- dbGetQuery(DSSTox,paste0("select gs.dsstox_substance_id from generic_substances gs left join generic_substance_compounds gsc on gsc.fk_generic_substance_id = gs.id left join compounds c on c.id = gsc.fk_compound_id where c.jchem_inchi_key = '",inchikey,"'"))
  if (nrow(new_sid) == 1){
    return(new_sid)
  } else if (nrow(new_sid) == 0){
    #try query on <inchikey first block>-UHFFFAOYSA-N:
    new_sid <- dbGetQuery(DSSTox,paste0("select gs.dsstox_substance_id from generic_substances gs left join generic_substance_compounds gsc on gsc.fk_generic_substance_id = gs.id left join compounds c on c.id = gsc.fk_compound_id where c.jchem_inchi_key = '",str_split(inchikey,"-")[[1]][1],"-UHFFFAOYSA-N'"))}
  if (nrow(new_sid) == 1){
    return(new_sid)
  } else if (nrow(new_sid) == 0){
    #try query on <inchikey first block>-UHFFFAOYNA-N:
    new_sid <- dbGetQuery(DSSTox,paste0("select gs.dsstox_substance_id from generic_substances gs left join generic_substance_compounds gsc on gsc.fk_generic_substance_id = gs.id left join compounds c on c.id = gsc.fk_compound_id where c.jchem_inchi_key = '",str_split(inchikey,"-")[[1]][1],"-UHFFFAOYNA-N'"))
    return(new_sid)
  }
}

fpath <- '~/Data/VSSTox/heparg/CTS ChemAxon Phase I'
metabolite_data <- fromJSON(txt = paste0(fpath,'/heparg_cts_human_3gens_9_complete.json'))
heparg <- read_xlsx(paste0(getwd(),'/heparg/heparg-chems-dashboard.xlsx'))
heparg <- heparg[which(!is.na(heparg$QSAR_READY_SMILES)),]
#for CTS, drop rows where API call failed/returned no result
metabolite_data <- metabolite_data[which(metabolite_data$api_success),]
#determine if any sids are NA before processing through for loop, replace those
metab_sids_na <- metabolite_data[is.na(metabolite_data$input$dtxsid),]
if (nrow(metab_sids_na) > 0) {
  for (i in 1:nrow(metab_sids_na)) {
    #pull DTXSID from prod_dsstox by matching on inchikey
    new_sid <- sid_from_inchikey(metab_sids_na[i,]$input$inchikey)
    if (nrow(new_sid) == 1){
      metabolite_data[(metabolite_data$input$inchikey == metab_sids_na[i,]$input$inchikey), ]$input$dtxsid <- new_sid
    }
  }
  #pull new_sid from heparg data if prod_dsstox query doesn't yield new_sid
  metab_sids_na <- metabolite_data[is.na(metabolite_data$input$dtxsid),]
  if (nrow(metab_sids_na) > 0){
    for (s in 1:length(metab_sids_na$input$smiles)){
      match_sid <- heparg[heparg$QSAR_READY_SMILES == metab_sids_na$input$smiles[[s]],]$DTXSID
      metabolite_data$input[metabolite_data$input$smiles == metab_sids_na$input$smiles[[s]],]$dtxsid <- match_sid
    }
  metab_sids_na <- metabolite_data[is.na(metabolite_data$input$dtxsid),]} 
  if (nrow(metab_sids_na) == 0) remove(metab_sids_na)
}


#cross-check that all NA cases were taken care of, otherwise drop them.
metabolite_data <- metabolite_data[which(!is.na(metabolite_data$input$dtxsid)),]
sids <- metabolite_data$input$dtxsid

for (sid in sids) {
  index <- dbGetQuery(VSSTox,paste0("select count(ps.id) from parent_structures ps"))[[1]][1] + 1
  parent_structures <- data.frame("id" = index,
                                  "vsstox_pid" = NA,
                                  "dsstox_substance_id" = NA,
                                  "dsstox_compound_id" = NA,
                                  "mrv_file" = NA,
                                  "mol_file" = NA,
                                  "smiles" = NA,
                                  "indigo_inchi" = NA,
                                  "indigo_inchikey" = NA,
                                  "jchem_inchi" = NA,
                                  "jchem_inchikey" = NA,
                                  "acd_iupac_name" = NA,
                                  "acd_index_name" = NA,
                                  "mol_formula" = NA,
                                  "mol_weight" = NA,
                                  "monoisotopic_mass" = NA,
                                  "has_stereo" = NA,
                                  "created_at" = NA,
                                  "updated_at" = NA,
                                  "created_by" = NA,
                                  "updated_by" = NA)
  check <- dbGetQuery(VSSTox,paste0("select dsstox_substance_id from parent_structures where dsstox_substance_id = '",sid,"'"))
  if (nrow(check) > 0) next #Don't insert the same DTXSID more than once....
  #Should add INCHIKEY check, too...
  parent_structures[1,"vsstox_pid"] <- paste0("VTXPID",checksumFx(index),"0",index)
  parent_structures[1,"dsstox_substance_id"] <- sid
  # url <- paste0("https://api-ccte.epa.gov/chemical/detail/search/by-dtxsid/",sid)
  # response <- GET(url, add_headers('x-api-key' = key), content_type("application/octet-stream"),flatten = T)
  # details <- content(response)
  details <- dbGetQuery(DSSTox,paste0("select gs.dsstox_substance_id, c.dsstox_compound_id, c.mrv_file, c.mol_file, c.smiles, c.inchi as 'jchem_inchi', c.jchem_inchi_key, c.indigo_inchi, c.indigo_inchi_key, c.acd_iupac_name,c.acd_index_name, c.mol_formula, c.mol_weight, c.monoisotopic_mass, c.has_stereochemistry as 'has_stereo' from generic_substances gs left join generic_substance_compounds gsc on gs.id = gsc.fk_generic_substance_id left join compounds c on gsc.fk_compound_id = c.id where gs.dsstox_substance_id = '",sid,"'"))
  if (is.null(details$indigo_inchi_key)){
    new_sid <- NULL
    cur_inchikey <- metabolite_data[metabolite_data$input$dtxsid == sid,]$input$inchikey
    new_sid <- sid_from_inchikey(cur_inchikey)}
    if (!is.null(new_sid[[1]]) | nrow(new_sid) == 0) next #skip if no updated sid found
  check <- dbGetQuery(VSSTox,paste0("select dsstox_substance_id from parent_structures where dsstox_substance_id = '",new_sid[[1]],"'"))
  if (nrow(check) > 0) next #Don't insert the same DTXSID more than once....
  parent_structures[1,"dsstox_substance_id"] <- new_sid[[1]]
  details <- dbGetQuery(DSSTox,paste0("select gs.dsstox_substance_id, c.dsstox_compound_id, c.mrv_file, c.mol_file, c.smiles, c.inchi as 'jchem_inchi', c.jchem_inchi_key, c.indigo_inchi, c.indigo_inchi_key, c.acd_iupac_name,c.acd_index_name, c.mol_formula, c.mol_weight, c.monoisotopic_mass, c.has_stereochemistry as 'has_stereo' from generic_substances gs left join generic_substance_compounds gsc on gs.id = gsc.fk_generic_substance_id left join compounds c on gsc.fk_compound_id = c.id where gs.dsstox_substance_id = '",new_sid[[1]],"'"))
  #skip to next sid if indigo_inchikey is still NULL after attempt to fill gap from metabolism_data$input$inchikey
  if (is.null(details$indigo_inchi_key)) next
  parent_structures[1,"dsstox_compound_id"] <- details$dsstox_compound_id
  if (is.null(details$smiles) == F) parent_structures[1,"smiles"] <- details$smiles
  if (is.null(details$acd_iupac_name) == F) parent_structures[1,"acd_iupac_name"] <- details$acd_iupac_name
  if (is.null(details$acd_index_name) == F) parent_structures[1,"acd_index_name"] <- details$acd_index_name
  if (is.null(details$has_stereo) == F) parent_structures[1,"has_stereo"] <- details$has_stereo
  if (is.null(details$indigo_inchi) == F) parent_structures[1,"indigo_inchi"] <- details$indigo_inchi
  if (is.null(details$indigo_inchi_key) == F) parent_structures[1,"indigo_inchikey"] <- details$indigo_inchi_key
  if (is.null(details$jchem_inchi) == F) parent_structures[1,"jchem_inchi"] <- details$jchem_inchi
  if (is.null(details$jchem_inchi_key) == F) parent_structures[1,"jchem_inchikey"] <- details$jchem_inchi_key
  if (is.null(details$mol_file) == F) parent_structures[1,"mol_file"] <- details$mol_file
  if (is.null(details$mrv_file) == F) parent_structures[1,"mrv_file"] <- details$mrv_file
  if (is.null(details$mol_formula) == F) parent_structures[1,"mol_formula"] <- details$mol_formula
  if (is.null(details$mol_weight) == F) parent_structures[1,"mol_weight"] <- details$mol_weight
  if (is.null(details$monoisotopic_mass) == F) parent_structures[1,"monoisotopic_mass"] <- details$monoisotopic_mass
  parent_structures[1,"created_at"] <- format(Sys.time(), '%Y-%m-%d %H:%M:%S')
  parent_structures[1,"updated_at"] <- format(Sys.time(), '%Y-%m-%d %H:%M:%S')
  parent_structures[1,"created_by"] <- "lgroff"
  parent_structures[1,"updated_by"] <- "lgroff"
  
  dbExecute(VSSTox,"start transaction;")
  dbExecute(VSSTox,"drop table if exists myTempTable")
  dbWriteTable(VSSTox,"myTempTable",parent_structures,row.names = FALSE)
  dbExecute(VSSTox,"insert into parent_structures select * from myTempTable")
  dbExecute(VSSTox,"drop table if exists myTempTable")
  dbExecute(VSSTox,"commit;")
}
details$dsstox_substance_id
#Add relationship types to parent_compound_relationship_types
# index <- 1
index <- dbGetQuery(VSSTox,paste0("select count(pcrt.id) from parent_compound_relationship_types pcrt"))[[1]][1] + 1
#Don't create multiples of the same thing...
if (index < 4) {
  parent_compound_relationship_types <- data.frame("id" = NA,
                                                   "name" = NA,
                                                   "description" = NA,
                                                   "metab_software_name" = NA,
                                                   "metab_software_version" = NA,
                                                   "metab_depth" = NA,
                                                   "metab_organism" = NA,
                                                   "metab_site_of_metabolism" = NA,
                                                   "metab_model" = NA,
                                                   # "metab_enzyme" = NA,
                                                   # "metab_mechanism" = NA,
                                                   # "metab_generation" = NA,
                                                   "created_at" = NA,
                                                   "updated_at" = NA,
                                                   "created_by" = NA,
                                                   "updated_by" = NA)
  parent_compound_relationship_types[1,"id"] <- as.integer(index)
  parent_compound_relationship_types[1,"name"] <- "Metabolism via EPA's Chemical Transformation Simulator (CTS)"
  parent_compound_relationship_types[1,"description"] <- "This relationship denotes the connection between a Parent structure and its metabolites."
  parent_compound_relationship_types[1,"metab_software_name"] <- metabolite_data$software[[1]][1]
  parent_compound_relationship_types[1,"metab_software_version"] <- metabolite_data$version[[1]][1]
  parent_compound_relationship_types[1,"metab_depth"] <- metabolite_data$params$depth[[1]][1]
  parent_compound_relationship_types[1,"metab_organism"] <- metabolite_data$params$organism[[1]][1]
  parent_compound_relationship_types[1,"metab_site_of_metabolism"] <- metabolite_data$params$site_of_metabolism[[1]][1]
  parent_compound_relationship_types[1,"metab_model"] <- paste(metabolite_data$params$model[[1]][1])
  # parent_compound_relationship_types[1,"metab_enzyme"] <- metabolite_data$output
  # parent_compound_relationship_types[1,"metab_mechanism"] <- metabolite_data$params
  # parent_compound_relationship_types[1,"metab_generation"] <- metabolite_data$params
  parent_compound_relationship_types[1,"created_at"] <- format(Sys.time(), '%Y-%m-%d %H:%M:%S')
  parent_compound_relationship_types[1,"updated_at"] <- format(Sys.time(), '%Y-%m-%d %H:%M:%S')
  parent_compound_relationship_types[1,"created_by"] <- "lgroff"
  parent_compound_relationship_types[1,"updated_by"] <- "lgroff"

  dbExecute(VSSTox,"start transaction;")
  dbExecute(VSSTox,"drop table if exists myTempTable")
  dbWriteTable(VSSTox,"myTempTable",parent_compound_relationship_types,row.names = FALSE)
  dbExecute(VSSTox,"insert into parent_compound_relationship_types select * from myTempTable")
  dbExecute(VSSTox,"drop table if exists myTempTable")
  dbExecute(VSSTox,"commit;")
}
i<-1
if (exists('i')){l <- i} else {l<-1}
#Load virtual compounds
for (i  in l:length(metabolite_data$output)) {
  for (j in 1:length(metabolite_data$output[[i]]$successors)) {
    index <- dbGetQuery(VSSTox,paste0("select max(c.id) from compounds c"))[[1]][1] + 1
    # check <- dbGetQuery(VSSTox,paste0("select dsstox_substance_id from markush_structures where dsstox_substance_id = '",sid,"'"))
    # if (nrow(check) > 0) next #No point in redoing existing entries... for now...

    smiles <- metabolite_data$output[[i]]$successors[[j]]$metabolite$smiles
    # temp <- nm[which(str_detect(nm,sid) == TRUE)]
    # smiles <- temp[which(str_detect(temp,"smiles") == TRUE)]
    # smiles <- fromJSON(txt = paste0("~/R Data/QACs/markush_jsons/",smiles))
    if (length(smiles) > 10000) next #The really large  sets can be done later/never...
    # smiles <- sapply(smiles, str_replace_all, pattern = "\\-\\*", replacement = "[H]")
    if (length(smiles) < 2) next #Pointless to store Markush that didn't enumerate properly.
    compounds <- data.frame("id" = NA,
                            "vsstox_vid" = NA,
                            "dsstox_compound_id" = NA,
                            "mrv_file" = NA,
                            "mol_file" = NA,
                            "smiles" = NA,
                            "indigo_inchi" = NA,
                            "indigo_inchikey" = NA,
                            "jchem_inchi" = NA,
                            "jchem_inchikey" = NA,
                            "acd_iupac_name" = NA,
                            "acd_index_name" = NA,
                            "mol_formula" = NA,
                            "mol_weight" = NA,
                            "monoisotopic_mass" = NA,
                            "has_stereo" = NA,
                            "created_at" = NA,
                            "updated_at" = NA,
                            "created_by" = NA,
                            "updated_by" = NA)
    compounds <- compounds[which(is.na(compounds$id)==F),]
    for (k in 1:length(smiles)) {
      # molecule <- rdkit$MolFromSmiles(smiles[i],sanitize=TRUE)
      tryCatch(molecule <- rdkit$MolFromSmiles(smiles[k],sanitize=TRUE), error = function(e) {molecule <<- NA})
      if (is.null(molecule)) next
      indigo_inchikey <- rdkit$MolToInchiKey(molecule)
      if (indigo_inchikey == "") next
      cid <- dbGetQuery(DSSTox,paste0("select c.dsstox_compound_id from compounds c where c.jchem_inchi_key = '",indigo_inchikey,"'"))
      check <- dbGetQuery(VSSTox,paste0("select c.indigo_inchikey from compounds c where c.indigo_inchikey = '",indigo_inchikey,"'"))
      # if (nrow(check) > 0) index <- index + 1
      if (nrow(check) > 0) next #No point in redoing existing entries... for now...
      if (indigo_inchikey %in% compounds$indigo_inchikey) next #Unique SMILES does not mean unique InChIKey...
      if (nrow(cid) > 0) compounds[k,"dsstox_compound_id"] <- cid
      # compounds[i,"id"] <- index

      compounds[k,"smiles"] <- smiles[k]
      compounds[k,"indigo_inchi"] <- rdkit$MolToInchi(molecule)
      compounds[k,"indigo_inchikey"] <- rdkit$MolToInchiKey(molecule)
      compounds[k,"mol_file"] <- rdkit$MolToV3KMolBlock(molecule)
      compounds[k,"mrv_file"] <- rdkit$MolToMrvBlock(molecule)
      compounds[k,"mol_formula"] <- rdMolDescriptors$CalcMolFormula(molecule)
      compounds[k,"mol_weight"] <- Descriptors$MolWt(molecule)
      compounds[k,"monoisotopic_mass"] <- Descriptors$ExactMolWt(molecule)
      compounds[k,"created_at"] <- format(Sys.time(), '%Y-%m-%d %H:%M:%S')
      compounds[k,"updated_at"] <- format(Sys.time(), '%Y-%m-%d %H:%M:%S')
      compounds[k,"created_by"] <- "lgroff"
      compounds[k,"updated_by"] <- "lgroff"
      # index <- index + 1


    }

    if (nrow(compounds) > 0) {
      compounds <- compounds[which(is.na(compounds$indigo_inchikey)==F),]
      if (nrow(compounds) == 1) {
        compounds$id <- index
        compounds$vsstox_vid <- paste0("VTXVID",checksumFx(index),"0",index)
      }
      if (nrow(compounds) > 1) {
        compounds$id <- seq(from = index,to = (index + nrow(compounds) - 1))
        compounds$vsstox_vid <- paste0("VTXVID",checksumFx(index),"0",seq(from = index,to = (index + nrow(compounds) - 1)))

      }
      dbExecute(VSSTox,"start transaction;")
      dbExecute(VSSTox,"drop table if exists myTempTable")
      dbWriteTable(VSSTox,"myTempTable",compounds,row.names = FALSE)
      dbExecute(VSSTox,"insert into compounds select * from myTempTable")
      dbExecute(VSSTox,"drop table if exists myTempTable")
      dbExecute(VSSTox,"commit;")
    }
    
    #parent_compound_relationships
    #add case where index is zero...
    # relationship_index <- 1
    relationship_index <- dbGetQuery(VSSTox,paste0("select max(pcr.id) from parent_compound_relationships pcr"))[[1]][1] + 1
    parent_index <- dbGetQuery(VSSTox,paste0("select ps.id from parent_structures ps where ps.dsstox_substance_id = '",sids[i],"'"))[[1]][1]
    if (is.na(parent_index)) 
      new_sid <- sid_from_inchikey(metabolite_data$input[i,]$inchikey) #NA parent index will cause assignment error
    parent_index <- dbGetQuery(VSSTox,paste0("select ps.id from parent_structures ps where ps.dsstox_substance_id = '",new_sid[[1]],"'"))[[1]][1]
    if (is.na(parent_index)) next #skip if parent index is still NA
    parent_compound_relationships <- data.frame("id" = NA,
                                                "fk_parent_id" = NA,
                                                "fk_compound_id" = NA,
                                                "fk_relationship_type" = NA,
                                                "metab_enzyme" = NA,
                                                "metab_mechanism" = NA,
                                                "metab_generation" = NA,
                                                "created_at" = NA,
                                                "updated_at" = NA,
                                                "created_by" = NA,
                                                "updated_by" = NA)
    parent_compound_relationships <- parent_compound_relationships[which(is.na(parent_compound_relationships$fk_compound_id)==F),]
    for (k in 1:length(smiles)) {
      # molecule <- rdkit$MolFromSmiles(smiles[i],sanitize=TRUE)
      tryCatch(molecule <- rdkit$MolFromSmiles(smiles[k],sanitize=TRUE), error = function(e) {molecule <<- NA})
      if (is.null(molecule)) next
      # parent_compound_relationships[i,"id"] <- relationship_index
      # relationship_index <- relationship_index + 1
      parent_compound_relationships[k,"fk_parent_id"] <- parent_index
      fk_compound_id <- dbGetQuery(VSSTox,paste0("select c.id from compounds c where c.indigo_inchikey = '",rdkit$MolToInchiKey(molecule),"'"))
      parent_compound_relationships[k,"fk_compound_id"] <- fk_compound_id
      parent_compound_relationships[k,"fk_relationship_type"] <- 3
      parent_compound_relationships[k,"metab_enzyme"] <- paste(metabolite_data$output[[i]]$successors[[j]]$enzyme[k])
      parent_compound_relationships[k,"metab_mechanism"] <- metabolite_data$output[[i]]$successors[[j]]$mechanism[k]
      parent_compound_relationships[k,"metab_generation"] <- metabolite_data$output[[i]]$successors[[j]]$generation[k]
      parent_compound_relationships[k,"created_at"] <- format(Sys.time(), '%Y-%m-%d %H:%M:%S')
      parent_compound_relationships[k,"updated_at"] <- format(Sys.time(), '%Y-%m-%d %H:%M:%S')
      parent_compound_relationships[k,"created_by"] <- "lgroff"
      parent_compound_relationships[k,"updated_by"] <- "lgroff"
    }
    if (nrow(parent_compound_relationships) > 0) {
      parent_compound_relationships <- parent_compound_relationships[which(is.na(parent_compound_relationships$fk_parent_id)==F),]
      parent_compound_relationships$id <- seq(from = relationship_index,to = (relationship_index + nrow(parent_compound_relationships) - 1))
      dbExecute(VSSTox,"start transaction;")
      dbExecute(VSSTox,"drop table if exists myTempTable")
      dbWriteTable(VSSTox,"myTempTable",parent_compound_relationships,row.names = FALSE)
      dbExecute(VSSTox,"insert into parent_compound_relationships select * from myTempTable")
      dbExecute(VSSTox,"drop table if exists myTempTable")
      dbExecute(VSSTox,"commit;")
    }
  }
}
