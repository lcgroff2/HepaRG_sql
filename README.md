# HepaRG_sql
SQL scripting in R for insertion of MetSim metabolic predictions into MySQL database.

* [Full R Script](https://github.com/lcgroff2/HepaRG_sql/blob/main/parents_load_script_lgroff_redacted.R) for conversion of JSON formatted unstructured MongoDB datasets of predicted metabolites into MySQL relational database format, followed by insertion into database.
* Example SQL Test Query Script for retrieving Distributed Structure-Searchable Toxicity (DSSTox) Database Substance Identifier (DTXSID, sid in code) based on InChIKey structural identifiers with commands from MariaDB R Package
```r
sid_from_inchikey <- function(inchikey) {
  #try MySQL query on the inchikey derived from SMILES result of in silico tool:
  new_sid <- dbGetQuery(DSSTox,paste0("select gs.dsstox_substance_id
                                       from generic_substances gs
                                       left join generic_substance_compounds gsc on gsc.fk_generic_substance_id = gs.id
                                       left join compounds c on c.id = gsc.fk_compound_id where c.jchem_inchi_key = '",inchikey,"'"))
  if (nrow(new_sid) == 1){
    return(new_sid)
  } else if (nrow(new_sid) == 0)
{
    #try MySQL query on <inchikey first block>-UHFFFAOYSA-N (standard InChIKey mapping for no stereochemistry - neutral charge):
    new_sid <- dbGetQuery(DSSTox,paste0("select gs.dsstox_substance_id
                                         from generic_substances gs
                                         left join generic_substance_compounds gsc on gsc.fk_generic_substance_id = gs.id
                                         left join compounds c on c.id = gsc.fk_compound_id
                                         where c.jchem_inchi_key = '",str_split(inchikey,"-")[[1]][1],"-UHFFFAOYSA-N'"))
}
  if (nrow(new_sid) == 1){
    return(new_sid)
  } else if (nrow(new_sid) == 0)
{
    #try MySQL query on <inchikey first block>-UHFFFAOYNA-N (Nonstandard InChIKey mapping for no stereochemistry - neutral charge):
    new_sid <- dbGetQuery(DSSTox,paste0("select gs.dsstox_substance_id from generic_substances gs
                                         left join generic_substance_compounds gsc on gsc.fk_generic_substance_id = gs.id
                                         left join compounds c on c.id = gsc.fk_compound_id
                                         where c.jchem_inchi_key = '",str_split(inchikey,"-")[[1]][1],"-UHFFFAOYNA-N'"))
    return(new_sid)
  }

}
```
