select gs.dsstox_substance_id, c.acd_iupac_name, c.jchem_inchi_key, c.indigo_inchi_key
from generic_substances gs
left join generic_substance_compounds gsc on gsc.fk_generic_substance_id = gs.id  
left join compounds c on c.id = gsc.fk_compound_id
where c.jchem_inchi_key like 'RFDAIACWWDREDC%';

