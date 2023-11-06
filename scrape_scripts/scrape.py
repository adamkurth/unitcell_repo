from rcsbsearchapi.search import TextQuery, AttributeQuery, SequenceQuery, SeqMotifQuery, StructSimilarityQuery, Attr
from rcsbsearchapi import rcsb_attributes as attrs
from rcsbsearchapi.const import STRUCTURE_ATTRIBUTE_SEARCH_SERVICE

# From rcsbsearchapi 
# comment and add as needed, for reference check the quickstart.ipynb in rcsbsearchapi repo 

def get_pdb_ids(space_group, str_type, limit=50):
    # By default, service is set to "text" for structural attribute search
    q1 = AttributeQuery("symmetry.cell_setting", "exact_match", str_type, STRUCTURE_ATTRIBUTE_SEARCH_SERVICE)
    q2 = AttributeQuery("symmetry.space_group_name_H_M", "exact_match", space_group, STRUCTURE_ATTRIBUTE_SEARCH_SERVICE)
    query = q2  # combining queries use & | operators
    pdb_ids = list(query())
    pdb_lim = truncate(pdb_ids, limit)
    print(f"From search: {space_group} \n", pdb_lim)
    return pdb_lim

def truncate(pdb_ids, limit):
    if len(pdb_ids) > limit:
        pdb_ids = pdb_ids[:limit]
    return pdb_ids

# get_pdb_ids("P 21 21 21", "triclinic", limit=10)