from astropy.table import Table, vstack, join
from astroquery.gaia import Gaia
from astroquery.simbad import Simbad
import numpy as np
import time

def get_dr3_sources_from_simbad(identifiers):
    """
    Query SIMBAD for a list of identifiers and return the 
    Gaia DR3 source IDs according to SIMBAD.

    Parameters:
    -----------
    identifiers : np.array
        List of SIMBAD-compatible identifier strings
    """
    simbad = Simbad()
    simbad.add_votable_fields("ids")
    result_table = simbad.query_objects(identifiers)
    dr3_ids = np.zeros(len(identifiers), dtype=np.int64)
    for i,row in enumerate(result_table):
        all_ids = row['ids'].split('|')
        dr3_idx = np.flatnonzero(np.char.find(all_ids, 'Gaia DR3') != -1)
        if len(dr3_idx == 1):
            dr3_ids[i] = int(all_ids[dr3_idx[0]].split(' ')[2])
    return dr3_ids  
    

def query_gaia_astrophysical_parameters(table, source_id_name='source_id', batch_size=500):
    """
    Query Gaia DR3 astrophysical parameters for a list of source_ids and
    add them to the provided Astropy Table in-place.

    Parameters:
    -----------
    table : astropy.table.Table
        Input table containing a column 'source_id'
    batch_size : int
        Number of sources to query in each batch
    """
    source_ids = table[source_id_name]
    all_results = []

    for i in range(0, len(source_ids), batch_size):
        batch = source_ids[i:i+batch_size]
        id_list_str = ",".join(str(sid) for sid in batch)
        query = f"""
        SELECT *
        FROM gaiadr3.astrophysical_parameters
        WHERE source_id IN ({id_list_str})
        """
        try:
            job = Gaia.launch_job_async(query, dump_to_file=False)
            batch_result = job.get_results()
            all_results.append(batch_result)
        except Exception as e:
            print(f"Query failed for batch {i//batch_size + 1}: {e}")
            continue

        time.sleep(1)

    if all_results:
        gaia_results = vstack(all_results)

        # Join results into the input table based on source_id
        updated_table = join(table, gaia_results, keys_left=source_id_name, keys_right="source_id", join_type="left")

        return updated_table
    else:
        print("No results returned from Gaia archive.")

# Example usage
if __name__ == "__main__":
    # Dummy example table
    input_table = Table()
    input_table["source_id"] = [587735153956000000, 587735153956000001, 587735153956000002]

    query_gaia_astrophysical_parameters(input_table)

    print(input_table)
