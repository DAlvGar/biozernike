import os
from rdkit import Chem
import simplejson
import psycopg2
from psycopg2 import sql
import time

##### SET VARIABLES #######
folder_path = 'splits'
system = 'hs90a'
moltablename = 'hs90a_mol'
sdf_extension='.sdf'
###########################

def extract_molecule_parts(molecule_name):
    parts = molecule_name.split('-')
    molid = parts[0]
    add = 0
    activity_parts = parts[1].split('.')
    activity = activity_parts[0]
    active = True if 'act' in activity else False
    if '_' in activity:
        # one of those shitty cases, add 100 to the conformer number, 
        # tricky trick i used also during import to not clash conformer names
        add = 100
    conformer = int(activity_parts[1]) + add
    return molid, active, conformer
    
# Function to create a new table for molecule data
def create_molecule_table(system, moltablename, cursor):
    create_table_query = f"""
        CREATE TABLE IF NOT EXISTS {moltablename} (
            id SERIAL PRIMARY KEY,
            name VARCHAR,
            molid VARCHAR,
            active BOOL,
            conformer INT,
            properties JSONB,
            ctab TEXT,
            UNIQUE (molid, conformer)
        );
    """
    cursor.execute(create_table_query)

# Function to create a lookup table
def create_molecule_lookup_table(system, moltablename, cursor):
    create_table_query = f"""
        CREATE TABLE IF NOT EXISTS {system}_molecule_lookup (
            id SERIAL PRIMARY KEY,
            molecule_id INTEGER REFERENCES {moltablename}(id),
            volume_id INTEGER REFERENCES {system}(id),
            UNIQUE (molecule_id, volume_id)
        );
    """
    cursor.execute(create_table_query)

# Function to insert molecule data into the database
def insert_molecule_data(cursor, system, moltablename, name,  molid, active, conformer, properties, ctab):
    insert_query_molecule = sql.SQL(f"""
        INSERT INTO {moltablename} (name, molid, active, conformer, properties, ctab)
        VALUES (%s, %s, %s, %s, %s, %s)
        ON CONFLICT (molid, conformer) DO NOTHING
        RETURNING id;
    """)
    cursor.execute(insert_query_molecule, (name, molid, active, conformer, properties, ctab))
    molecule_id = cursor.fetchone()
    return molecule_id[0] if molecule_id else None


# Function to insert molecule data into the database and populate the lookup table
def insert_molecule_data_link(cursor, system, moltablename, name,  molid, activity, conformer, properties, ctab, volume_ids):
    molecule_id = insert_molecule_data(cursor, system, moltablename, name,  molid, activity, conformer, properties, ctab)
    if not molecule_id: return 
        
    # Insert into the lookup table
    insert_query_lookup = sql.SQL(f"""
        INSERT INTO {system}_molecule_lookup (molecule_id, volume_id)
        VALUES (%s, unnest(%s))
        ON CONFLICT (molecule_id, volume_id) DO NOTHING;
    """)
    cursor.execute(insert_query_lookup, (molecule_id, volume_ids))
    return molecule_id

# Function to get volume primary ids based on molid and conformer
def get_volume_ids(cursor, system, molid, conformer):
    query = sql.SQL(f"""
        SELECT id FROM {system}
        WHERE molid = %s AND conformer = %s;
    """)
    cursor.execute(query, (molid, conformer))
    result = cursor.fetchall()
    return result

# Function to read SDF files and store data in the database
def process_sdf_files(folder_path, system, moltablename, conn):
    with conn.cursor() as cursor:
        create_molecule_table(system, moltablename, cursor)
        create_molecule_lookup_table(system, moltablename, cursor)
        conn.commit()
        start_time = time.time()
        for filename in os.listdir(folder_path):
            if filename.endswith(sdf_extension):
                sdf_file_path = os.path.join(folder_path, filename)
                print(f"Processing file {sdf_file_path}")
                file_time = time.time()
                suppl = Chem.SDMolSupplier(sdf_file_path, removeHs=False)
                for mol in suppl:
                    if mol is not None:
                        name = mol.GetProp("_Name")  # Assumes the molecule name is stored in the _Name property
                        molid, active, conformer = extract_molecule_parts(name)
                        properties = simplejson.dumps(mol.GetPropsAsDict(), ignore_nan=True)
                        ctab = Chem.MolToMolBlock(mol)

                        # Get original file primary keys associated to molid+conformer pair
                        volume_ids = get_volume_ids(cursor, system, molid, conformer)
                        
                        if volume_ids is not None and len(volume_ids) > 0:
                            # Insert data into the database and populate the lookup table
                            molecule_id = insert_molecule_data_link(cursor, system, moltablename, name,  
                                                                    molid, active, conformer, properties, ctab, volume_ids)
                            #print(f"Inserted molecule '{name}' with ID {molecule_id} and linked to volume_ids {volume_ids}")
                        else:
                            # INSERT WITHOUT LINKING
                            molecule_id = insert_molecule_data(cursor, system, moltablename, name,  molid, active, conformer, properties, ctab)
                            #print(f"Inserted molecule '{name}' with ID {molecule_id}. NOT LINKED.")
                            
                conn.commit()
                file_end_time = time.time()
                print(f"DONE in {file_end_time - file_time:.2f} seconds with {sdf_file_path}")
                
        end_time = time.time()
        print(f"Processing time: {end_time - start_time:.2f} seconds")


# Connect to the PostgreSQL database
db_params = {
    'dbname': 'dalvarez',
    'user': 'dalvarez',
    'password': 'psqldalvarez',
    'host': 'localhost',
    'port': 5432
}
conn = psycopg2.connect(**db_params)

# Process SDF files and store data in the database
process_sdf_files(folder_path, system, moltablename, conn)

# Close the database connection
conn.close()
