import numpy as np
import psycopg2
from psycopg2 import sql
import pandas as pd

def get_ctab_by_mol_id_and_conformer(cursor, moltable, molid, conformer):
    query = sql.SQL(f"""
        SELECT ctab
        FROM {moltable}
        WHERE molid = %s AND conformer = %s;
    """)
    cursor.execute(query, (molid, conformer))
    result = cursor.fetchone()
    return result[0] if result else None

def get_data_by_mol_id_and_conformer(cursor, system, moltable, molid, conformer):
    query = sql.SQL(f"""
        SELECT m.ctab, b.original, b.active, b.map, b.fp, b.moments_real, b.moments_imag
        FROM {moltable} m
        JOIN {system}_molecule_lookup bm ON m.id = bm.molecule_id
        JOIN {system} b ON bm.volume_id = b.id
        WHERE m.molid = %s AND m.conformer = %s;
    """)
    cursor.execute(query, (molid, conformer))
    result = cursor.fetchall()
    if len(result):
        d = pd.DataFrame(result, columns=['ctab','original','active','map','fp','moments_real','moments_imag'])
        d['moments'] = d.apply(lambda x: np.vectorize(complex)(x.moments_real, x.moments_imag), axis=1)
        d = d.drop(columns=['moments_real','moments_imag'])
        return d
    return None

def write_complex_data_to_file(file_path, real_values, imag_values):
    # Combine real and imaginary parts to form complex numbers
    complex_numbers = np.vectorize(complex)(real_values, imag_values)

    # Write complex numbers to a text file, one per line
    with open(file_path, 'w') as file:
        for complex_number in complex_numbers:
            file.write(f"{complex_number}\n")
