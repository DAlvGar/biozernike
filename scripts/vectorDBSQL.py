import psycopg2
import numpy as np
from io import BytesIO
import pickle

class VectorDatabase:
    def __init__(self, tablename, dbname, user, password, host='localhost', port=5432):
        self.conn = psycopg2.connect(dbname=dbname, user=user, password=password, host=host, port=port)
        self.table=tablename
        self.cursor = self.conn.cursor()
        self.create_table()

    def create_table(self):
        self.cursor.execute(f'''
            CREATE EXTENSION IF NOT EXISTS vector;

            -- Create the molecules table
            CREATE TABLE IF NOT EXISTS {self.table}_mol (
                id SERIAL PRIMARY KEY,
                molid VARCHAR(150),
                activity VARCHAR(150),
                conformer int,
                properties JSONB,
                ctab TEXT,
                UNIQUE (molid, conformer)
            );

            -- Create the system table
            CREATE TABLE IF NOT EXISTS {self.table} (
                id SERIAL PRIMARY KEY,
                original VARCHAR(150),
                molid VARCHAR(150),
                active BOOLEAN,
                conformer int,
                map VARCHAR(150),
                fp vector(72),
                moments_real float[],
                moments_imag float[]
            );

            -- Create the system_molecule_lookup table
            CREATE TABLE IF NOT EXISTS {self.table}_molecule_lookup (
                id SERIAL PRIMARY KEY,
                volume_id INTEGER REFERENCES {self.table}(id),
                molecule_id INTEGER REFERENCES  {self.table}_mol(id),
                UNIQUE (volume_id, molecule_id)
            );

            CREATE INDEX IF NOT EXISTS idx_{self.table}_molid_conformer_mapname ON {self.table} (molid, conformer, map);
            CREATE INDEX IF NOT EXISTS idx_{self.table}_molid_conformer ON {self.table} (molid, conformer);
            CREATE INDEX IF NOT EXISTS idx_{self.table}_molid ON {self.table} (molid);
        ''')
        self.conn.commit()

    def exists(self, original):
        self.cursor.execute(f'SELECT 1 FROM {self.table} WHERE original = %s;', (original,))
        result = self.cursor.fetchone()
        return bool(result)
        
    def sql(self, sql, data=None, commit=False):
        self.cursor.execute(sql, data)
        rows = self.cursor.fetchall()
        if commit: self.conn.commit()
        return rows
    
    
    def insert(self, data, commit=True):
        molid, conformer = data[1], data[3]
        
        # Check if the molid and conformer combination already exists in the molecules table
        self.cursor.execute(f'SELECT id FROM {self.table}_mol WHERE molid = %s AND conformer = %s;', (molid, conformer))
        result = self.cursor.fetchone()
        
        # Insert into the vector table
        insert_query_vector = f'''
            INSERT INTO {self.table} (original, molid, active, conformer, map, fp, moments_real, moments_imag)
            VALUES (%s, %s, %s, %s, %s, %s, %s, %s)
            RETURNING id;
        '''
        self.cursor.execute(insert_query_vector, data)
        vector_id = self.cursor.fetchone()[0]
        
        if result:
            # If the entry exists in the molecules table, use the existing ID
            molecule_id = result[0]
            # Insert into the lookup table
            insert_query_lookup = f'''
                INSERT INTO {self.table}_molecule_lookup (volume_id, molecule_id)
                VALUES (%s, %s);
            '''
            self.cursor.execute(insert_query_lookup, (vector_id, molecule_id))
        # else:
        #     # If not, insert into the molecules table and retrieve the new ID
        #     insert_query_molecule = f'''
        #         INSERT INTO {self.table}_mol (molid, activity, conformer, properties, ctab)
        #         VALUES (%s, %s, %s, %s, %s)
        #         RETURNING id;
        #     '''
        #     self.cursor.execute(insert_query_molecule, (molid, None, conformer, None, None))
        #     molecule_id = self.cursor.fetchone()[0]
        if commit:
            self.conn.commit()

        return vector_id
        
    # def insert(self, data, commit=True):
    #     sql = f"INSERT INTO {self.table} (original, molid, active, conformer, map, fp, moments_real, moments_imag) \
    #                 VALUES(%s, %s, %s, %s, %s, %s, %s, %s) \
    #                 RETURNING id;"
    #     self.cursor.execute(sql, data)
    #         #rows = cursor.fetchall()
    #     if commit: self.conn.commit()
    #     #return rows
    
    # def insert_many(self, data, commit=True):
    #     sql = f"INSERT INTO {self.table}(original, molid, active, conformer, map, fp, moments_real, moments_imag) \
    #                 VALUES(%s, %s, %s, %s, %s, %s, %s, %s) \
    #                 RETURNING id;"
    #     self.cursor.executemany(sql, data)
    #         #rows = cursor.fetchall()
    #     if commit: self.conn.commit()
        #return rows
        
    def insert_many(self, data, commit=True):
        for entry in data:
            self.insert(entry, commit=False)
        if commit:
            self.conn.commit()
            
    def commit(self):
        self.conn.commit()

    def query_all(self):
        self.cursor.execute(f'SELECT fp FROM {self.table} ')
        rows = self.cursor.fetchall()
        vectors = [row[0] for row in rows]
        return vectors

    def close_connection(self):
        self.cursor.close()
        self.conn.close()
        
    def clear_table(self):
        self.cursor.execute(f'TRUNCATE TABLE {self.table} ')
        self.conn.commit()
    def drop_table(self):
        self.cursor.execute(f'DROP TABLE IF EXISTS  {self.table} CASCADE')
        self.conn.commit()

# Example Usage:
if __name__ == "__main__":
    # Replace with your PostgreSQL credentials and database information
    db = VectorDatabase('andr',dbname='dalvarez', user='dalvarez', password='psqldalvarez')

    # Example vectors
    # vector1 = np.array([1.0, 2.0, 3.0]).tolist()
    # vector2 = np.array([4.0, 5.0, 6.0]).tolist()

    # Insert vectors into the database
    # db.insert(['mol123123/¡', 'mol1', True, 2, 'mapA', vector1, vector1*2, vector1*3])
    # db.insert(['molsdfc223123/¡', 'mol2', True, 1, 'mapB', vector2, vector2*2, vector2*3])

    # Query vectors from the database
    # retrieved_vectors = db.query_all()
    # print(retrieved_vectors[:10])
    print(f"tmp_dx_1/100_C00109011-dec.4_hbond_Acceptors.dx EXISTS: {db.exists('tmp_dx_1/100_C00109011-dec.4_hbond_Acceptors.dx')}")
    print(f"tmp_dx_1/100_C00109011-dec.4_hbond_Donors.dx EXISTS: {db.exists('tmp_dx_1/100_C00109011-dec.4_hbond_Donors.dx')}")

    # print("Inserted Vectors:", [vector.tolist() for vector in [vector1, vector2]])
    # print("Retrieved Vectors:", [vector.tolist() for vector in retrieved_vectors])

    # Close the database connection
    db.close_connection()
