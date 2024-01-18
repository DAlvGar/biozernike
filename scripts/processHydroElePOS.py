from vectorDBSQL import VectorDatabase
import rocksdb
import re
import numpy as np
import time
import sys

N_dbs=int(sys.argv[2])
target=sys.argv[1]
start_at = 1
check_db = False
if len(sys.argv) > 3:
    start_at = int(sys.argv[3])
if len(sys.argv) > 4:
    check_db=bool(sys.argv[4])
psql = VectorDatabase(target,dbname='dalvarez', user='dalvarez', password='psqldalvarez')
#psql.clear_table()
#psql.drop_table()

# check if its hydroele and contains also POS or pos
def isPosMap(input_string):
   return ('hydroele') in input_string and (('_POS' in input_string) or ('_pos' in input_string))

def parse_string(input_string):
    pattern = r".*/\d+_(?P<molid>[A-Za-z0-9]+)-(?P<active>[A-Za-z]+[A-Za-z0-9_]*)\.(?P<conformer>\d+)_(?P<map>[A-Za-z_]+)\.dx"
    match = re.match(pattern, input_string)
    if match:
        result_dict = match.groupdict()
        result_dict["conformer"] = int(result_dict["conformer"])
        result_dict["original"] = input_string
        if '_' in result_dict['active']:
            result_dict['active'] = result_dict['active'].split('_')[0]
            result_dict["conformer"] += 100 # avoid repeated conformer numbers
        if 'dec' in result_dict['active']: result_dict['active'] = False
        else: result_dict['active'] = True
        
        # special treatment for hydroele_pos
        if isPosMap(input_string): result_dict['map'] = 'hydroele_pos'
        
        return result_dict
    else:
        print("ERROR WITH: "+input_string)
        return None

keys = ['original','molid', 'active', 'conformer','map','fp','moments_real','moments_imag']


def parseKey(k):
    k = k.decode('utf-8')
    return parse_string(k)

def parseFp(v):
    f = list(map(float, v.decode('utf-8').split(',')))
    if np.isnan(f[0]): return None
    return f

def parseMoments(v):
    clist = list(map(complex, v.decode('utf-8').split(',')))
    real = list(map(lambda c: c.real, clist))
    imag = list(map(lambda c: c.imag, clist))
    del clist
    return (real, imag)
    
def parseEntry(e):
    info = parseKey(e[0])
    info['fp'] = parseFp(e[1])
    if not info['fp']: return None
    return info

keys=['original', 'molid', 'active', 'conformer', 'map', 'fp', 'moments_real', 'moments_imag']
failed = open("FailedCodes.txt", "w")
start=time.time()
for i in range(start_at, N_dbs+1):
    chunk=0
    dbstart=time.time()
    f = f'pos_fp_{i}.db'
    m = f'pos_moments_{i}.db'
    print('parsing '+f)
    db=rocksdb.DB(f, rocksdb.Options(create_if_missing=False))
    mdb=rocksdb.DB(m, rocksdb.Options(create_if_missing=False))
    it=db.iteritems()
    it.seek_to_first()
    data = []
    i=0
    while (x := next(it, None)) is not None:
        k = x[0].decode('utf8')
        if not isPosMap(k): continue # focus only on the ones we need to process now :)
        if (check_db and psql.exists(k)): continue
        e = parseEntry(x)
        if not e:
            print('.', end='', flush=True)
            #print("Failed "+x[0].decode('utf8'), x[1])
            failed.write(x[0].decode('utf8')+"\n")
            continue    
        i += 1
        e['moments_real'], e['moments_imag'] = parseMoments(mdb.get(x[0]))
        data.append([e[j] for j in keys])
        if i > 1000:
            print(f'Commiting 1000 Chunk {chunk} - {time.time() - dbstart}')
            psql.insert_many(data)
            # psql.commit()
            data = []
            chunk+=1
            i = 0
    print(f'File time: {time.time() - dbstart} Total time: {time.time()-start}')
    
psql.close_connection()
failed.close()
print('done')
#np.savetxt('psql_fail_import.log',failed, fmt='%s')
