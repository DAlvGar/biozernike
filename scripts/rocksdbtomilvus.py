import rocksdb
import glob
import re
import numpy as np
import time
from pymilvus import (
    connections,
    utility,
    FieldSchema, CollectionSchema, DataType,
    Collection,
)

target="braf"
N_dbs=4

index = {
    "index_type": "IVF_FLAT",
    "metric_type": "L2",
    "params": {"nlist": 128},
    }

print("Connecting to Milvus")
connections.connect("default", host="localhost", port="19530")
has = utility.has_collection(target)
if has:
   utility.drop_collection(target)
   has=False

# class ConformerMapMoments(BaseDoc):
#   molid: str = ''
#   active: bool
#   fp: NdArray[72]
#   moments: NdArray[444]
#   map: str = ''
#   conformer: int
#   original: str

if not has:
    # Create schema
    print("Create Milvus collection "+target)
    # fields = [
    #     FieldSchema(name="original", dtype=DataType.VARCHAR, is_primary=True, auto_id=False, max_length=150),
    #     FieldSchema(name="molid", dtype=DataType.VARCHAR, max_length=150),
    #     FieldSchema(name="active", dtype=DataType.BOOL),
    #     FieldSchema(name="conformer", dtype=DataType.INT16),
    #     FieldSchema(name="map", dtype=DataType.VARCHAR, max_length=150),
    #     FieldSchema(name="fp", dtype=DataType.FLOAT_VECTOR, dim=72),
    #     FieldSchema(name="moments_real", dtype=DataType.ARRAY, element_type=DataType.FLOAT, max_capacity=444),
    #     FieldSchema(name="moments_imag", dtype=DataType.ARRAY, element_type=DataType.FLOAT, max_capacity=444)
    # ]
    fields = [
        FieldSchema(name="original", dtype=DataType.VARCHAR, is_primary=True, auto_id=False, max_length=100),
        FieldSchema(name="fp", dtype=DataType.FLOAT_VECTOR, dim=72),
        ]
    schema = CollectionSchema(fields, "Zernike moments DB")
    target_col = Collection(target, schema, consistency_level="Strong")
else:
    target_col = Collection(target)
    target_col.drop_index()
    
# INSERT DATA

def parse_string(input_string):
    pattern = r".*/\d+_(?P<molid>[A-Za-z0-9]+)-(?P<active>[A-Za-z]+_*[0-9]*)\.(?P<conformer>\d+)_(?P<map>[A-Za-z_]+)\.dx"
    match = re.match(pattern, input_string)
    if match:
        result_dict = match.groupdict()
        result_dict["conformer"] = int(result_dict["conformer"])
        result_dict["original"] = input_string
        if '_' in result_dict['active']:
            result_dict['active'] = result_dict['active'].split('_')[0]
        if 'dec' in result_dict['active']: result_dict['active'] = False
        else: result_dict['active'] = True
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
    return (real, imag)
    
def parseEntry(e):
    info = parseKey(e[0])
    info['fp'] = parseFp(e[1])
    if not info['fp']: return None
    return info

failed=[]
start=time.time()
for i in range(1, N_dbs+1):
    dbstart=time.time()
    f = f'fp_{i}.db'
    # m = f'moments_{i}.db'
    print('parsing '+f)
    db=rocksdb.DB(f, rocksdb.Options(create_if_missing=False))
    # mdb=rocksdb.DB(m, rocksdb.Options(create_if_missing=False))
    it=db.iteritems()
    it.seek_to_first()
    # data = dict(zip(keys, [[]]*len(keys)))
    data = []
    i=0
    while (x := next(it, None)) is not None:
        i += 1
        k=x[0].decode('utf-8')
        v = parseFp(x[1])
        if not v:
            failed.append(k)
            continue    
        data.append({'original':k,'fp':v})
        if i == 200:
            print('Inserting 100...')
            target_col.insert(data)
            data = []
            i = 0
        #e['moments_real'], e['moments_imag'] = parseMoments(mdb.get(x[0]))
        # data.append(e)
        #print((k,v))
        # target_col.insert()
        #print(e)
        #for k in keys:
        #    data[k].append(e[k])
    print(f'File time: {time.time() - dbstart} Total time: {time.time()-start}')
    # add data to milvus
    # entries = {"rows":data}
    # target_col.insert(entries)

target_col.flush()
    
# ADD BACK INDEX
target_col.create_index("fp", index)
connections.disconnect("default")
print('done')
np.savetxt('milvus_fail_import.log',failed, fmt='%s')

# import rocksdb
# k=[]
# for i in range(1,5):
#      d = rocksdb.DB(f'fp_{i}.db', rocksdb.Options(create_if_missing=False))
#      it = d.iterkeys()
#      it.seek_to_first()
#      k+=list(it)