from .checks import fields,array,integer,ident,unique,fraction
from .affine import key,determinant,compose,IDENTITY

def parse_operations(data):
    operations=[]
    for o in array(data['operations'],1,24):
        fields(o,('id','matrix','translation'));name=ident(o['id'])
        m=tuple(tuple(integer(v,-2,2) for v in array(row,3,3)) for row in array(o['matrix'],3,3))
        if abs(determinant(m))!=1:raise ValueError('unimodular operation required')
        t=tuple(fraction(v) for v in array(o['translation'],3,3));operations.append((name,key(m,t)))
    unique([o[0] for o in operations]);unique([o[1] for o in operations]);keys={o[1] for o in operations}
    if IDENTITY not in keys:raise ValueError('identity operation missing')
    for _,a in operations:
        for _,b in operations:
            if compose(a,b) not in keys:raise ValueError('operation set is not closed')
    return operations

def parse(data):
    fields(data,('operations','sites'));operations=parse_operations(data);sites=[]
    for s in array(data['sites'],1,200):
        fields(s,('id','species','point','occupancy'));name=ident(s['id']);species=ident(s['species'])
        p=tuple(fraction(v)%1 for v in array(s['point'],3,3));occ=fraction(s['occupancy'],True)
        if occ>1:raise ValueError('occupancy exceeds one')
        sites.append({'id':name,'species':species,'point':p,'occupancy':occ})
    unique([s['id'] for s in sites]);return operations,sites
