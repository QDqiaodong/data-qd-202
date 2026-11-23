from fractions import Fraction
from .affine import apply

def expand(operations,sites):
    positions={};orbits=[]
    for site in sites:
        orbit={}
        for name,op in operations:orbit.setdefault(apply(op,site['point']),[]).append(name)
        orbits.append({'site':site['id'],'multiplicity':len(orbit),'stabilizer_size':len(operations)//len(orbit)})
        if len(operations)%len(orbit):raise ValueError('orbit size invariant')
        for p,op_ids in orbit.items():
            positions.setdefault(p,[]).append({'site':site['id'],'species':site['species'],'occupancy':site['occupancy'],'operations':sorted(op_ids)})
    records=[];composition={}
    for p,contributions in sorted(positions.items()):
        total=sum((c['occupancy'] for c in contributions),Fraction(0))
        if total>1:raise ValueError('overoccupied expanded position')
        for c in contributions:composition[c['species']]=composition.get(c['species'],Fraction(0))+c['occupancy']
        records.append({'point':[str(v) for v in p],'occupancy':str(total),'vacancy':str(1-total),'contributions':[{**c,'occupancy':str(c['occupancy'])} for c in sorted(contributions,key=lambda c:c['site'])]})
    return records,orbits,{k:str(v) for k,v in sorted(composition.items())}
