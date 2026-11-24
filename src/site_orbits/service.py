from .model import parse
from .expand import expand
from .recover import recover_asymmetric_unit

def process(data):
    operations,sites=parse(data);positions,orbits,composition=expand(operations,sites)
    return {'positions':positions,'site_orbits':orbits,'cell_composition':composition,'summary':{'operations':len(operations),'source_sites':len(sites),'expanded_positions':len(positions),'contributions':sum(len(p['contributions']) for p in positions)}}
