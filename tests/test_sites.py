import copy,json,unittest
from pathlib import Path
from site_orbits import process

class SiteTests(unittest.TestCase):
    def setUp(self):self.d=json.loads(Path('fixtures/input.json').read_text())
    def test_special_position_not_double_counted(self):
        r=process(self.d);self.assertEqual(r['cell_composition'],{'X':'1','Y':'1'})
        self.assertEqual(r['site_orbits'],[{'site':'a','multiplicity':2,'stabilizer_size':1},{'site':'b','multiplicity':1,'stabilizer_size':2}])
        self.assertEqual(r['positions'][0]['contributions'][0]['operations'],['identity','invert'])
    def test_translation_group_and_wrapping(self):
        self.d['operations'][1].update(matrix=[[1,0,0],[0,1,0],[0,0,1]],translation=[[1,2],[0,1],[0,1]])
        self.d['sites']=self.d['sites'][:1];self.d['sites'][0]['point'][0]=[-1,4]
        r=process(self.d);self.assertEqual([p['point'][0] for p in r['positions']],['1/4','3/4'])
    def test_disorder_and_overoccupancy(self):
        s=copy.deepcopy(self.d['sites'][0]);s.update(id='alternate',species='Z');self.d['sites'].append(s)
        r=process(self.d);self.assertEqual(r['positions'][1]['occupancy'],'1');self.assertEqual(r['cell_composition']['Z'],'1')
        s['occupancy']=[3,4]
        with self.assertRaisesRegex(ValueError,'overoccupied'):process(self.d)
    def test_group_closure(self):
        self.d['operations'][1].update(matrix=[[1,0,0],[0,1,0],[0,0,1]],translation=[[1,3],[0,1],[0,1]])
        with self.assertRaisesRegex(ValueError,'closed'):process(self.d)
    def test_missing_identity_and_duplicate_operation(self):
        for mutate in (lambda d:d['operations'].pop(0),lambda d:d['operations'][1].update(matrix=[[1,0,0],[0,1,0],[0,0,1]])):
            d=copy.deepcopy(self.d);mutate(d)
            with self.assertRaises(ValueError):process(d)
    def test_unimodular_and_fraction_validation(self):
        for mutate in (lambda d:d['operations'][1]['matrix'][0].__setitem__(0,0),lambda d:d['sites'][0].update(occupancy=[True,2])):
            d=copy.deepcopy(self.d);mutate(d)
            with self.assertRaises(ValueError):process(d)
