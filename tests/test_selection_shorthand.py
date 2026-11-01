"""Consistent residue shorthand at prompt, config and planning boundaries."""
import contextlib
import io
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
from mdwb.selections import normalize_mask
from mdwb.config import normalize_config, load_config, resolve_selections
from mdwb.analysis import resolve_mask, build_pooled_batches
from mdwb.cli import prompt_mask, _quick_geometry
from mdwb.runner import make_plan
from test_config_topology import amber_metadata


class SelectionShorthandTests(unittest.TestCase):
    def test_ranges_and_lists_are_residues(self):
        for raw, expected in [('6-8',':6-8'), ('7',':7'), ('6-8,12,15-17',':6-8,12,15-17'),
                              (' 6 - 8 , 12 ',':6-8,12')]:
            self.assertEqual(normalize_mask(raw),expected)

    def test_explicit_atom_masks_and_named_selections_are_preserved(self):
        for value in ('@6-8',':6-8@CA','(:6-8)&!@/H','backbone','site6-8',":6@C4'"):
            self.assertEqual(normalize_mask(value),value)
        self.assertEqual(resolve_mask('site',{'site':'6-8'}),':6-8')

    def test_malformed_ranges_are_not_silently_reinterpreted(self):
        for value in ('8-6','0','6--8','6,','6 8','-6','6-0'):
            with self.subTest(value=value), self.assertRaises(ValueError):
                normalize_mask(value)

    def test_saved_config_normalizes_every_selection_field(self):
        data={'topology':'top.prmtop','replicas':[{'name':'A','trajectories':['a.nc']}],
              'fit_mask':'6-8','advanced':{'matrix_mask':'6-8'},'imaging':{'anchor':'6-8'},
              'sections':[{'name':'site','mask':'6-8'}],
              'interactions':[{'name':'contact','mask1':'6-8','mask2':'9-10'}],
              'monitors':[{'name':'distance','type':'distance','masks':['6-8','9-10']}],
              'analysis':{'additional':['watershell']},
              'solvation':[{'name':'shell','type':'watershell','mask1':'6-8','mask2':':WAT','lower':3,'upper':5}]}
        with tempfile.TemporaryDirectory() as temp:
            path=Path(temp)/'saved.json'; path.write_text(json.dumps(data))
            cfg=load_config(path,check_files=False)
        self.assertEqual(cfg['fit_mask'],':6-8')
        self.assertEqual(cfg['advanced']['matrix_mask'],':6-8')
        self.assertEqual(cfg['imaging']['anchor'],':6-8')
        self.assertEqual(cfg['sections'][0]['mask'],':6-8')
        self.assertEqual(cfg['interactions'][0]['mask1'],':6-8')
        self.assertEqual(cfg['monitors'][0]['masks'],[':6-8',':9-10'])
        self.assertEqual(cfg['solvation'][0]['mask1'],':6-8')
        self.assertEqual(data['fit_mask'],'6-8')

    def test_prompt_and_geometry_shorthand_share_normalization(self):
        with patch('builtins.input',return_value='6-8'), contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(prompt_mask('Global alignment mask','backbone'),':6-8')
        self.assertEqual(_quick_geometry('6-8;9:CA',2),[':6-8',':9@CA'])

    def test_advanced_plan_accepts_old_config_without_rerunning_wizard(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp); top=root/'top.prmtop'
            top.write_text(amber_metadata([('ALA',['N','CA','C','O'])]*10))
            cfg=normalize_config({'topology':str(top),'replicas':[{'name':'A','trajectories':['a.nc']}],
                'fit_mask':'6-8','sections':[{'name':'site','mask':'6-8'}],
                'advanced':{'matrix_mask':'6-8'},'analysis':{'level':'advanced'}},check_files=False)
            # Direct resolution also accepts raw caller data that bypass normalization.
            cfg['advanced']['matrix_mask']='6-8'
            selections,_=resolve_selections(cfg)
            self.assertEqual(cfg['advanced']['matrix_mask'],':6-8')
            batches=build_pooled_batches(cfg,selections)
            self.assertTrue(all(b['coordinate_mask']==':6-8' for b in batches))
            plan=make_plan(cfg,root/'plan')
            self.assertIn(':6-8',(plan/'PLAN.md').read_text())


if __name__ == '__main__':
    unittest.main()
