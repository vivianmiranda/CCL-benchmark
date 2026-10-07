"""Regression checks showing that quadrature refinement detects the narrow-input error."""
import argparse
import json
from pathlib import Path
import unittest
import numpy as np
from common import digest, write_json


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--case',required=True)
    a=p.parse_args()
    source=a.output/'cocoa_integration_ladder'/f'{a.case}.json'
    data=json.loads(source.read_text())
    assert data['status']=='completed' and data['integration_ladder']
    values=np.array([r['values'] for r in data['pairs']]).T
    assert np.isfinite(values).all() and (values>0).all()
    assert [r['integration_accuracy'] for r in data['controls']]==list(range(5))
    assert data['nodes']==[96,128,256,512,1024]
    class NarrowOverlapConvergence(unittest.TestCase):
        def test_default_fails_one_part_in_ten_thousand(self):
            with self.assertRaises(AssertionError):
                np.testing.assert_allclose(values[0],values[4],rtol=1e-4,atol=0)
        def test_level_two_resolves_all_fifteen_pairs(self):
            np.testing.assert_allclose(values[2],values[4],rtol=1e-4,atol=0)
        def test_reference_agrees_with_previous_level(self):
            np.testing.assert_allclose(values[3],values[4],rtol=1e-5,atol=0)
    result=unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(NarrowOverlapConvergence))
    record={'source_sha256':digest(source),'script_sha256':digest(__file__),
        'tests_run':result.testsRun,'successful':result.wasSuccessful(),
        'convergence_rtol':1e-4,'reference_rtol':1e-5,'atol':0,
        'max_fractional_change_all_pairs_by_level':np.max(abs(values/values[-1]-1),axis=1).tolist(),
        'interpretation':'The expected default rejection is a regression diagnostic, not a passing convergence result. Level 2 passes 0.01% against level 4 on all 15 pairs; levels 3/4 also agree within 0.001%. This tests quadrature only.'}
    write_json(source.parent/'regression_check.json',record)
    if not result.wasSuccessful():raise SystemExit(1)


if __name__=='__main__':main()
