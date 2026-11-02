import unittest
from casmi.validate.unopened_population import select_unopened


class UnopenedTests(unittest.TestCase):
    def test_transitive_protection_precedes_representative_selection(self):
        rows = [dict(key=str(i), family=f, metric=m) for i,(f,m) in enumerate([
            ('opened','bridge'),('linked','bridge'),('linked','second'),('third','second'),
            ('new1','m1'),('new2','m2'),('new3','m3')])]
        chosen = select_unopened(rows, {'opened'}, set(), {'dev':1,'test':2}, 'fixed')
        self.assertEqual({r['family'] for v in chosen.values() for r in v}, {'new1','new2','new3'})
        self.assertEqual(chosen, select_unopened(list(reversed(rows)), {'opened'}, set(), {'dev':1,'test':2}, 'fixed'))

    def test_protected_metric_and_insufficient_population_fail(self):
        rows = [dict(key='a',family='f',metric='used'), dict(key='b',family='g',metric='fresh')]
        with self.assertRaises(ValueError): select_unopened(rows, set(), {'used'}, {'dev':1,'test':1}, 's')
