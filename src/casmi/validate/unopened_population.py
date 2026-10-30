"""Select unopened family components, including transitive metric protection."""
from __future__ import annotations

import hashlib


def select_unopened(records, protected_families, protected_metrics, counts, salt):
    if set(counts) != {'dev', 'test'} or any(type(n) is not int or n < 1 for n in counts.values()):
        raise ValueError('Positive development and test counts required')
    parents, first = {}, {}
    def find(family):
        parents.setdefault(family, family)
        while parents[family] != family:
            parents[family] = parents[parents[family]]
            family = parents[family]
        return family
    for row in records:
        if not all(isinstance(row.get(k), str) and row[k] for k in ('key', 'family', 'metric')):
            raise ValueError('Valid key/family/metric records required')
        a, b = find(row['family']), find(first.setdefault(row['metric'], row['family']))
        if a != b: parents[max(a, b)] = min(a, b)
    components = {}
    for row in records: components.setdefault(find(row['family']), []).append(row)
    eligible = {family: rows for family, rows in components.items()
        if not any(r['family'] in protected_families or r['metric'] in protected_metrics for r in rows)}
    ordered = sorted(eligible, key=lambda family: hashlib.sha256((salt+'|'+family).encode()).hexdigest())
    if len(ordered) < sum(counts.values()): raise ValueError('Insufficient unopened families')
    result, offset = {}, 0
    for name in ('dev', 'test'):
        result[name] = [{**min(eligible[f], key=lambda r: r['key']), 'group': f} for f in ordered[offset:offset+counts[name]]]
        offset += counts[name]
    return result
