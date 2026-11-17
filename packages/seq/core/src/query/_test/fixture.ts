/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A small model the query tests read: a spec tree and the systems, numberings, annotations and dynamic nodes it
// names. The static model implements the same interface later; this one is written out by hand.

import type { Atom, Node } from '@molstar/core/repr/tree';
import { type Location, location, partOf } from '../location.js';
import type { ItemsRead, Numbering, QueryAnnotation, QueryItem, QueryModel, QuerySystem } from '../model.js';

/** A numbering from a label per position (undefined: a hole); a repeated label names its first position. */
export function labelsNumbering(id: string, labels: readonly (string | undefined)[]): Numbering {
  const first = new Map<string, number>();
  labels.forEach((l, i) => {
    if (l !== undefined && !first.has(l)) first.set(l, i);
  });
  return { id, length: labels.length, labelAt: (i) => labels[i], indexOf: (l) => first.get(l) ?? -1 };
}

/** A numbering by offset: position i is `i + by`, skipping zero when `skipZero` (…, -1, 1, …). */
export function offsetNumbering(id: string, length: number, by: number, skipZero = false): Numbering {
  const labels: string[] = [];
  for (let i = 0; i < length; i++) {
    const n = i + by;
    labels.push(`${skipZero && n >= 0 ? n + 1 : n}`);
  }
  return labelsNumbering(id, labels);
}

function system(
  address: string,
  length: number,
  numberings: readonly Numbering[] = [],
  more: Partial<Pick<QuerySystem, 'cyclic' | 'parent'>> = {},
): QuerySystem {
  const by = new Map(numberings.map((n) => [n.id, n]));
  return { address, length, numbering: (id) => by.get(id), numberingIds: () => [...by.keys()], ...more };
}

/** An annotation whose items are loaded only over `loaded` (all of them when undefined). */
function annotation(address: string, target: string, rows: readonly QueryItem[], loaded?: [number, number]) {
  const fields = [...new Set(rows.flatMap((r) => Object.keys(r.fields ?? {})))];
  const a: QueryAnnotation & { reads: number } = {
    address,
    target,
    fields,
    reads: 0,
    items(over?: Location): ItemsRead {
      a.reads++;
      const ranges = over && partOf(over, target)?.ranges;
      const asked: [number, number][] = [];
      if (ranges) for (let i = 0; i < ranges.length; i += 2) asked.push([ranges[i], ranges[i + 1]]);
      else asked.push([0, 2 ** 31 - 1]);
      if (!loaded) return { items: rows };
      const [lo, hi] = loaded;
      const pending: number[] = [];
      for (const [s, e] of asked) {
        if (s < lo) pending.push(s, Math.min(e, lo));
        if (e > hi) pending.push(Math.max(s, hi), e);
      }
      const items = rows.filter((r) => r.start < hi && r.end > lo);
      return { items, pending: location({ parts: [{ system: target, ranges: pending }] }) };
    },
  };
  return a;
}

// Kabat-like: 82, 82A, 82B, 82C, 83 with a hole; IMGT-like: CDR3 labels not in index order (111, 111A, 112A, 112).
const AbKabat = labelsNumbering('kabat', [
  '80',
  '81',
  '82',
  '82A',
  '82B',
  '82C',
  '83',
  undefined,
  '85',
  '86',
  '87',
  '88',
]);
const AbImgt = labelsNumbering('imgt', [
  '108',
  '109',
  '110',
  '111',
  '111A',
  '112A',
  '112',
  '113',
  '114',
  '115',
  '116',
  '117',
]);
const AbAuth = offsetNumbering('auth', 12, 200);
/** The unit's own `auth` beats the molecule's; it has none of `kabat`, so the molecule's is used, shifted by 2. */
const UnitAuth = offsetNumbering('auth', 8, 1);

const Ab = system('/data/ab', 12, [AbKabat, AbImgt, AbAuth]);

export const Systems = {
  p53: system('/data/p53', 393, [offsetNumbering('uniprot', 393, 1), labelsNumbering('dup', ['1', '2', '2', '3'])]),
  // `gapped`: 1-based, with no label at positions 5 and 96
  plasmid: system(
    '/data/plasmid',
    100,
    [
      labelsNumbering(
        'gapped',
        Array.from({ length: 100 }, (_, i) => (i === 4 || i === 95 ? undefined : `${i + 1}`)),
      ),
    ],
    { cyclic: true },
  ),
  // a chain is a system data makes, addressed by its form
  chainA: system('/data/1tup/(chain :id A)', 6, [labelsNumbering('auth', ['92', '93', '94', '94A', '95', '96'])]),
  ab: Ab,
  unit: system('/data/abunit', 8, [UnitAuth], { parent: { system: Ab, start: 2 } }),
  empty: system('/data/empty', 0),
  offset: system('/data/offset', 6, [
    offsetNumbering('skip', 6, -3, true),
    // `.n` labels, `100.10` distinct from `100.1`
    labelsNumbering('dot', ['100', '100.1', '100.10', '101', '102', '103']),
  ]),
};

export const Annotations = {
  domains: annotation('/data/p53/domains', '/data/p53', [
    { key: 'TAD', start: 0, end: 61, fields: { label: 'transactivation', type: 'domain' } },
    { key: 'DBD', start: 93, end: 292, fields: { label: 'DNA-binding', type: 'domain' } },
    { key: 'TET', start: 324, end: 355, fields: { label: 'tetramerization', type: 'domain' } },
    { key: 1, start: 300, end: 300, fields: { label: 'site', type: 'site' } },
    { key: '1', start: 100, end: 101, fields: { label: 'one', type: 'site' } },
  ]),
  // a wrapped feature is two rows of one item
  features: annotation('/data/plasmid/features', '/data/plasmid', [
    { key: 'ori', start: 90, end: 100 },
    { key: 'ori', start: 0, end: 5 },
    { key: 'gene', start: 20, end: 40 },
    // a boundary at the system's end
    { key: 'end', start: 100, end: 100 },
  ]),
  sites: annotation('/data/1tup/sites', '/data/1tup/(chain :id A)', [{ key: 's1', start: 3, end: 4 }]),
  // read by region: loaded over [0, 50) only
  peaks: annotation(
    '/data/p53/peaks',
    '/data/p53',
    [
      { key: 'a', start: 10, end: 20 },
      { key: 'b', start: 40, end: 45 },
    ],
    [0, 50],
  ),
};

const node = (kind: string, id: string, children?: Node[]): Node => (children ? { kind, id, children } : { kind, id });

export const Spec: Node = {
  kind: 'spec',
  children: [
    node('data', 'data', [
      node('sequence', 'p53', [node('ranges', 'domains'), node('ranges', 'peaks')]),
      node('sequence', 'plasmid', [node('ranges', 'features')]),
      node('sequence', 'ab'),
      node('sequence', 'abunit'),
      node('sequence', 'empty'),
      node('sequence', 'offset'),
      node('genes', 'genes'),
      node('structure', '1tup', [node('ranges', 'sites')]),
    ]),
    {
      kind: 'refs',
      id: 'refs',
      children: [{ kind: 'ref', id: 'dbd', params: { to: '/data/p53/(range :from 94 :to 292)' } }],
    },
    node('views', 'views', [node('view', 'main', [node('section', 'protein'), node('section', 'gene')])]),
  ],
};

const SystemsByAddress = new Map(Object.values(Systems).map((s) => [s.address, s]));
const AnnotationsByAddress = new Map<string, QueryAnnotation>(Object.values(Annotations).map((a) => [a.address, a]));

export const Model: QueryModel = {
  system: (a) => SystemsByAddress.get(a),
  annotation: (a) => AnnotationsByAddress.get(a),
  node: (owner: string, key: Atom) =>
    owner === '/data/genes' && (key === 'T1' || key === 'T2') ? { kind: 'gene', params: { key } } : undefined,
  chain: (owner: string, id: Atom) =>
    owner === '/data/1tup' && id === 'A' ? { kind: 'chain', params: { id } } : undefined,
};
