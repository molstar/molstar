/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { readFileSync } from 'node:fs';
import { gunzipSync } from 'node:zlib';
import type { Node } from '@molstar/core/repr/tree';
import { describe, expect, it } from 'vitest';
import { Data, entryLigands } from '../complex/data.js';
import { build, transplantUnit } from '../complex/units.js';
import { PackageFolder } from '../node-assets.js';

// 1FIN as a complex: the units' author numbers the entry's (the whole mmCIF in the medium lane;
// the spec's units, ATP's atoms and the transplants' scores quick), the contacts' rows per pair of units.

const find = (n: Node, id: string): Node | undefined =>
  n.id === id ? n : n.children?.reduce<Node | undefined>((a, c) => a ?? find(c, id), undefined);
type Category = { dict: string[]; codes: ArrayLike<number> };
const strings = (c: unknown) => {
  const x = c as Category | string[];
  return Array.isArray(x) ? x : Array.from(x.codes, (k) => x.dict[k]);
};
const labelsOf = (spec: Node, unit: string) => strings(find(find(spec, '1fin')!, unit)!.children![0].params!.labels);

describe('1FIN as a complex', () => {
  const spec = build();

  it('its units: CDK2 A and C, cyclin A B and D, the ATPs, the transplants; ATP with its 31 heavy atoms', () => {
    const units = find(spec, '1fin')!.children!.filter((c) => c.kind === 'unit');
    expect(units.map((u) => u.id)).toEqual([
      'A',
      'B',
      'C',
      'D',
      'atp-a',
      'atp-c',
      ...entryLigands().map((l) => transplantUnit(l.entry)),
    ]);
    expect(strings((find(spec, 'atp')!.params!.atoms as { name: unknown }).name)).toHaveLength(31);
    expect(labelsOf(spec, 'B')[0]).toBe(String(Data.cyclinFirst));
    // each transplant its resolution as a score, a field of its row
    const score = find(spec, `${transplantUnit(entryLigands()[0].entry)}-score`)!.params!.rows as Record<
      string,
      unknown
    >;
    expect(Array.from(score.resolution as ArrayLike<number>)).toEqual([entryLigands()[0].resolution]);
  });

  // the whole entry parsed: slow by nature, the medium lane
  it('every modeled residue’s author number is its unit’s label', { tags: ['medium'] }, async () => {
    // the parser imported here: only this medium test needs it (the quick lane does not load @molstar/io)
    const { readCif } = await import('../../scripts/cif.js');
    const cif = await readCif(gunzipSync(readFileSync(`${PackageFolder}fixtures/structures/1FIN.cif.gz`)));
    for (const ch of Data.chains) {
      const labels = labelsOf(spec, ch.id);
      let n = 0;
      for (const r of cif.rows('pdbx_poly_seq_scheme'))
        if (r.asym_id === ch.id && r.auth_seq_num !== '?') {
          expect(labels[Number(r.seq_id) - 1], `${ch.id} ${r.seq_id}`).toBe(r.auth_seq_num);
          n++;
        }
      expect(n).toBeGreaterThan(200);
    }
  });
});
