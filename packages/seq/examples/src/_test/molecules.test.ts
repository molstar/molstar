/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The molecules example's data says what its text says: the disulfide joins two cysteines, the glycan sits on the N of
// an N-x-S/T sequon, the sodium on an aspartate's OD1 and ATP, the modified residues on their parents' letters.

import type { Node } from '@molstar/core/repr/tree';
import { describe, expect, it } from 'vitest';
import { complexView } from '../molecules/molecules.js';

const find = (n: Node, id: string): Node | undefined =>
  n.id === id ? n : n.children?.reduce<Node | undefined>((a, c) => a ?? find(c, id), undefined);

/** A rows column as written: numbers, or a dictionary-coded string column. */
type Column = ArrayLike<number> | { dict: string[]; codes: ArrayLike<number> };
const values = (c: Column): unknown[] =>
  'dict' in c ? Array.from(c.codes, (i) => c.dict[i]) : Array.from(c as ArrayLike<number>);

describe('the molecules example', () => {
  it('its links join the residues its text names', () => {
    const root = complexView();
    const letters = find(root, 'fc')!.params!.letters as string;
    const raw = find(root, 'links')!.params!.rows as Record<string, Column>;
    const rows = Object.fromEntries(Object.entries(raw).map(([k, v]) => [k, values(v)]));
    const at = (unit: unknown, start: unknown) => (unit === 'A' || unit === 'B' ? letters[start as number] : unit);
    const links = rows.type.map((t, i) => [
      t,
      at(rows.a_unit[i], rows.a_start[i]),
      rows.a_atom[i],
      at(rows.b_unit[i], rows.b_start[i]),
      rows.b_atom[i],
    ]);
    expect(links).toEqual([
      ['covalent', 'N', 'ND2', 'glycan', 'C1'],
      ['disulfide', 'C', 'SG', 'C', 'SG'],
      ['metal', 'na', 'NA', 'D', 'OD1'],
      ['metal', 'na', 'NA', 'atp', 'O1G'],
    ]);
    // the sequon: N-x-S/T, x not P
    const n = rows.a_start[0] as number;
    expect(letters.slice(n, n + 3)).toMatch(/^N[^P][ST]$/);
    // MSE on M, SEP on S, CRO over T Y G
    expect([letters[0], letters[20], letters.slice(30, 33)]).toEqual(['M', 'S', 'TYG']);
  });
});
