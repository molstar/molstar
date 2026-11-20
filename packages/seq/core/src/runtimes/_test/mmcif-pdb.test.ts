/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The mmCIF runtime's conformance on the larger PDB entries (molecules.md, Reading): 1FIN, 4HHB and 3SGJ read against
// the counts recorded from each file, counted again from the file's categories, and valid by the language. Medium
// lane: each file is a few hundred kilobytes parsed; the quick representative is 1EHZ (`mmcif.test.ts`).

import { describe, expect, it } from 'vitest';
import { countedFromFile, type EntryId, expectValid, frameOf, openEntry, Recorded, summary } from './mmcif-entries.js';

describe('the mmCIF runtime on PDB entries', { tags: ['medium'] }, () => {
  for (const id of ['1FIN', '4HHB', '3SGJ'] as const satisfies readonly EntryId[]) {
    it(`${id}: every entity, chain, link and author number as recorded, counted again from the file (EC-118, EC-65)`, async () => {
      const source = await openEntry(id);
      expect(summary(source)).toEqual(Recorded[id]);
      const counted = countedFromFile(await frameOf(id));
      const s = Recorded[id];
      expect(Object.keys(s.sequences).length + Object.keys(s.molecules).length).toBe(counted.entities);
      expect(Object.keys(s.units).length).toBe(counted.chains);
      expect(s.links).toEqual(counted.links);
      for (const [e, n] of Object.entries(counted.lengths)) expect(s.sequences[e]).toMatch(new RegExp(`^\\w+ ${n} `));
      expectValid(source.coverage.nodes!);
    });
  }

  it('3SGJ: a glycan’s identifiers, its composition as a name, and the Asn–NAG links glycosidic', async () => {
    const source = await openEntry('3SGJ', { name: 'composition' });
    const glycan = source.coverage.nodes!.find((n) => n.key === '3')!.node;
    expect(glycan.params?.name).toBe('NAG2 BMA MAN NAG MAN NAG FUC');
    const ids = glycan.params?.identifiers as { value: string[] };
    expect(ids.value.some((v) => v.startsWith('WURCS=2.0/4,8,7/'))).toBe(true);
    // each N-glycan's first NAG on an Asn of the Fc (A, B) or of the receptor (C), its own links within the unit
    const links = source.coverage
      .nodes!.find((n) => n.key === 'complex')!
      .node.children!.find((c) => c.id === 'links')!;
    const rows = links.params!.rows as Record<string, any>;
    const pairs = Array.from({ length: rows.a_start.length }, (_, i) => `${rows.a_unit.get(i)}:${rows.b_unit.get(i)}`);
    expect(pairs.filter((p) => p === 'A:D' || p === 'B:E' || p === 'C:F' || p === 'C:G').length).toBe(4);
    expect(pairs.filter((p) => p === 'D:D').length).toBe(7);
  });

  it('4HHB: B-factors by unit, the default field when no pLDDT is there', async () => {
    const source = await openEntry('4HHB');
    expect(source.coverage.shapes[0].fields?.map((f) => f.name)).toEqual(['b_factor', 'occupancy']);
    const answer = await source.read(
      { system: 'A', start: 0, end: 141, level: 1, shapes: ['values'] },
      new AbortController().signal,
    );
    if (!answer.ok || answer.tables.values?.form !== 'positions') throw new Error('no values');
    expect(answer.tables.values.value.length).toBe(141);
    expect(answer.tables.values.value.every((v) => Number.isFinite(v) && v > 0)).toBe(true);
  });

  it('4HHB assembly 1: the deposited chains under its one operator, each a unit keyed by chain and operator', async () => {
    const source = await openEntry('4HHB', { assembly: '1' });
    expect(source.coverage.systems.map((s) => s.name).slice(0, 4)).toEqual([
      'A-ASM-1',
      'B-ASM-1',
      'C-ASM-1',
      'D-ASM-1',
    ]);
    expect(summary(source).links).toEqual({ metal: 4 });
  });
});
