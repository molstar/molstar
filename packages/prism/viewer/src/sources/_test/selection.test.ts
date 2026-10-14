/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { readFileSync } from 'node:fs';
import { OrderedSet } from '@molstar/core/data/int';
import { StructureElement } from '@molstar/model/model/structure';
import { type ChainInstance, chainsOf, readStructure, structureOf } from '@molstar/seq-molstar-poc';
import { describe, expect, it } from 'vitest';
import { selectionSummary } from '../selection.js';

// The selection bar's summary over a recorded structure (1TUP: p53 and DNA chains, zinc ions): one thing is named, more
// are counted by kind, nothing is nothing.

const FILE = '../../../../../sequences/molstar-poc/_test/fixtures/structures/1tup.bcif';

describe('the selection summary', () => {
  const setup = async () => {
    const bytes = new Uint8Array(readFileSync(new URL(FILE, import.meta.url)));
    const structure = await structureOf(await readStructure(bytes, 'bcif'), 0, null);
    const { chains } = chainsOf(structure);
    const loci = (c: ChainInstance, from: number, to: number) => {
      const at: number[] = [];
      for (let e = from; e < to; e++) for (let i = c.start[e]; i < c.end[e]; i++) at.push(i);
      const indices = OrderedSet.ofSortedArray(Int32Array.from(at)) as OrderedSet<StructureElement.UnitIndex>;
      return StructureElement.Loci(structure, [{ unit: c.unit, indices }]);
    };
    return {
      loci,
      protein: chains.find((c) => c.kind === 'polymer' && c.entityId === '3') as ChainInstance,
      ligand: chains.find((c) => c.kind === 'ligand') as ChainInstance,
    };
  };

  it('names a residue, a chain and a ligand; counts a range', async () => {
    const { loci, protein, ligand } = await setup();
    expect(selectionSummary([loci(protein, 3, 4)])).toMatch(/^[A-Z]{3} \d+$/);
    expect(selectionSummary([loci(protein, 3, 15)])).toBe('12 residues');
    expect(selectionSummary([loci(protein, 0, protein.count)])).toBe(protein.labelAsymId);
    expect(selectionSummary([loci(ligand, 0, ligand.count)])).toBe('ZN 951');
  });

  it('counts a mix, and says nothing for nothing', async () => {
    const { loci, protein, ligand } = await setup();
    expect(selectionSummary([loci(protein, 3, 15), loci(ligand, 0, ligand.count)])).toBe('1 chain + 12 residues');
    expect(selectionSummary([])).toBeUndefined();
    expect(selectionSummary([StructureElement.Loci(loci(protein, 0, 1).structure, [])])).toBeUndefined();
  });
});
