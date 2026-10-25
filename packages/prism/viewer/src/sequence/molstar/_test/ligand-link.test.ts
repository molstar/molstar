/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

// A ligand in the link (seqspec:libraries#the-mol-adapter): its Seq* positions are its heavy atoms, but what a person
// points at is the ligand: a position of a non-polymer unit lights, and selects, its whole residue in 3D (hydrogens
// included), not the one atom its position names.

import { readFileSync } from 'node:fs';
import { createSeqViewSpec } from '@molstar/seq-core/builder';
import { compile } from '@molstar/seq-core/format';
import { Extents } from '@molstar/seq-core/model';
import { Translator } from '@molstar/seq-core/translate';
import { StructureElement } from '@molstar/model/model/structure';
import { describe, expect, it } from 'vitest';
import { complexInputs, readStructure, StructureLink, structureOf } from '../index.js';
import { registries } from './fixtures.js';

/** ins1.cif (a nine-residue chain) with a ligand B: three heavy atoms C1, C2, O1 and a hydrogen. */
function withLigand(): Uint8Array {
  let text = readFileSync(new URL('./fixtures/structures/ins1.cif', import.meta.url), 'utf8').replace(/\r\n/g, '\n');
  text = text.replace(
    "1 polymer 'Heavy chain fragment'\n",
    "1 polymer 'Heavy chain fragment'\n2 non-polymer 'A small ligand'\n",
  );
  text = text.replace(
    '_struct_asym.id A\n_struct_asym.entity_id 1\n',
    'loop_\n_struct_asym.id\n_struct_asym.entity_id\nA 1\nB 2\n',
  );
  const atoms = [
    'HETATM 9 C C1 . LIG B 2 . ? 40.000 0.000 0.000 1.00 10.00 101 LIG L C1 1',
    'HETATM 10 C C2 . LIG B 2 . ? 41.500 0.000 0.000 1.00 10.00 101 LIG L C2 1',
    'HETATM 11 O O1 . LIG B 2 . ? 43.000 0.000 0.000 1.00 10.00 101 LIG L O1 1',
    'HETATM 12 H H1 . LIG B 2 . ? 40.000 1.000 0.000 1.00 10.00 101 LIG L H1 1',
  ];
  text = `${text.replace(/#\s*$/, '')}${atoms.join('\n')}\n#\n`;
  return new TextEncoder().encode(text);
}

describe('a ligand unit in the link', () => {
  it("a position names the whole ligand in 3D, not its one atom; a polymer's residue stays itself", async () => {
    const bytes = withLigand();
    const b = createSeqViewSpec();
    b.data.complex({
      ref: 'ab',
      provider: 'molstar:structure',
      source: { uri: 'https://example.org/lig.cif', format: 'mmcif' },
    } as never);
    b.view({ ref: 'v', label: 'V' }).track({ label: '' }).ruler();
    const c = await compile(b.getDocument(), registries(), { load: async () => bytes } as never);
    expect(c.issues.filter((i) => i.severity === 'error')).toEqual([]);
    const structure = await structureOf(await readStructure(bytes, 'mmcif'), 0, null);
    const link = StructureLink.attach(new Translator(), '3d', structure, c.set, 'ab');
    const ligand = link.chains.find((m) => m.chain.kind === 'ligand');
    const polymer = link.chains.find((m) => m.chain.kind === 'polymer');
    expect(ligand?.chain.count).toBe(3);
    const at = (m: NonNullable<typeof ligand>, start: number, stop: number) => {
      const runs = Extents.create(1);
      Extents.add(runs, m.system, start, stop);
      return link.toLoci(runs, { wholeResidue: true });
    };
    // One heavy atom's position, with the linked views' `wholeResidue` reaction: the ligand's four atoms (hydrogen
    // included); the link's own mapping stays exact (whole-residue.test.ts).
    expect(StructureElement.Loci.size(at(ligand as NonNullable<typeof ligand>, 0, 1))).toBe(4);
    expect(StructureElement.Loci.size(at(ligand as NonNullable<typeof ligand>, 2, 3))).toBe(4);
    // A polymer position is its residue, as it was.
    expect(StructureElement.Loci.size(at(polymer as NonNullable<typeof polymer>, 2, 3))).toBe(1);
  });
});

describe('the ligands a complex includes', () => {
  it('`ligand_comps` keeps the ligand units of those chemical components, and a polymer stays', async () => {
    const src = await readStructure(withLigand(), 'mmcif');
    const structure = await structureOf(src, 0, null);
    const units = (ligand_comps?: string[]) =>
      complexInputs([{ structure, model: null }], src.frame, ligand_comps ? { ligand_comps } : {}, () => {}).units.map(
        (u) => u.key,
      );
    expect(units()).toEqual(['A', 'B']);
    expect(units(['LIG'])).toEqual(['A', 'B']);
    expect(units(['HEM'])).toEqual(['A']);
    expect(units([])).toEqual(['A']);
  });
});
