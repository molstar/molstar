/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A large set's summary in 3D (`wholeChains`, link-views.ts): the whole chains its runs touch. Mol* merges consecutive
// one-atom chains (ions) into one unit: the chains between two touched ones are not lit.

import { readFileSync } from 'node:fs';
import { createSeqViewSpec } from '@molstar/seq-core/builder';
import { compile } from '@molstar/seq-core/format';
import { Extents } from '@molstar/seq-core/model';
import { Translator } from '@molstar/seq-core/translate';
import { StructureElement } from '@molstar/model/model/structure';
import { describe, expect, it } from 'vitest';
import { readStructure, StructureLink, structureOf } from '../index.js';
import { wholeChains } from '../link-views.js';
import { registries } from './fixtures.js';

/** ins1.cif (a nine-residue chain A) with three zinc ions, chains B, C and D of one author chain (one Mol* unit). */
function withIons(): Uint8Array {
  let text = readFileSync(new URL('./fixtures/structures/ins1.cif', import.meta.url), 'utf8');
  text = text.replace(
    "1 polymer 'Heavy chain fragment'\n",
    "1 polymer 'Heavy chain fragment'\n2 non-polymer 'ZINC ION'\n",
  );
  text = text.replace(
    '_struct_asym.id A\n_struct_asym.entity_id 1\n',
    'loop_\n_struct_asym.id\n_struct_asym.entity_id\nA 1\nB 2\nC 2\nD 2\n',
  );
  const atoms = [
    'HETATM 9 ZN ZN . ZN B 2 . ? 40.000 0.000 0.000 1.00 10.00 101 ZN Z ZN 1',
    'HETATM 10 ZN ZN . ZN C 2 . ? 45.000 0.000 0.000 1.00 10.00 102 ZN Z ZN 1',
    'HETATM 11 ZN ZN . ZN D 2 . ? 50.000 0.000 0.000 1.00 10.00 103 ZN Z ZN 1',
  ];
  text = `${text.replace(/#\s*$/, '')}${atoms.join('\n')}\n#\n`;
  return new TextEncoder().encode(text);
}

/** ins1.cif with a ligand, chain B: its heavy atoms first, its hydrogens after them (as wwPDB lists them). */
function withLigand(): Uint8Array {
  let text = readFileSync(new URL('./fixtures/structures/ins1.cif', import.meta.url), 'utf8');
  text = text.replace(
    "1 polymer 'Heavy chain fragment'\n",
    "1 polymer 'Heavy chain fragment'\n2 non-polymer 'LIGAND'\n",
  );
  text = text.replace(
    '_struct_asym.id A\n_struct_asym.entity_id 1\n',
    'loop_\n_struct_asym.id\n_struct_asym.entity_id\nA 1\nB 2\n',
  );
  const atoms = [
    'HETATM 9 C C1 . LIG B 2 . ? 40.000 0.000 0.000 1.00 10.00 101 LIG L C1 1',
    'HETATM 10 O O1 . LIG B 2 . ? 41.400 0.000 0.000 1.00 10.00 101 LIG L O1 1',
    'HETATM 11 H H1 . LIG B 2 . ? 39.500 0.900 0.000 1.00 10.00 101 LIG L H1 1',
    'HETATM 12 H H2 . LIG B 2 . ? 39.500 -0.900 0.000 1.00 10.00 101 LIG L H2 1',
  ];
  text = `${text.replace(/#\s*$/, '')}${atoms.join('\n')}\n#\n`;
  return new TextEncoder().encode(text);
}

async function linkOf(bytes: Uint8Array, ref: string) {
  const b = createSeqViewSpec();
  b.data.complex({
    ref,
    provider: 'molstar:structure',
    source: { uri: `https://example.org/${ref}.cif`, format: 'mmcif' },
  } as never);
  b.view({ ref: 'v', label: 'V' }).track({ label: '' }).ruler();
  const c = await compile(b.getDocument(), registries(), { load: async () => bytes } as never);
  expect(c.issues.filter((i) => i.severity === 'error')).toEqual([]);
  const structure = await structureOf(await readStructure(bytes, 'mmcif'), 0, null);
  return StructureLink.attach(new Translator(), '3d', structure, c.set, ref);
}

describe('wholeChains', () => {
  it("takes a ligand's hydrogens after its last heavy atom (perf-audit-2 10/B5)", async () => {
    const link = await linkOf(withLigand(), 'lig');
    const B = link.chains.find((m) => m.chain.labelAsymId === 'B')!;
    // two heavy atoms are its elements
    expect(B.chain.count).toBe(2);

    const runs = Extents.create(1);
    Extents.add(runs, B.system, 0, 2);
    const names: string[] = [];
    StructureElement.Loci.forEachLocation(wholeChains(link, runs), (l) => {
      names.push(l.unit.model.atomicHierarchy.atoms.label_atom_id.value(l.element));
    });
    expect(names).toEqual(['C1', 'O1', 'H1', 'H2']);
  });

  it('lights the touched ions of a merged unit, not the one between them', async () => {
    const link = await linkOf(withIons(), 'zn');
    const ion = (asym: string) => link.chains.find((m) => m.chain.labelAsymId === asym)!;
    const [B, C, D] = [ion('B'), ion('C'), ion('D')];
    expect(B.chain.unit).toBe(D.chain.unit);
    expect(C.chain.unit).toBe(B.chain.unit);

    const runs = Extents.create(2);
    Extents.add(runs, B.system, 0, 1);
    Extents.add(runs, D.system, 0, 1);
    const asyms = new Set<string>();
    StructureElement.Loci.forEachLocation(wholeChains(link, runs), (l) => {
      const h = l.unit.model.atomicHierarchy;
      asyms.add(h.chains.label_asym_id.value(h.chainAtomSegments.index[l.element]));
    });
    expect([...asyms].sort()).toEqual(['B', 'D']);
  });
});
