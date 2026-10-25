/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

// structureMvs's `ligandComps` (seqspec:libraries#the-mol-adapter): the named ligands as ball and stick, never a
// polymer's residue of the same component (a modified residue) and never an ion, which the ion component draws once.

import { readFileSync } from 'node:fs';
import { MVSData } from '@molstar/mvs-builder/mvs-data';
import { QueryContext, StructureElement, StructureSelection } from '@molstar/model/model/structure';
import { compile } from '@molstar/model/script/runtime/query/compiler';
import { MolScriptBuilder as MS } from '@molstar/query-language/language/builder';
import type { Expression } from '@molstar/query-language/language/expression';
import { describe, expect, it } from 'vitest';
import { readStructure, structureOf } from '@molstar/prism-viewer/sequence/molstar';
import { structureMvs } from '../../src/_stories/support/seqstar-mvs.js';

/** ins1.cif with its first residue a selenomethionine (MSE, a modified residue), a ligand B (LIG) and a zinc C. */
function withLigandAndIon(): Uint8Array {
  let text = readFileSync(
    new URL('../../../../prism/viewer/src/sequence/molstar/_test/fixtures/structures/ins1.cif', import.meta.url),
    'utf8',
  ).replace(/\r\n/g, '\n');
  text = text
    .replace('1 1 TRP n\n', '1 1 MSE n\n')
    .replace(' TRP A 1 1 ', ' MSE A 1 1 ')
    .replace(' TRP H CA ', ' MSE H CA ');
  text = text.replace(
    "1 polymer 'Heavy chain fragment'\n",
    "1 polymer 'Heavy chain fragment'\n2 non-polymer 'A small ligand'\n3 non-polymer 'ZINC ION'\n",
  );
  text = text.replace(
    '_struct_asym.id A\n_struct_asym.entity_id 1\n',
    'loop_\n_struct_asym.id\n_struct_asym.entity_id\nA 1\nB 2\nC 3\n',
  );
  const atoms = [
    'HETATM 9 C C1 . LIG B 2 . ? 40.000 0.000 0.000 1.00 10.00 101 LIG L C1 1',
    'HETATM 10 O O1 . LIG B 2 . ? 41.400 0.000 0.000 1.00 10.00 101 LIG L O1 1',
    'HETATM 11 ZN ZN . ZN C 3 . ? 50.000 0.000 0.000 1.00 10.00 201 ZN L ZN 1',
  ];
  text = `${text.replace(/#\s*$/, '')}${atoms.join('\n')}\n#\n`;
  return new TextEncoder().encode(text);
}

type Node = { kind: string; params?: Record<string, unknown>; children?: Node[] };
const components = (n: Node, out: Node[] = []): Node[] => {
  if (n.kind === 'component') out.push(n);
  for (const c of n.children ?? []) components(c, out);
  return out;
};

describe("structureMvs's ligandComps", () => {
  it('draws the named ligands only: not a modified residue of the polymer, not an ion a second time', async () => {
    const structure = await structureOf(await readStructure(withLigandAndIon(), 'mmcif'), 0, null);
    const state = structureMvs({
      url: 'https://example.org/x.cif',
      format: 'mmcif',
      ligandComps: ['LIG', 'ZN', 'MSE'],
    });
    expect(MVSData.validationIssues(state)).toBeUndefined();
    const selectors = components(state.root as Node).map((c) => c.params?.selector);
    // The polymer, the named ligands (one MolQL selector), the ions.
    expect(selectors.filter((s) => s === 'ion')).toHaveLength(1);
    const named = selectors.filter(
      (s): s is { molql: Expression } => typeof s === 'object' && s !== null && 'molql' in s,
    );
    expect(named).toHaveLength(1);
    const compsOf = (molql: Expression) => {
      const picked = compile<StructureSelection>(molql)(new QueryContext(structure));
      const comps = new Set<string>();
      StructureElement.Loci.forEachLocation(StructureSelection.toLociWithSourceUnits(picked), (l) => {
        comps.add(l.unit.model.atomicHierarchy.atoms.label_comp_id.value(l.element));
      });
      return [...comps].sort();
    };
    expect(compsOf(named[0].molql)).toEqual(['LIG']);
    // The control: by component id alone (what MVS's `{ label_comp_id }` selects), all three.
    const byComp = MS.struct.generator.atomGroups({
      'residue-test': MS.core.set.has([MS.set('LIG', 'ZN', 'MSE'), MS.ammp('label_comp_id')]),
    });
    expect(compsOf(byComp)).toEqual(['LIG', 'MSE', 'ZN']);
  });

  it('draws every ligand without the list, and none with an empty one; ions either way', () => {
    const selectorsOf = (ligandComps?: string[]) =>
      components(structureMvs({ url: 'u', format: 'mmcif', ...(ligandComps ? { ligandComps } : {}) }).root as Node).map(
        (c) => c.params?.selector,
      );
    expect(selectorsOf()).toEqual(['polymer', 'ligand', 'ion']);
    expect(selectorsOf([])).toEqual(['polymer', 'ion']);
  });
});
