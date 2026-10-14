/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The link (seqspec:libraries#the-mol-adapter): a structure's chains as external systems over a Translator, Seq* runs
// to Mol* loci and back, through P04637's SIFTS mapping, an assembly's copies onto one deposited unit, and detach;
// pLDDT both ways
//.

import { compile } from '@molstar/seq-core/format';
import { type CoordinateSystem, CoordinateSystemKind, Extents, KeyKind, type SeqSet } from '@molstar/seq-core/model';
import { Translator } from '@molstar/seq-core/translate';
import { type Structure, StructureElement } from '@molstar/model/model/structure';
import { describe, expect, it } from 'vitest';
import { readStructure, StructureLink, structureOf } from '../src/index.js';
import { bytesOf, loaders, p53Doc, registries } from './fixtures.js';

const UNIT = CoordinateSystemKind.unit;
function unitOf(set: SeqSet, ref: string): CoordinateSystem {
  const t = { kind: 0, id: 0, row: -1 };
  if (!set.byName(ref, t) || t.kind !== KeyKind.unit) throw new Error(`no unit ${ref}`);
  return { kind: UNIT, id: t.id, row: -1 };
}
const runs = (...rs: [CoordinateSystem, number, number][]) => {
  const e = Extents.create(4);
  for (const [s, a, b] of rs) Extents.add(e, s, a, b);
  return e;
};
/** Runs as [unit ref, start, stop], sorted. */
const named = (set: SeqSet, e: Extents) =>
  Array.from({ length: e.count }, (_, i): [string, number, number] => [
    e.kind[i] === UNIT ? set.strings.get(set.units.key[e.id[i]]) : `${e.kind[i]}:${e.id[i]}`,
    e.start[i],
    e.stop[i],
  ]).sort((a, b) => a[0].localeCompare(b[0]) || a[1] - b[1]);
/** A Mol* loci's residues as `<label_asym_id> <auth_seq_id> <comp>`, and its atom count. */
function residues(loci: StructureElement.Loci) {
  const out = new Set<string>();
  let atoms = 0;
  StructureElement.Loci.forEachLocation(loci, (l) => {
    const h = l.unit.model.atomicHierarchy;
    const chain = h.chains.label_asym_id.value(h.chainAtomSegments.index[l.element]);
    const r = h.residueAtomSegments.index[l.element];
    out.add(`${chain} ${h.residues.auth_seq_id.value(r)} ${h.atoms.label_comp_id.value(l.element)}`);
    atoms++;
  });
  return { residues: [...out].sort(), atoms };
}

async function linked(pdb: string, entity: string, assembly: string | null = null) {
  const ref = pdb.toLowerCase();
  const { set, issues } = await compile(p53Doc([{ ref, pdb, entity }]), registries(), loaders());
  expect(issues.filter((i) => i.severity === 'error')).toEqual([]);
  const structure: Structure = await structureOf(await readStructure(bytesOf(`${ref}.bcif`), 'bcif'), 0, assembly);
  const t = new Translator();
  return { set, t, link: StructureLink.attach(t, '3d-1', structure, set, ref) };
}

describe('StructureLink', () => {
  it('registers the structure and its chains, each mapped onto its Seq* unit', async () => {
    const { t, link } = await linked('1TUP', '3');
    expect(link.key).toBe('molstar:3d-1');
    expect(link.chains.map((c) => [c.key.slice('molstar:3d-1/'.length), c.unitRef, c.chain.count])).toEqual([
      ['A@1_555', '1tup/A', 21],
      ['B@1_555', '1tup/B', 21],
      ['C@1_555', '1tup/C', 196],
      ['D@1_555', '1tup/D', 194],
      ['E@1_555', '1tup/E', 195],
      ['F@1_555', '1tup/F', 1],
      ['G@1_555', '1tup/G', 1],
      ['H@1_555', '1tup/H', 1],
    ]);
    expect(t.descriptor('molstar:3d-1')?.members?.length).toBe(8);
    expect(t.descriptor('molstar:3d-1/C@1_555')?.structure).toMatchObject({ labelAsymId: 'C', authAsymId: 'A' });
  });

  it('Seq* runs to Mol* loci and back: a residue, a chain, an ion, an unobserved stretch', async () => {
    const { set, link } = await linked('1TUP', '3');
    const c = unitOf(set, '1tup/C');
    const r175 = link.toLoci(runs([c, 81, 82]));
    expect(residues(r175).residues).toEqual(['C 175 ARG']);
    expect(named(set, link.seqRunsOf(r175, Extents.create(1)))).toEqual([['1tup/C', 81, 82]]);
    const chain = link.toLoci(runs([c, 0, 196]));
    expect(residues(chain).residues.length).toBe(196);
    expect(named(set, link.seqRunsOf(chain, Extents.create(1)))).toEqual([['1tup/C', 0, 196]]);
    const zn = link.toLoci(runs([unitOf(set, '1tup/F'), 0, 1]));
    expect(residues(zn)).toEqual({ residues: ['F 951 ZN'], atoms: 1 });
    expect(named(set, link.seqRunsOf(zn, Extents.create(1)))).toEqual([['1tup/F', 0, 1]]);
    expect(link.toLoci(runs([c, 200, 210])).elements.length).toBe(0);
  });

  it('a residue of P04637 reaches every copy through the SIFTS mapping, and comes back to each', async () => {
    const { set, link } = await linked('1TUP', '3');
    const loci = link.toLoci(runs([unitOf(set, 'p53/unit'), 174, 175]));
    expect(residues(loci).residues).toEqual(['C 175 ARG', 'D 175 ARG', 'E 175 ARG']);
    expect(named(set, link.seqRunsOf(loci, Extents.create(1)))).toEqual([
      ['1tup/C', 81, 82],
      ['1tup/D', 81, 82],
      ['1tup/E', 81, 82],
    ]);
  });

  it("an assembly's copies of a deposited chain all pair with its one unit", async () => {
    const { set, link } = await linked('1C26', '1', '1');
    expect(link.chains.length).toBe(4);
    expect(link.chains.every((c) => c.unitRef === '1c26/A')).toBe(true);
    const loci = link.toLoci(runs([unitOf(set, '1c26/A'), 10, 11]));
    expect(loci.elements.length).toBe(4);
    expect(named(set, link.seqRunsOf(loci, Extents.create(1)))).toEqual([['1c26/A', 10, 11]]);
  });

  it('detach unregisters; attaching again names the same systems', async () => {
    const { set, t, link } = await linked('2OCJ', '1');
    const ids = link.chains.map((c) => c.system.id);
    link.detach();
    expect(t.descriptor(link.key)).toBe(null);
    const out = Extents.create(1);
    t.translate(runs([link.chains[0].system, 0, 5]), link.chains[0].system, unitOf(set, '2ocj/A'), out);
    expect(out.count).toBe(0);
    const again = StructureLink.attach(t, '3d-1', link.structure, set, '2ocj');
    expect(again.chains.map((c) => c.system.id)).toEqual(ids);
  });
});

describe('pLDDT', () => {
  it("molstar:plddt reads the B-factors; the model's ma_qa_metric_local category gives the same values", async () => {
    const { set, issues } = await compile(p53Doc([{ ref: 'af' }]), registries(), loaders());
    expect(issues.filter((i) => i.severity === 'error')).toEqual([]);
    const byKey = (k: string) => set.annotations.find((a) => a && set.strings.get(a.key) === k);
    const [b, qa] = [byKey('af-plddt'), byKey('af-qa')];
    expect([b?.count, qa?.count]).toEqual([393, 393]);
    const values = (a: typeof b) => {
      const v = new Float64Array(393).fill(Number.NaN);
      for (let r = 0; a && r < a.count; r++) v[a.start[r]] = a.sparse?.[r] ?? Number.NaN;
      return Array.from(v);
    };
    expect(values(b)).toEqual(values(qa));
    expect(values(b).every((x) => x >= 0 && x <= 100)).toBe(true);
    expect(set.strings.get(b?.semanticType ?? -1)).toBe('plddt');
  });
});
