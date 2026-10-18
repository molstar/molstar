/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// SeqQL over Mol* (query.ts; seqspec:libraries#the-mol-adapter): `molql` on 1TUP (p53's DNA-binding domain on DNA, a
// zinc in each p53 chain) onto P04637, Mol*'s selection as the `structure` input, a result back as a loci, and the
// errors.

import { compile } from '@molstar/seq-core/format';
import { Extents, type SeqSet } from '@molstar/seq-core/model';
import { compileQuery, createQueryEnv, type QueryEnv, SeqSelection } from '@molstar/seq-core/model/query';
import { bindQuery, readQuery, renderIssue } from '@molstar/seq-core/query';
import { Translator } from '@molstar/seq-core/translate';
import { StructureElement } from '@molstar/model/model/structure';
import { beforeAll, describe, expect, it } from 'vitest';
import {
  MOLSTAR_QUERY_SYMBOLS,
  molstarQueryExtension,
  readStructure,
  StructureLink,
  selectionLoci,
  structureInput,
  structureOf,
  structureReach,
} from '../index.js';
import { bytesOf, loaders, p53Doc, registries } from './fixtures.js';

let set: SeqSet;
let link: StructureLink;
let shown: StructureLink | null;
let env: QueryEnv;

beforeAll(async () => {
  const r = await compile(p53Doc([{ ref: '1tup', pdb: '1TUP', entity: '3' }]), registries(), loaders());
  set = r.set;
  const structure = await structureOf(await readStructure(bytesOf('1tup.bcif'), 'bcif'), 0, null);
  link = StructureLink.attach(new Translator(), '3d', structure, set, '1tup');
  env = createQueryEnv([molstarQueryExtension(() => shown)]);
});

function run(text: string, o: { env?: QueryEnv; structure?: SeqSelection } = {}) {
  const e = o.env ?? env;
  const read = readQuery(text);
  const bound = bindQuery(read.node as never, e.language);
  if (!bound.expr) return { selection: null, said: bound.issues.map((i) => renderIssue(i).split('\n')[0]) };
  const r = compileQuery(bound.expr, e, bound.source).run({
    set,
    inputs: o.structure ? { structure: o.structure } : {},
  });
  return { selection: r.selection, said: r.issues.map((i) => renderIssue(i).split('\n')[0]) };
}
const summary = (text: string) => {
  shown = link;
  const r = run(text);
  expect(r.said).toEqual([]);
  return SeqSelection.summary(r.selection as SeqSelection, { maxParts: 20 });
};

describe('molstar.select (molql)', () => {
  it("the zinc site: Mol*'s residues within 3 Å of each zinc (C176, H179, C238, C242, and N239's backbone), on P04637 through SIFTS", () => {
    expect(summary('(onto (molql "polymer and (resn ZN around 3)") p53)')).toBe(
      '1 group, 5 positions: p53/unit 176, 179, 238-239, 242',
    );
  });

  it('the residues touching DNA, on the chains and on P04637', () => {
    const onP53 = summary('(onto (molql "polymer.protein and (polymer.nucleic around 4)") p53)');
    // K120 in the major groove, R248 in the minor, R273 and R280 on the backbone and the bases
    expect(onP53).toBe(
      '1 group, 24 positions: p53/unit 101, 118-122, 165, 167-168, 239-241, 243, 247-249, 273-277, 280-281, 283',
    );
    // and on the three p53 chains
    expect(summary('(molql "polymer.protein and (polymer.nucleic around 4)")')).toMatch(
      /^1 group, 39 positions: 1tup\/C 101, .*; 1tup\/D 118-122, .*; 1tup\/E 119-121, /,
    );
  });

  it('reads MolScript with :lang', () => {
    const text = '(sel.atom.atom-groups :residue-test (= atom.label_comp_id ZN))';
    expect(summary(`(molql "${text}" :lang mol-script)`)).toBe('1 group, 3 positions: 1tup/F ZN; 1tup/G ZN; 1tup/H ZN');
  });

  it('an unreadable selection, nothing selected, no structure, no runtime, no extension', () => {
    shown = link;
    expect(run('(molql "resn ZN aroundd 3")').said[0]).toMatch(/^1:8 molstar.select.syntax: Mol\* cannot read/);
    expect(run('(molql "resn HEM")').said).toEqual([
      '1:1 seqql.empty: (molql "resn HEM") selected nothing: none in 1TUP',
    ]);
    shown = null;
    expect(run('(molql "polymer")').said).toEqual([
      '1:1 seqql.symbol.unavailable: molql needs a structure shown in 3D over these sequences; it runs where a structure is shown in 3D and linked',
    ]);
    expect(run('(molql "polymer")', { env: createQueryEnv([MOLSTAR_QUERY_SYMBOLS]) }).said[0]).toMatch(
      /seqql.symbol.unavailable: molql cannot run here: this host has no runtime for molstar.select/,
    );
    expect(run('(molql "polymer")', { env: createQueryEnv() }).said[0]).toMatch(/seqql.symbol.unknown/);
  });
});

describe('the structure input and back', () => {
  it("Mol*'s selection is (input structure); a result goes back to 3D as a loci", () => {
    shown = link;
    const c = link.chains.find((x) => x.unitRef === '1tup/C');
    const runs = Extents.create(1);
    Extents.add(runs, c?.unit as never, 81, 82);
    const picked = link.toLoci(runs);
    const structure = structureInput(link, picked);
    expect(SeqSelection.summary(structure)).toBe('1 group, 1 position: 1tup/C 175');
    const r = run('(onto (input structure) p53)', { structure });
    expect(SeqSelection.summary(r.selection as SeqSelection)).toBe('1 group, 1 position: p53/unit 175');

    const zn = run('(onto (molql "polymer and (resn ZN around 3)") p53)').selection as SeqSelection;
    const back = selectionLoci(link, zn);
    const residues = new Set<string>();
    StructureElement.Loci.forEachLocation(back, (l) => {
      const h = l.unit.model.atomicHierarchy;
      residues.add(String(h.residues.auth_seq_id.value(h.residueAtomSegments.index[l.element])));
    });
    expect([...residues].sort()).toEqual(['176', '179', '238', '239', '242']);
  });
});

describe('instance space (QM1b)', () => {
  it('a result on the chains stays on them in 3D; one on P04637 reaches every chain', () => {
    shown = link;
    const reach = (text: string) => structureReach(link, run(text).selection as SeqSelection);
    expect(reach('(molql "polymer and (resn ZN around 3)")')).toEqual({ chains: ['A', 'B', 'C'], exact: true });
    // the DNA contacts of chain A alone stay chain A's; onto P04637 they are any chain's
    const a = '(and 1tup/C (molql "polymer.protein and (polymer.nucleic around 4)"))';
    expect(reach(a)).toEqual({ chains: ['A'], exact: true });
    expect(reach(`(onto ${a} p53)`)).toEqual({ chains: ['A', 'B', 'C'], exact: false });
    const d = link.chains.find((x) => x.unitRef === '1tup/D');
    const runs = Extents.create(1);
    Extents.add(runs, d?.unit as never, 81, 82);
    expect(SeqSelection.summary(structureInput(link, link.toLoci(runs)))).toBe('1 group, 1 position: 1tup/D 175');
  });
});
