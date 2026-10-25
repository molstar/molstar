/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

// The link's whole-residue reaction (seqspec:libraries#the-mol-adapter): `linkViews`' hover and selection light a
// non-polymer position's whole residue in 3D (`LinkReaction.wholeResidue`, default true), a ligand as a whole and a
// glycan's one monosaccharide; with it off, the atoms the positions name. `StructureLink.toLoci` itself, and so a
// SeqQL result back in 3D (`selectionLoci`), is exact. A polymer's position is its residue either way.

import { readFileSync } from 'node:fs';
import { createSeqViewSpec } from '@molstar/seq-core/builder';
import { compile, type SeqViewSpecDoc } from '@molstar/seq-core/format';
import { Extents, type SeqSet } from '@molstar/seq-core/model';
import { compileQuery, createQueryEnv, SeqSelection } from '@molstar/seq-core/model/query';
import { bindQuery, readQuery } from '@molstar/seq-core/query';
import { Translator } from '@molstar/seq-core/translate';
import { createSeqViewer, DefaultRegistry, DefaultSeqViewerSpec, manualScheduler } from '@molstar/seq-viewer';
import { StructureElement } from '@molstar/model/model/structure';
import { describe, expect, it } from 'vitest';
import { fixtureMeasurer } from '@molstar/seq-viewer/testing';
import {
  MolstarSvsEntry,
  molstarQueryExtension,
  readStructure,
  StructureLink,
  selectionLoci,
  structureOf,
} from '../index.js';
import { type LinkViewsOptions, linkViews } from '../link-views.js';
import { registries } from './fixtures.js';

/**
 * ins1.cif (a nine-residue chain A, a CA per residue; VAL 2 also a CB) with a ligand B (LIG: C1, C2, O1 in two
 * alternate locations, a hydrogen) and a glycan C (branched: NAG 1 with C1, O1 and a hydrogen, NAG 2 with C1, O1).
 */
function mixed(): Uint8Array {
  let text = readFileSync(new URL('./fixtures/structures/ins1.cif', import.meta.url), 'utf8').replace(/\r\n/g, '\n');
  text = text.replace(
    "1 polymer 'Heavy chain fragment'\n",
    "1 polymer 'Heavy chain fragment'\n2 non-polymer 'A small ligand'\n3 branched 'A glycan'\n",
  );
  text = text.replace(
    '_struct_asym.id A\n_struct_asym.entity_id 1\n',
    'loop_\n_struct_asym.id\n_struct_asym.entity_id\nA 1\nB 2\nC 3\n',
  );
  text = text.replace(
    'ATOM 2 C CA . VAL A 1 2 ? 3.800 0.000 0.000 1.00 20.00 51 VAL H CA 1\n',
    'ATOM 2 C CA . VAL A 1 2 ? 3.800 0.000 0.000 1.00 20.00 51 VAL H CA 1\n' +
      'ATOM 20 C CB . VAL A 1 2 ? 3.800 1.500 0.000 1.00 20.00 51 VAL H CB 1\n',
  );
  const atoms = [
    'HETATM 9 C C1 . LIG B 2 . ? 40.000 0.000 0.000 1.00 10.00 101 LIG L C1 1',
    'HETATM 10 C C2 . LIG B 2 . ? 41.500 0.000 0.000 1.00 10.00 101 LIG L C2 1',
    'HETATM 11 O O1 A LIG B 2 . ? 43.000 0.000 0.000 0.50 10.00 101 LIG L O1 1',
    'HETATM 12 O O1 B LIG B 2 . ? 43.000 0.500 0.000 0.50 10.00 101 LIG L O1 1',
    'HETATM 13 H H1 . LIG B 2 . ? 40.000 1.000 0.000 1.00 10.00 101 LIG L H1 1',
    'HETATM 14 C C1 . NAG C 3 . ? 60.000 0.000 0.000 1.00 10.00 201 NAG G C1 1',
    'HETATM 15 O O1 . NAG C 3 . ? 61.400 0.000 0.000 1.00 10.00 201 NAG G O1 1',
    'HETATM 16 H H1 . NAG C 3 . ? 60.000 1.000 0.000 1.00 10.00 201 NAG G H1 1',
    'HETATM 17 C C1 . NAG C 3 . ? 62.800 0.000 0.000 1.00 10.00 202 NAG G C1 1',
    'HETATM 18 O O1 . NAG C 3 . ? 64.200 0.000 0.000 1.00 10.00 202 NAG G O1 1',
  ];
  text = `${text.replace(/#\s*$/, '')}${atoms.join('\n')}\n#\n`;
  return new TextEncoder().encode(text);
}

function doc(): SeqViewSpecDoc {
  const b = createSeqViewSpec();
  b.data.complex({
    ref: 'ab',
    provider: 'molstar:structure',
    source: { uri: 'https://example.org/mixed.cif', format: 'mmcif' },
  } as never);
  b.view({ ref: 'v', label: 'V' }).track({ label: '' }).ruler();
  return b.getDocument();
}

/** The atoms of a loci as `<asym> <auth_seq_id> <atom><alt>`, sorted. */
function atoms(loci: StructureElement.Loci | null): string[] {
  const out: string[] = [];
  if (loci)
    StructureElement.Loci.forEachLocation(loci, (l) => {
      const h = l.unit.model.atomicHierarchy;
      const r = h.residueAtomSegments.index[l.element];
      const asym = h.chains.label_asym_id.value(h.chainAtomSegments.index[l.element]);
      const alt = h.atoms.label_alt_id.value(l.element);
      out.push(`${asym} ${h.residues.auth_seq_id.value(r)} ${h.atoms.label_atom_id.value(l.element)}${alt}`);
    });
  return out.sort();
}

async function linked() {
  const bytes = mixed();
  const c = await compile(doc(), registries(), { load: async () => bytes } as never);
  expect(c.issues.filter((i) => i.severity === 'error')).toEqual([]);
  const structure = await structureOf(await readStructure(bytes, 'mmcif'), 0, null);
  const link = StructureLink.attach(new Translator(), '3d', structure, c.set, 'ab');
  const member = (kind: string) => {
    const m = link.chains.find((x) => x.chain.kind === kind);
    if (!m) throw new Error(`no ${kind} chain`);
    return m;
  };
  const at = (kind: string, start: number, stop: number, o?: { wholeResidue?: boolean }) => {
    const runs = Extents.create(1);
    Extents.add(runs, member(kind).system, start, stop);
    return atoms(link.toLoci(runs, o));
  };
  return { set: c.set, link, member, at };
}

const LIGAND = ['B 101 C1', 'B 101 C2', 'B 101 H1', 'B 101 O1A', 'B 101 O1B'];

describe("StructureLink.toLoci's wholeResidue", () => {
  it('a ligand position: its whole residue with the option (hydrogens, both alt locs), its atom without', async () => {
    const { member, at } = await linked();
    // Heavy atoms, one per name: C1, C2, O1.
    expect(member('ligand').chain.count).toBe(3);
    expect(at('ligand', 0, 1)).toEqual(['B 101 C1']);
    expect(at('ligand', 2, 3)).toEqual(['B 101 O1A']);
    expect(at('ligand', 0, 1, { wholeResidue: true })).toEqual(LIGAND);
    expect(at('ligand', 0, 3, { wholeResidue: true })).toEqual(LIGAND);
  });

  it("a branched chain's position lights its monosaccharide, not the whole glycan", async () => {
    const { member, at } = await linked();
    expect(member('branched').chain.count).toBe(4);
    expect(at('branched', 0, 1)).toEqual(['C 201 C1']);
    expect(at('branched', 0, 1, { wholeResidue: true })).toEqual(['C 201 C1', 'C 201 H1', 'C 201 O1']);
    expect(at('branched', 2, 3, { wholeResidue: true })).toEqual(['C 202 C1', 'C 202 O1']);
    expect(at('branched', 1, 3, { wholeResidue: true })).toEqual([
      'C 201 C1',
      'C 201 H1',
      'C 201 O1',
      'C 202 C1',
      'C 202 O1',
    ]);
  });

  it("a polymer position is its residue's atoms either way", async () => {
    const { at } = await linked();
    expect(at('polymer', 1, 2)).toEqual(['A 51 CA', 'A 51 CB']);
    expect(at('polymer', 1, 2, { wholeResidue: true })).toEqual(['A 51 CA', 'A 51 CB']);
  });
});

describe('a SeqQL result back in 3D is exact', () => {
  it('selectionLoci of a ligand unit, and of a molql result, names only their atoms', async () => {
    const { set, link } = await linked();
    const env = createQueryEnv([molstarQueryExtension(() => link)]);
    const run = (text: string) => {
      const bound = bindQuery(readQuery(text).node as never, env.language);
      if (!bound.expr) throw new Error(`cannot bind ${text}`);
      return compileQuery(bound.expr, env, bound.source).run({ set, inputs: {} }).selection as SeqSelection;
    };
    const unit = link.chains.find((x) => x.chain.kind === 'ligand')?.unitRef as string;
    // The ligand's unit: its three heavy-atom positions, no hydrogen, no second alternate location.
    expect(atoms(selectionLoci(link, run(unit)))).toEqual(['B 101 C1', 'B 101 C2', 'B 101 O1A']);
    // Mol*'s own selection of one ligand atom, through Seq* and back: that atom.
    expect(atoms(selectionLoci(link, run('(molql "resn LIG and name C2")')))).toEqual(['B 101 C2']);
  });
});

describe("linkViews' hover and selection", () => {
  async function viewer(options: LinkViewsOptions = {}) {
    const bytes = mixed();
    const scheduler = manualScheduler();
    const ctx = createSeqViewer(
      { ...DefaultSeqViewerSpec(), registry: [DefaultRegistry, { svs: MolstarSvsEntry }] },
      { scheduler, loaders: { load: async () => bytes } as never, text: fixtureMeasurer() },
    );
    const done = ctx.setDocument(doc());
    for (const ref of ctx.viewports.sections.value) ctx.viewports.get(ref)?.setSize(800, 400, 1);
    expect((await done).ok).toBe(true);
    await ctx.graph.update();
    scheduler.advance(50);
    const structure = await structureOf(await readStructure(bytes, 'mmcif'), 0, null);
    const link = StructureLink.attach(ctx.systems, '3d', structure, ctx.graph.assembled.value?.set as SeqSet, 'ab');
    const lit: string[][] = [];
    const selected: string[][] = [];
    const unlink = linkViews(
      ctx,
      {
        highlight: (l) => lit.push(atoms(l)),
        select: (l) => selected.push(atoms(l)),
        focus: () => {},
        onHover: () => () => {},
        onSelect: () => () => {},
      },
      link,
      options,
    );
    const unitOf = (kind: string) => link.chains.find((x) => x.chain.kind === kind)?.unitRef as string;
    const hover = (kind: string, start: number) => {
      ctx.hover.setLoci(ctx.refs.toLoci(unitOf(kind), { start, stop: start + 1 }), 'pointer');
      return lit.at(-1);
    };
    const select = (kind: string, start: number) => {
      ctx.selection.set(ctx.refs.toLoci(unitOf(kind), { start, stop: start + 1 }), 'pointer');
      return selected.at(-1);
    };
    const close = () => {
      unlink();
      link.detach();
      ctx.dispose();
    };
    return { hover, select, close };
  }

  it('light and select the whole ligand by default; a polymer residue stays its residue', async () => {
    const { hover, select, close } = await viewer();
    expect(hover('ligand', 1)).toEqual(LIGAND);
    expect(select('ligand', 2)).toEqual(LIGAND);
    expect(hover('branched', 3)).toEqual(['C 202 C1', 'C 202 O1']);
    expect(hover('polymer', 1)).toEqual(['A 51 CA', 'A 51 CB']);
    expect(select('polymer', 1)).toEqual(['A 51 CA', 'A 51 CB']);
    close();
  }, 30_000);

  it('with wholeResidue off, the atoms the positions name', async () => {
    const { hover, select, close } = await viewer({ hover: { wholeResidue: false }, select: { wholeResidue: false } });
    expect(hover('ligand', 1)).toEqual(['B 101 C2']);
    expect(select('ligand', 2)).toEqual(['B 101 O1A']);
    expect(hover('branched', 3)).toEqual(['C 202 O1']);
    expect(hover('polymer', 1)).toEqual(['A 51 CA', 'A 51 CB']);
    close();
  }, 30_000);
});
