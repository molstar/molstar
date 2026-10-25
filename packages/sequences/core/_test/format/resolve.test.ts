/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { createSeqViewSpec } from '@molstar/seq-core/builder';
import {
  type ComposedPath,
  Coords,
  type Issue,
  KeyKind,
  type KeyTarget,
  Numbering,
  OrderedSet,
  Paths,
  type PositionLabel,
  type SeqSet,
} from '@molstar/seq-core/model';
import { describe, expect, it } from 'vitest';
import {
  applyPatch,
  ColorLists,
  compile,
  composeView,
  indexDoc,
  layerPalette,
  normalize,
  type Palette,
  type Patch,
  type ResolvedView,
  resolveLayer,
  resolveSection,
  resolveView,
  type SeqViewSpecDoc,
  snapshotDoc,
  splitRows,
  supportedAll,
  type ThemeRef,
  validate,
} from '../../src/format/index.js';
import { pathsSearched } from '../../src/format/resolve/reach.js';
import { pf00042 } from '../model/alignment-fixtures.js';
import { CASES, case3, case6, case7, case9, case10, case12, case13 } from './cases.js';
import { files, fixtureLoaders, registries } from './load-fixtures.js';

const reg = registries();
const all = supportedAll(reg);
const codes = (issues: Issue[]) => [...new Set(issues.map((i) => i.code))].sort();
const label = (set: SeqSet, id: number) => set.strings.get(id);
const t: KeyTarget = { kind: -1, id: -1, row: -1 };
const id = (set: SeqSet, key: string) => (set.byName(key, t) ? t.id : -1);
/** A path's from index whose image is `i` (the old `backward`), through its coordinate map. */
const back = (p: ComposedPath, i: number) => {
  const hit = { element: -1, phase: 0 };
  return Coords.forEachFromAt(p.map, i, hit, () => {}, null) > 0 ? hit.element : -1;
};

async function resolve(doc: SeqViewSpecDoc, view: string, supported: ReadonlySet<string> = all, r = reg, fs = files()) {
  const { set, issues } = await compile(doc, r, fixtureLoaders(fs));
  expect(issues.filter((i) => i.severity === 'error')).toEqual([]);
  return { set, ...resolveView(doc, set, view, r, supported) };
}
/** A layer's placements: [block, member, lane, path] per placement. */
function placed(v: ResolvedView, set: SeqSet, layer: string) {
  const l = Array.from(v.layers.ref).indexOf(set.strings.find(layer));
  const out: number[][] = [];
  for (let k = v.layers.placementStart[l]; k < v.layers.placementStart[l] + v.layers.placementCount[l]; k++) {
    out.push([v.placements.block[k], v.placements.member[k], v.placements.lane[k], v.placements.path[k]]);
  }
  return out;
}
const patched = (doc: SeqViewSpecDoc, ops: Patch['ops']) => {
  const r = applyPatch(doc, { format: 'seq_view_spec_patch', ops }, reg);
  expect(r.issues).toEqual([]);
  return r.doc;
};

describe('every case resolves', () => {
  it.each(Object.entries(CASES))('%s: each view, no errors', async (_, make) => {
    const doc = make();
    const { set } = await compile(doc, reg, fixtureLoaders());
    for (const v of indexDoc(normalize(doc, reg)).order.filter((n) => n.kind === 'view')) {
      const r = resolveView(doc, set, v.ref, reg, all);
      expect(r.issues.filter((i) => i.severity === 'error')).toEqual([]);
      expect(r.view.sections.count).toBeGreaterThan(0);
    }
  });
});

describe('resolving the cases', () => {
  it("7: UniProt features on six 1BRS chains through two mappings; copies merged at the first one's place", async () => {
    const doc = case7();
    const each = await resolve(doc, 'complex-each');
    expect(Array.from(each.view.blocks.label, (l) => label(each.set, l))).toEqual(['A', 'B', 'C', 'D', 'E', 'F']);
    const feats = placed(each.view, each.set, 'complex-each.ranges');
    expect(feats.map((p) => p[0])).toEqual([0, 1, 2, 3, 4, 5]);
    expect(feats.every((p) => p[1] === -1 && p[3] >= 0)).toBe(true);
    // Barnase P00648 48 is chain A's index 0: the path's forward carries the precursor into the chain.
    const path = each.view.paths[feats[0][3]];
    expect([Paths.at(path, 47), Paths.at(path, 46)]).toEqual([0, -1]);
    expect(codes(each.issues)).toEqual(['svs.layer.unplaced']); // the hydropathy provider is not registered

    const one = await resolve(doc, 'complex-one');
    expect(Array.from(one.view.blocks.label, (l) => label(one.set, l))).toEqual([
      'Barnase ×3 (A, B, C)',
      'Barstar ×3 (D, E, F)',
    ]);
    expect(one.view.blocks.members[0]).toEqual(Int32Array.from([0, 1, 2].map((k) => id(one.set, `1brs/${'ABC'[k]}`))));
    expect(placed(one.view, one.set, 'complex-one.ranges').map((p) => p.slice(0, 2))).toEqual([
      [0, -1],
      [1, -1],
    ]);
    // Unit-level data lands per member: chain A's secondary structure never reaches its copies.
    const ss = placed(one.view, one.set, 'complex-one.secondary');
    expect(ss.filter((p) => p[0] === 0).map((p) => p[1])).toEqual(['A', 'B', 'C'].map((c) => id(one.set, `1brs/${c}`)));
    expect(ss.every((p) => p[3] === -1)).toBe(true);
    // The variant effects matrix only on the barnase block.
    expect(placed(one.view, one.set, 'complex-one.matrix').map((p) => p[0])).toEqual([0]);
    expect(Array.from(one.view.groups.collapsed)).toEqual([0, 1]);
  });

  it("7 without the provider: the sequences' own units, nothing cascading", async () => {
    const r = await resolve(case7(), 'complex-one', all, registries(false));
    expect(Array.from(r.view.blocks.unit, (u) => label(r.set, r.set.units.key[u]))).toEqual([
      'barnase/unit',
      'barstar/unit',
    ]);
    expect(r.issues.filter((i) => i.severity === 'error')).toEqual([]);
  });

  it('12: a complex block, a copy group of two chains, the heme as one position, statistics on both chains', async () => {
    const r = await resolve(case12(), 'design-view');
    expect(Array.from(r.view.blocks.label, (l) => label(r.set, l))).toEqual(['Design 7 ×2 (A, B)', 'dimer.heme']);
    expect([Array.from(r.view.blocks.length), Array.from(r.view.blocks.ligand)]).toEqual([
      [55, 1],
      [0, 1],
    ]);
    const cons = placed(r.view, r.set, 'design-view.values');
    expect(cons.map((p) => p.slice(0, 2))).toEqual([[0, -1]]);
    const each = await resolve(
      patched(case12(), [{ op: 'set', ref: 'design-view.main', params: { copies: 'each' } }]),
      'design-view',
    );
    expect(each.view.blocks.count).toBe(3);
    expect(placed(each.view, each.set, 'design-view.values').map((p) => p[0])).toEqual([0, 1]);
  });

  it("3: members' sites on their rows' lanes in tree order; the filter keeps phosphosites", async () => {
    const r = await resolve(case3(), 'family');
    const order = Array.from(r.set.alignments[0]?.tree?.leafOrder ?? []);
    expect(Array.from(r.view.lanes.row)).toEqual(order);
    const sites = placed(r.view, r.set, 'family.points-2');
    expect(sites.map((p) => r.view.lanes.row[p[2]]).sort()).toEqual([0, 1, 2]);
    const l = Array.from(r.view.layers.ref).indexOf(r.set.strings.find('family.points-2'));
    expect(Array.from(r.view.kept[r.view.layers.dataStart[l]] ?? [])).toEqual([1, 1, 0]);
    expect(r.view.layers.domainMin[1]).toBeGreaterThanOrEqual(0); // the conservation layer's domain over its data
    expect(r.view.layers.domainMax[1]).toBeLessThanOrEqual(1);
  });

  it('6: a numbering through the path chain → 1a → 1b; a ligand block keeps its own; the pocket cuts the section', async () => {
    const doc = patched(case6(), [{ op: 'set', ref: 'chains.main', params: { numbering: 'uniprot:P00519-2' } }]);
    const r = await resolve(doc, 'chains');
    const b = r.view.blocks;
    expect([label(r.set, b.numbering[0]?.name ?? -1), label(r.set, b.numbering[1]?.name ?? -1)]).toEqual([
      'uniprot:P00519-2',
      'atom',
    ]);
    // Auth 315 is chain index 90 (the chain starts at 225); 1a 315; 1b 334 (its first exon is 19 longer).
    const path = r.view.paths[b.numberingPath[0]];
    const out: PositionLabel = { num: 0, suffix: -1 as never, name: -1 as never };
    expect(Numbering.labelOf(b.numbering[0] as Numbering, back(path, 90) as never, out) && out.num).toBe(334);
    expect(b.numberingPath[1]).toBe(-1);

    const poses = await resolve(case6(), 'poses-view');
    const shown = poses.view.blocks.positions[0];
    expect(shown && OrderedSet.has(shown, 314)).toBe(true); // T315 (1a index 314) lines the pocket
    expect(shown && OrderedSet.size(shown)).toBeLessThan(40);
  });

  it('10: the gatekeeper snapshot enters the strip, focused on auth 305–325', async () => {
    const { doc } = snapshotDoc(case10(), 's3', reg);
    const r = await resolve(doc, 'strip');
    expect(r.view.sections.focus[0]).toEqual({ block: 0, start: 80, stop: 101 });
    const byItems = patched(case10(), [
      {
        op: 'insert',
        parent: 'strip.main',
        node: { kind: 'focus', ref: 'f', params: { annotation: 'site.gatekeeper' } },
      },
    ]);
    expect((await resolve(byItems, 'strip')).view.sections.focus[0]).toEqual({ block: 0, start: 90, stop: 91 });
    const unknown = patched(case10(), [{ op: 'set', ref: 'strip.main', params: { numbering: 'kabat' } }]);
    expect(codes((await resolve(unknown, 'strip')).issues)).toEqual(['svs.numbering.fallback']);
    const both = await resolve(case10(), 'both');
    expect([both.view.links.count, both.issues]).toEqual([1, []]);
    const bad = patched(case10(), [{ op: 'set', ref: 'both.ann', params: { positions: 'msa-gatekeeper' } }]);
    expect(codes((await resolve(bad, 'both')).issues)).toEqual(['svs.section.scroll_group']);
  });
});

/** PF00042's 32-row seed as an aligned FASTA, HBA_HUMAN (row 0) bound to 1A3N's α entity, and α helices. */
function seed() {
  const rows = pf00042();
  const fasta = rows
    .map((r) => `>${r.label}/${r.start}-${r.start + r.gapped.replace(/-/g, '').length - 1}\n${r.gapped}\n`)
    .join('');
  const b = createSeqViewSpec({ title: 'PF00042 seed' });
  const aln = b.data.alignment({ ref: 'seed', source: { host_file: 'PF00042.seed.fasta', format: 'fasta' } });
  const hb = b.data.complex({
    ref: '1a3n',
    provider: 'molstar:structure',
    source: { service: 'pdbe', resource: 'structure', id: '1a3n' },
  });
  aln.bind(0, hb.sequence('entity:1'), { method: 'align', ref: 'seed-hba' });
  const helices = b.data.annotation({
    ref: 'helices',
    geometry: 'ranges',
    target: hb.unit('A'),
    numbering: 'auth',
    value: { beg: [3, 20, 52, 94], end: [18, 35, 71, 112], category: ['helix', 'helix', 'helix', 'helix'] },
  });
  const v = b.view({ ref: 'msa', label: 'Alignment' });
  v.block({ alignment: aln });
  v.track({ label: 'Helices (1A3N α)' }).secondary({ ref: 'msa.helices', data: helices });
  v.section({ ref: 'msa.chains', label: '1A3N' })
    .block({ unit: hb.unit('A') })
    .block({ unit: hb.unit('B') })
    .track({ label: 'Helices' })
    .secondary({ data: helices });
  return { doc: b.getDocument(), files: { ...files(), 'PF00042.seed.fasta': fasta } };
}

describe('placement across sections', () => {
  it("HBA helices on the PF00042 seed's columns, through the row's binding", async () => {
    const { doc, files: fs } = seed();
    const r = await resolve(doc, 'msa', all, reg, fs);
    expect(r.issues).toEqual([]);
    const [p] = placed(r.view, r.set, 'msa.helices');
    expect(p.slice(0, 3)).toEqual([0, -1, -1]);
    const path = r.view.paths[p[3]];
    const rows = r.set.alignments[0]?.rows;
    // Helix 52–71 (chain index 51–70) lands on row 0's columns: chain index i is row residue i - 25 (UniProt 27).
    const col = (i: number) => Paths.at(path, i);
    expect([col(51), col(70)]).toEqual([rows?.posMap(0 as never)?.[26], rows?.posMap(0 as never)?.[45]]);
    expect(col(10)).toBe(-1); // helix 3–18 lies before the row's first residue
  });

  it('a column reaches the α chain in the other section, β nothing (links)', async () => {
    const { doc, files: fs } = seed();
    const r = await resolve(doc, 'msa', all, reg, fs);
    const l = r.view.links;
    expect([l.count, l.fromBlock[0], l.toBlock[0]]).toEqual([1, 0, 1]);
    const path = r.view.paths[l.path[0]];
    const rows = r.set.alignments[0]?.rows;
    const column = rows?.posMap(0 as never)?.[30] as number;
    expect(Paths.at(path, column)).toBe(55);
    expect(back(path, 55)).toBe(column);
  });

  /** Case 9's chains section: per placement, its block's chain and its slot's chain (and whether through a path). */
  async function chains(doc: SeqViewSpecDoc) {
    const { set } = await compile(doc, reg, fixtureLoaders());
    const { view: v } = resolveView(doc, set, 'msa', reg, all);
    const chain = (u: number) => label(set, set.units.label[u]);
    const out: [string, string, boolean][] = [];
    for (let k = 0; k < v.placements.count; k++) {
      const b = v.placements.block[k];
      if (v.blocks.section[b] !== 1) continue;
      const a = set.annotations[v.placements.annotation[k]];
      out.push([chain(v.blocks.unit[b]), chain(a?.slots.unit[v.placements.slot[k]] ?? -1), v.placements.path[k] >= 0]);
    }
    return out;
  }

  it("case 9: each chain's helices only on its own block, though the alignment's rows reach the other's", async () => {
    expect(await chains(case9())).toEqual([
      ['A', 'A', false],
      ['B', 'B', false],
    ]);
  });

  it('data whose own chain the section does not show is still placed through the alignment', async () => {
    const hba = indexDoc(case9()).order.find((n) => n.kind === 'block' && n.params.unit === '1a3n/A')?.ref as string;
    expect(await chains(patched(case9(), [{ op: 'remove', ref: hba }]))).toEqual([
      ['B', 'A', true],
      ['B', 'B', false],
    ]);
  });

  it('a filter view or a numbering switch keeps the axis key, and every layer equal', async () => {
    const doc = case7();
    const { set } = await compile(doc, reg, fixtureLoaders());
    const layers = (d: SeqViewSpecDoc) => {
      const index = indexDoc(normalize(d, reg));
      const s = resolveSection(index, set, 'complex-one.main', reg).section;
      return { s, layers: s.layers.map((x) => resolveLayer(s, x.ref, index, set, reg, all).layer) };
    };
    const before = layers(doc);
    for (const params of [{ positions: 'iface' }, { numbering: 'label' }]) {
      const after = layers(patched(doc, [{ op: 'set', ref: 'complex-one.main', params }]));
      expect(after.s.axisKey).toBe(before.s.axisKey);
      expect(after.layers).toEqual(before.layers);
    }
    const filtered = layers(patched(doc, [{ op: 'set', ref: 'complex-one.main', params: { positions: 'iface' } }]));
    expect(filtered.s.blocks.map((b) => b.positions && OrderedSet.size(b.positions) > 0)).toEqual([true, true]);
  });
});

describe('fallbacks and resolve-phase issues', () => {
  it("VL-28: a block past the Int32 positions the layout holds is left out with an issue; a section's blocks together may hold more", async () => {
    const b = createSeqViewSpec({ title: 'big' });
    const aln = b.data.alignment({ ref: 'a', alphabet: 'protein', value: ['MKV-L', 'MKVAL'] });
    const [s1, s2] = ['s1', 's2'].map((ref) => b.data.sequence({ ref, value: 'MKVLAAG' }));
    b.view({ ref: 'one', label: 'One' }).block({ alignment: aln });
    const two = b.view({ ref: 'two', label: 'Two' });
    two.block({ sequence: s1 });
    two.block({ sequence: s2 });
    const doc = b.getDocument() as SeqViewSpecDoc;
    const { set, issues } = await compile(doc, reg, fixtureLoaders(files()));
    expect(issues).toEqual([]);
    const rangeIssues = (r: { issues: Issue[] }) => r.issues.filter((i) => i.code === 'svs.position.range');
    const column = set.alignments[0] as { columnCount: number } | null;

    (column as { columnCount: number }).columnCount = 3_000_000_000;
    const huge = resolveView(doc, set, 'one', reg, all);
    expect(huge.view.sections.ref.length).toBe(1);
    expect(rangeIssues(huge).map((i) => i.message)).toEqual([
      'the block holds 3,000,000,000 positions, past the 2,147,483,647 a block lays out: left out; show a `range` of it',
    ]);
    (column as { columnCount: number }).columnCount = 2_147_483_647;
    expect(rangeIssues(resolveView(doc, set, 'one', reg, all))).toEqual([]);

    // A section's blocks together may hold more (a whole genome): only a block is held to the Int32.
    set.units.rangeStop.fill(1_500_000_000);
    expect(rangeIssues(resolveView(doc, set, 'two', reg, all))).toEqual([]);
  });

  it('types, themes and glyphs fall back against what the viewer supports', async () => {
    const supported = new Set([...all].filter((x) => x !== 'rows' && x !== 'char'));
    const r = await resolve(case3(), 'family', supported);
    const l = Array.from(r.view.layers.ref).indexOf(r.set.strings.find('family.rows'));
    expect(label(r.set, r.view.layers.type[l])).toBe('letters');
    const flvr = Array.from(r.view.layers.ref).indexOf(r.set.strings.find('family.points'));
    expect(r.view.layers.themes[flvr][1]?.params?.glyph).toBe('circle'); // the glyph channel
    expect(codes(r.issues)).toEqual(['svs.render.fallback', 'svs.render.glyph']);
    expect(r.issues.filter((i) => i.code === 'svs.render.glyph').length).toBe(1); // one uniform char glyph, once
    const doc = patched(CASES.case4(), [{ op: 'set', ref: 'rna.rows', params: { color: { theme: 'clustal_x' } } }]);
    const rna = await resolve(doc, 'rna');
    expect(codes(rna.issues)).toEqual(['svs.render.theme']);
  });

  it("a fallback chain is checked against the section's layout as well as its axis", async () => {
    const make = (layout: 'tracks' | 'text') => {
      const b = createSeqViewSpec({ title: 'layouts' });
      b.data.sequence({ ref: 's', value: 'MKVLAAGIVG' });
      const d = b.data.annotation({
        ref: 'd',
        geometry: 'values',
        target: 's',
        value: { index: [0, 1], value: [1, 2] },
      });
      const p = b.data.annotation({ ref: 'p', geometry: 'pairs', target: 's', value: { start: [0], b_start: [5] } });
      const v = b.view({ ref: 'v', label: 'V', layout });
      v.track({ label: 'R' }).ruler();
      v.track({ label: 'D' }).values({ data: d });
      v.track({ label: 'P' }).pairs({ data: p });
      return b.getDocument();
    };
    const types = (r: { view: ResolvedView; set: SeqSet }) =>
      Array.from(r.view.layers.type, (x) => (x < 0 ? 'omitted' : label(r.set, x)));
    const tracks = await resolve(make('tracks'), 'v');
    expect([types(tracks), tracks.issues]).toEqual([['ruler', 'values', 'pairs'], []]);
    // Text draws values (a strip lane) but no ruler (omitted) or pairs (its fallback, points).
    const text = await resolve(make('text'), 'v');
    expect(types(text)).toEqual(['omitted', 'values', 'points']);
    expect(text.issues.map((i) => i.code)).toEqual(['svs.render.omitted', 'svs.render.fallback']);
    expect(text.view.sections.ref.length).toBe(1);
  });

  it("an alignment's rows in text are drawn as rows; a frame's rows are its span's rows track's lanes", async () => {
    const make = (layout: 'tracks' | 'text', from = 'msa') => {
      const b = createSeqViewSpec({ title: 'rows in text' });
      const aln = b.data.alignment({ ref: 'a', value: ['MK-LV', 'MKALV', 'M--LI'] } as never);
      const f = b.data.annotation({
        ref: 'f',
        geometry: 'ranges',
        value: [
          { target: 'a', start: 0, stop: 2, id: 'f1', label: 'All' },
          { target: 'a/1', start: 1, stop: 3, id: 'f2', label: 'Rows' },
          { target: 'a/2', start: 0, stop: 1, id: 'f2', label: 'Rows' },
        ],
      } as never);
      const view = b.view({ ref: 'v', label: 'V', layout });
      const v = view.section({ ref: 's1' } as never).block({ alignment: aln } as never);
      v.track({ ref: 'frames', label: 'F' }).frame({ ref: 'fl', data: f, span: { from } } as never);
      v.track({ ref: 'msa', label: 'M' }).rows();
      view
        .section({ ref: 's2' } as never)
        .block({ alignment: aln } as never)
        .track({ ref: 'other' })
        .rows();
      return b.getDocument();
    };
    const text = await resolve(make('text'), 'v');
    expect(Array.from(text.view.layers.type, (x) => label(text.set, x))).toEqual(['frame', 'rows', 'rows']);
    expect(text.issues).toEqual([]);
    // Column data on no lane; a row's on its lane in the rows track (lanes 0..2).
    expect(
      placed(text.view, text.set, 'fl')
        .map((p) => p[2])
        .sort(),
    ).toEqual([-1, 1, 2]);
    // A span naming a track the section does not hold: said, and drawn on its own band.
    const off = await resolve(make('tracks', 'other'), 'v');
    expect(off.issues.map((i) => [i.code, i.path])).toEqual([['svs.layer.span', expect.stringMatching(/span\/from$/)]]);
  });

  it("a matrix's default palette follows its domain: sequential when non-negative, diverging at 0 when signed", async () => {
    const make = (value: (number | null)[][], color?: { theme: string; params: Record<string, unknown> }) => {
      const b = createSeqViewSpec({ title: 'matrix' });
      b.data.sequence({ ref: 's', value: 'MKVLAAGIVG' });
      const m = b.data.annotation({
        ref: 'm',
        geometry: 'matrix',
        target: 's',
        states: ['A', 'B'],
        value: { index: [0, 1], value },
      });
      b.view({ ref: 'v', label: 'V' }).track({ label: 'M' }).matrix({ data: m, color });
      return b.getDocument();
    };
    const colorOf = async (value: (number | null)[][]) => {
      const r = await resolve(make(value), 'v');
      expect(r.issues).toEqual([]);
      return r.view.layers.themes[0][0];
    };
    const palette = (colors: string, value_domain: number[]) => ({
      theme: 'field_palette',
      params: { palette: { kind: 'continuous', colors, value_domain } },
    });
    // Non-negative data (counts with missing cells): Blues from 0, light low.
    expect(
      await colorOf([
        [1, null],
        [3, 2],
      ]),
    ).toMatchObject(palette('Blues', [0, 3])); // normalized: every param
    // Signed data: RdBu centered at 0 over the largest magnitude.
    expect(
      await colorOf([
        [-0.5, 2],
        [1, 0],
      ]),
    ).toMatchObject(palette('RdBu', [-2, 2]));
    // A color given in the document is kept.
    const given = make(
      [
        [1, 2],
        [0, 1],
      ],
      { theme: 'uniform', params: { color: 'red' } },
    );
    expect((await resolve(given, 'v')).view.layers.themes[0][0]).toMatchObject({
      theme: 'uniform',
      params: { color: 'red' },
    });
  });

  it('a range inside an A3M insert lands after it (svs.layer.insert); a focus outside its target', async () => {
    const b = createSeqViewSpec();
    const aln = b.data.alignment({ ref: 'a', value: ['ACGD', 'ACcGD'], flavor: 'a3m' });
    const ins = b.data.annotation({
      ref: 'ins',
      geometry: 'ranges',
      target: aln.row(1),
      value: { start: [2], stop: [3] },
    });
    const v = b.view({ ref: 'v', label: 'V' }).block({ alignment: aln });
    v.track({ label: 'Inserts' }).ranges({ data: ins });
    v.focus({ target: aln, start: 2, stop: 9 });
    const r = await resolve(b.getDocument(), 'v');
    expect(codes(r.issues)).toEqual(['svs.focus.range', 'svs.layer.insert']);
  });

  it('a focus by labels written backwards says so; one naming a label not there names it (CF-01)', async () => {
    const focus = async (beg: string, end: string) => {
      const b = createSeqViewSpec();
      const h = b.data
        .sequence({ ref: 'h', value: 'WVARIYPTNGYT' })
        .numbering({ name: 'auth', method: 'labels', labels: [50, 51, 52, '52A', '52B', 53, 54, 55, 56, 57, 58, 59] });
      b.view({ ref: 'v', label: 'V' })
        .block({ sequence: h })
        .focus({ target: h, numbering: 'auth', beg, end } as never);
      const r = await resolve(b.getDocument(), 'v');
      return r.issues.map((i) => `${i.code} ${i.message}`);
    };
    expect(await focus('53', '52A')).toEqual([
      expect.stringMatching(/^svs\.focus\.range .*beg "53" comes after end "52A" in numbering "auth": write beg ≤ end/),
    ]);
    expect(await focus('70', '75')).toEqual([expect.stringMatching(/label 70 is not in numbering "auth"/)]);
    expect(await focus('52', '53')).toEqual([]);
  });

  it('a ligand where a sequence is named, at load and at resolve (svs.ref.wrong_kind)', async () => {
    const doc = case6();
    const data = doc.root.children?.find((c) => c.kind === 'data');
    const at = data?.children?.push({ kind: 'unit', ref: 'wrong', params: { ligand: '1iep/entity:1' } }) ?? 0;
    const main = doc.root.children?.find((c) => c.ref === 'chains')?.children?.[0];
    main?.children?.unshift({ kind: 'block', ref: 'sti-seq', params: { sequence: '1iep/entity:2' } });
    const { set, issues } = await compile(doc, reg, fixtureLoaders());
    expect(issues.filter((i) => i.code === 'svs.ref.wrong_kind').map((i) => i.path)).toEqual([
      `/root/children/0/children/${at - 1}/params/ligand`,
    ]);
    expect(set.byName('wrong', t)).toBe(false);
    const r = resolveView(doc, set, 'chains', reg, all);
    expect(r.issues.map((i) => [i.code, i.path])).toEqual([
      ['svs.ref.wrong_kind', '/root/children/1/children/0/children/0/params/sequence'],
    ]);
    expect(r.view.blocks.count).toBe(2);
    expect(KeyKind.unit).toBe(1);
  });
});

describe('split rows', () => {
  const trackOf = (v: ResolvedView, set: SeqSet, ref: string) =>
    Array.from(v.tracks.ref).indexOf(set.strings.find(ref));
  it("a track's split layer and its state at open: variants expanded, ranges collapsed, `expanded` wins", async () => {
    const r = await resolve(CASES.case13(), 'features');
    const nv = trackOf(r.view, r.set, 'features.natural-variants');
    expect(label(r.set, r.view.layers.ref[r.view.tracks.split[nv]])).toBe('features.variants-2'); // the group is "features.variants"
    expect(r.view.tracks.expanded[nv]).toBe(1);
    expect(r.view.tracks.split[trackOf(r.view, r.set, 'features.chain')]).toBe(-1);
    const doc = case7();
    const one = await resolve(doc, 'complex-one');
    const uni = trackOf(one.view, one.set, 'complex-one.uniprot');
    expect([label(one.set, one.view.layers.ref[one.view.tracks.split[uni]]), one.view.tracks.expanded[uni]]).toEqual([
      'complex-one.ranges',
      0,
    ]);
    const open = await resolve(
      patched(doc, [{ op: 'set', ref: 'complex-one.uniprot', params: { expanded: true } }]),
      'complex-one',
    );
    expect(open.view.tracks.expanded[uni]).toBe(1);
  });

  it("splitRows: order's values first, every one, then the others in data order; missing as -; the filter's rows -1", async () => {
    const b = createSeqViewSpec({ title: 'split' });
    const s = b.data.sequence({ ref: 's', value: 'MKVLAAGIV' });
    const value = { start: [0, 1, 2, 3, 4], alt: ['W', 'R', '', 'W', '*'], category: ['a', 'a', 'a', 'b', 'a'] };
    const a = b.data.annotation({ ref: 'a', geometry: 'ranges', target: s, value });
    b.view({ ref: 'v', label: 'V' })
      .track({ ref: 't' })
      .variants({ ref: 'l', data: a, split: 'alt', filter: { field: 'category', values: ['a'] } });
    const doc = b.getDocument();
    expect(doc.root.children?.[1].children?.[0].children?.[0].children?.[0].params?.split).toEqual({
      by: 'alt',
      order: [...'GAVLISTCMDNEQRKHFYWP', '-', '*'],
    });
    const r = await resolve(doc, 'v');
    const ann = r.set.annotations[r.view.data[0]];
    const split = { by: 'alt', order: ['*', 'Q'] };
    const { values, rowOf } = splitRows(r.set, ann ? [ann] : [], [r.view.kept[0]], split);
    expect(values).toEqual(['*', 'Q', 'W', 'R', '-']);
    expect(Array.from(rowOf[0])).toEqual([2, 3, 4, -1, 0]);
    // `empty: "hide"`: only the values a kept row has ("Q" none; "W"'s second row filtered out), in the same order.
    const hidden = splitRows(r.set, ann ? [ann] : [], [r.view.kept[0]], { ...split, empty: 'hide' });
    expect(hidden.values).toEqual(['*', 'W', 'R', '-']);
    expect(Array.from(hidden.rowOf[0])).toEqual([1, 2, 3, -1, 0]);
  });

  it("split's `empty` defaults to show and is kept as written", async () => {
    const b = createSeqViewSpec({ title: 'split' });
    const s = b.data.sequence({ ref: 's', value: 'MKVLAAGIV' });
    const a = b.data.annotation({ ref: 'a', geometry: 'ranges', target: s, value: { start: [0], alt: ['W'] } });
    const t = b.view({ ref: 'v', label: 'V' }).track({ ref: 't' });
    t.variants({ ref: 'l', data: a, split: 'alt' });
    t.variants({ ref: 'h', data: a, split: { by: 'alt', order: ['A', 'W'], empty: 'hide' } });
    const view = normalize(b.getDocument(), reg).root.children?.[1];
    const layers = view?.children?.[0]?.children?.[0]?.children ?? [];
    expect(layers.map((l) => (l.params?.split as { empty?: string } | undefined)?.empty)).toEqual(['show', 'hide']);
  });

  it('a split by a field no row of the data has is svs.layer.split when resolved', async () => {
    const doc = case7();
    const at = 'complex-one.ranges';
    const r = await resolve(patched(doc, [{ op: 'set', ref: at, params: { split: { by: 'alt' } } }]), 'complex-one');
    expect(r.issues.filter((i) => i.code === 'svs.layer.split').map((i) => i.path)).toEqual([
      expect.stringMatching(/\/params\/split$/),
    ]);
  });
});

/** A layer's index by ref. */
const layerOf = (v: ResolvedView, set: SeqSet, ref: string) => Array.from(v.layers.ref).indexOf(set.strings.find(ref));
/** One sequence, its data, a view built by `view`. */
function small(build: (b: ReturnType<typeof createSeqViewSpec>, s: string) => void) {
  const b = createSeqViewSpec();
  b.data.sequence({ ref: 's', value: 'MKVLAAGIVGLLLAQ' });
  build(b, 's');
  return b.getDocument();
}

describe('VW4: channels, named rule tables and one palette per layer', () => {
  it('a theme per channel of the type: color first, a glyph on ranges, points and pairs, then the size channels', async () => {
    const r = await resolve(case3(), 'family');
    for (let l = 0; l < r.view.layers.count; l++) {
      const type = label(r.set, r.view.layers.type[l]);
      const sizes: Record<string, number> = { ranges: 2, points: 1, pairs: 3, values: 1 };
      const glyph = ['ranges', 'points', 'pairs'].includes(type) ? 1 : 0;
      expect(r.view.layers.themes[l].length, type).toBe(1 + glyph + (sizes[type] ?? 0));
    }
    const values = layerOf(r.view, r.set, 'family.values');
    expect(r.view.layers.themes[values]).toEqual([null, null]); // values: the viewer's ink and its bar's thickness
  });

  it('named tables reach the viewer as rules, nested themes resolved alike', async () => {
    const r = await resolve(case10(), 'alignment');
    const rows = r.view.layers.themes[Array.from(r.view.layers.type, (x) => label(r.set, x)).indexOf('rows')][0];
    expect(rows?.theme).toBe('rules');
    expect(rows?.params?.rules).toMatchObject([
      { when: [{ measure: { equals_row: 0 }, min: 1 }], theme: { theme: 'letter' } },
    ]);
    // identity_to_row over clustal_x: both expanded; clustal_x's rules colored by share.
    const nested = patched(case10(), [
      {
        op: 'set',
        ref: 'alignment.rows',
        params: { color: { theme: 'identity_to_row', params: { base: { theme: 'clustal_x' } } } },
      },
    ]);
    const n = await resolve(nested, 'alignment');
    const t = n.view.layers.themes[layerOf(n.view, n.set, 'alignment.rows')][0] as ThemeRef;
    const base = ((t.params?.rules ?? []) as { theme: ThemeRef }[])[0].theme;
    expect([base.theme, ((base.params?.rules ?? []) as unknown[]).length > 30]).toEqual(['rules', true]);
  });

  it('clustal_x on a lazy alignment is letter, reported', async () => {
    const rows = indexDoc(case9()).order.find((n) => n.kind === 'layer' && n.params.type === 'rows')?.ref as string;
    const d = patched(case9(), [{ op: 'set', ref: rows, params: { color: { theme: 'clustal_x' } } }]);
    const { set } = await compile(d, reg, fixtureLoaders());
    const r = resolveView(d, set, 'msa', reg, all);
    expect(r.view.layers.themes[layerOf(r.view, set, rows)][0]).toMatchObject({ theme: 'letter' });
    expect(codes(r.issues)).toEqual(['svs.render.theme']);
  });

  it("the coverage provider's default palette; numeric defaults for color over values", async () => {
    const doc = small((b, s) => {
      const u = b.data.unit({ ref: 'u', sequence: s, observed: [[0, 5]] } as never);
      const c = b.data.computed({ ref: 'cov', provider: 'coverage', of: u } as never);
      const d = b.data.annotation({ ref: 'd', geometry: 'values', target: s, value: { index: [0, 1], value: [2, 4] } });
      const v = b.view({ ref: 'v', label: 'V' }).block({ unit: u });
      v.track({ label: 'C' }).ranges({ data: c });
      v.track({ label: 'D' }).color({ data: d });
    });
    const r = await resolve(doc, 'v');
    expect(r.view.layers.themes[0][0]?.params?.palette).toMatchObject({ colors: { unobserved: '#bab0ac' } });
    expect(r.view.layers.themes[1][0]?.params?.palette).toMatchObject({
      kind: 'continuous',
      colors: 'Blues',
      value_domain: [0, 4],
    });
  });

  it("a filtered layer's domain and colors over every row", async () => {
    const doc = small((b, s) => {
      const value = { start: [0, 1, 2], value: [1, 5, 9], category: ['SH3', 'SH2', 'Kinase'] };
      const a = b.data.annotation({ ref: 'a', geometry: 'ranges', target: s, value });
      const t = b.view({ ref: 'v', label: 'V' }).track({ label: 'T' });
      t.color({ ref: 'all', data: a });
      t.color({ ref: 'some', data: a, filter: { field: 'category', values: ['SH2'] } as never });
    });
    const r = await resolve(doc, 'v');
    const [all, some] = ['all', 'some'].map((x) => layerOf(r.view, r.set, x));
    expect([r.view.layers.domainMin[some], r.view.layers.domainMax[some]]).toEqual([1, 9]);
    expect(Array.from(r.view.kept[r.view.layers.dataStart[some]] ?? [])).toEqual([0, 1, 0]);
    const palette = (l: number) => {
      const data = Array.from(
        r.view.data.subarray(r.view.layers.dataStart[l], r.view.layers.dataStart[l] + 1),
        (id) => r.set.annotations[id] as never,
      );
      const p = r.view.layers.themes[l][0]?.params?.palette as Palette; // normalized: Category10
      return layerPalette(r.set, data, p, null, [r.view.layers.domainMin[l], r.view.layers.domainMax[l]]);
    };
    expect([palette(all)('SH2'), palette(some)('SH2')]).toEqual([ColorLists.Category10[1], ColorLists.Category10[1]]);
  });

  it("a matrix's domain is a param; data with other states are left out", async () => {
    const doc = small((b, s) => {
      const value = {
        index: [0, 1],
        value: [
          [1, 2],
          [3, 4],
        ],
      };
      const m1 = b.data.annotation({ ref: 'm1', geometry: 'matrix', target: s, states: ['x', 'y'], value });
      const m2 = b.data.annotation({ ref: 'm2', geometry: 'matrix', target: s, states: ['y', 'z'], value });
      b.view({ ref: 'v', label: 'V' })
        .track({ label: 'M' })
        .matrix({ data: [m1, m2], domain: [0, 10] } as never);
    });
    expect(validate(doc, reg)).toEqual([]);
    const r = await resolve(doc, 'v');
    // valid data the layer leaves out: a loss, a warning, not invalid input (GV-10)
    expect(r.issues.map((i) => `${i.severity} ${i.code}`)).toEqual(['warning svs.render.omitted']);
    expect([r.view.layers.dataCount[0], r.view.layers.placementCount[0]]).toEqual([1, 1]);
    expect(r.view.layers.themes[0][0]?.params?.palette).toMatchObject({ value_domain: [0, 10] });
  });
});

describe('VW4: placement and sections', () => {
  it('the first rows layer of a section wins; a later one is left out', async () => {
    const doc = small((b) => {
      const aln = b.data.alignment({ ref: 'a', value: ['MKVL', 'MKIL', 'MRVL'] });
      const v = b.view({ ref: 'v', label: 'V' }).block({ alignment: aln });
      v.track({ ref: 't1' }).rows({ row_subset: [0, 1] });
      v.track({ ref: 't2' }).rows({ row_subset: [2] });
    });
    const r = await resolve(doc, 'v');
    expect(r.issues.map((i) => [i.code, i.severity])).toEqual([['svs.layer.rows', 'warning']]);
    expect([r.view.layers.count, Array.from(r.view.lanes.row), Array.from(r.view.tracks.laneCount)]).toEqual([
      1,
      [0, 1],
      [2, 0],
    ]);
  });

  it("a track's split layer by ref when an earlier layer is left out", async () => {
    const doc = case13();
    const { set } = await compile(doc, reg, fixtureLoaders());
    const index = indexDoc(normalize(doc, reg));
    const s = resolveSection(index, set, 'features.main', reg).section;
    const layers = s.layers.map((x) => resolveLayer(s, x.ref, index, set, reg, all).layer);
    const v = composeView('features', [s], [layers.slice(1)], set);
    const nv = Array.from(v.tracks.ref).indexOf(set.strings.find('features.natural-variants'));
    expect(label(set, v.layers.ref[v.tracks.split[nv]])).toBe('features.variants-2');
  });

  it('a ranges row with no value counts in positions and a focus', async () => {
    const doc = small((b, s) => {
      const value = { start: [1, 6], stop: [3, 9], value: [0.5, null], id: ['a', 'b'] };
      const d = b.data.annotation({ ref: 'd', geometry: 'ranges', target: s, value });
      b.view({ ref: 'v', label: 'V', positions: d })
        .block({ sequence: s })
        .focus({ annotation: d, items: ['b'] })
        .track({ label: 'L' })
        .letters();
    });
    const r = await resolve(doc, 'v');
    const shown = r.view.blocks.positions[0];
    expect(shown && OrderedSet.size(shown)).toBe(5);
    expect(r.view.sections.focus[0]).toEqual({ block: 0, start: 6, stop: 9 });
  });

  it('`positions` entries: a ref, or `{ data, filter }` covering only the rows the filter keeps; a bad entry is reported', async () => {
    const make = (positions: unknown) =>
      small((b, s) => {
        const value = { start: [1, 6, 12], stop: [3, 9, 15], category: ['a', 'm', 'm'] };
        b.data.annotation({ ref: 'd', geometry: 'ranges', target: s, value });
        b.view({ ref: 'v', label: 'V', positions: positions as never })
          .block({ sequence: s })
          .track({ label: 'L' })
          .letters();
      });
    const shown = async (positions: unknown) => {
      const set = (await resolve(make(positions), 'v')).view.blocks.positions[0];
      return set && Array.from({ length: OrderedSet.size(set) }, (_, i) => OrderedSet.getAt(set, i));
    };
    expect(await shown('d')).toEqual([1, 2, 6, 7, 8, 12, 13, 14]);
    expect(await shown({ data: 'd', filter: { values: ['m'] } })).toEqual([6, 7, 8, 12, 13, 14]);
    expect(await shown([{ data: 'd', filter: { field: 'category', values: ['a'] } }, 'd'])).toHaveLength(8);
    expect(validate(make({ data: 'd' }), reg).length).toBeGreaterThan(0);
    expect(validate(make({ data: 'd', filter: { values: ['m'] } }), reg)).toEqual([]);
  });

  it('a focus on an unavailable target is quiet; on a derived ref its node did not yield, reported', async () => {
    const unknown = patched(case10(), [
      { op: 'insert', parent: 'strip.main', node: { kind: 'focus', ref: 'f', params: { target: '1iep/Z' } } },
    ]);
    const { set } = await compile(unknown, reg, fixtureLoaders());
    expect(codes(resolveView(unknown, set, 'strip', reg, all).issues)).toEqual(['svs.ref.derived_unknown']);
    // Without the structure extension, the complex is unavailable: s3's focus on "1iep/A" is ignored quietly.
    const focused = snapshotDoc(case10(), 's3', reg).doc;
    const none = registries(false);
    const quiet = await compile(focused, none, fixtureLoaders());
    expect(resolveView(focused, quiet.set, 'strip', none, supportedAll(none)).issues).toEqual([]);
  });

  it("letters of every protein mapped to an alignment's row land on the alignment block", async () => {
    const doc = small((b) => {
      const aln = b.data.alignment({ ref: 'a', value: ['ATGGCTTCC', 'ATGGCTTCC'], alphabet: 'dna' } as never);
      b.data.sequence({ ref: 'p', alphabet: 'protein', value: 'MAS' });
      b.data.mapping({ ref: 'm', from: 'p', to: 'a/0', method: 'codon', cds: [[0, 9]] } as never);
      b.view({ ref: 'v', label: 'V' }).block({ alignment: aln }).track({ label: 'P' }).letters({ of: 'mapped' });
    });
    const r = await resolve(doc, 'v');
    expect(r.issues).toEqual([]);
    expect(placed(r.view, r.set, r.set.strings.get(r.view.layers.ref[0]))).toEqual([[0, -1, -1, 0]]);
    expect(r.view.paths[0].span).toBe(3);
  });

  it("an empty table reaches nothing and is not unplaced (VW2's item)", async () => {
    const r = await resolve(
      small((b, s) => {
        const e = b.data.annotation({ ref: 'e', geometry: 'ranges', target: s, value: [] });
        b.view({ ref: 'v', label: 'V' }).track({ label: 'E' }).ranges({ data: e });
      }),
      'v',
    );
    expect(r.issues).toEqual([]);
  });
});

describe('VW4: links memoized per set (V7)', () => {
  it('a second compose searches no path', async () => {
    const doc = case10();
    const { set } = await compile(doc, reg, fixtureLoaders());
    const index = indexDoc(normalize(doc, reg));
    const sections = ['both.ann', 'both.msa'].map((ref) => resolveSection(index, set, ref, reg).section);
    const layers = sections.map((s) => s.layers.map((l) => resolveLayer(s, l.ref, index, set, reg, all).layer));
    const first = composeView('both', sections, layers, set);
    const before = pathsSearched();
    const again = composeView('both', sections, layers, set);
    expect([pathsSearched() - before, first.links.count, again.links.count]).toEqual([0, 1, 1]);
  });

  it('120 chains in two sections: copies are never searched, the same unit links once per pair', async () => {
    const b = createSeqViewSpec();
    const s = b.data.sequence({ ref: 's', value: 'MKVLAAGIVG' });
    const units = Array.from({ length: 120 }, (_, k) => b.data.unit({ ref: `u${k}`, sequence: s }));
    const v = b.view({ ref: 'v', label: 'V', copies: 'each' });
    for (const u of units) v.block({ unit: u });
    const w = v.section({ ref: 'v.two', copies: 'each' });
    for (const u of units) w.block({ unit: u });
    v.track({ label: 'L' }).letters();
    const doc = b.getDocument();
    const { set } = await compile(doc, reg, fixtureLoaders());
    const t0 = pathsSearched();
    const first = resolveView(doc, set, 'v', reg, all).view;
    const t1 = pathsSearched();
    resolveView(doc, set, 'v', reg, all);
    expect([first.links.count, t1 - t0, pathsSearched() - t1]).toEqual([120, 120, 0]);
  });

  // 3J3Q: 1,356 copies of one chain merged into one block, each with its own secondary structure, searched every copy
  // from every copy (1.8M dead-end searches, seconds). A unit reaches only itself among its copies, unsearched.
  const copies = (n: number, mapped: boolean) => {
    const b = createSeqViewSpec();
    const s = b.data.sequence({ ref: 's', value: 'MKVLAAGIVG' });
    const units = Array.from({ length: n }, (_, k) => b.data.unit({ ref: `u${k}`, sequence: s }));
    if (mapped) b.data.mapping({ from: units[0], to: units[1], method: 'identity' });
    const per = b.data.annotation({
      ref: 'per',
      geometry: 'ranges',
      value: { target: units, start: units.map((_, k) => k % 5), stop: units.map(() => 8) },
    });
    const v = b.view({ ref: 'v', label: 'V' });
    for (const u of units) v.block({ unit: u });
    v.track({ label: 'L' }).ranges({ ref: 'v.per', data: per });
    return b.getDocument();
  };

  it('1,000 merged copies with data of their own: each on its own member, no path searched', async () => {
    const doc = copies(1000, false);
    const { set } = await compile(doc, reg, fixtureLoaders());
    const t0 = pathsSearched();
    const { view } = resolveView(doc, set, 'v', reg, all);
    const ps = placed(view, set, 'v.per');
    expect([view.blocks.count, ps.length, pathsSearched() - t0]).toEqual([1, 1000, 0]);
    expect(ps.every((p, k) => p[1] === id(set, `u${k}`) && p[3] === -1)).toBe(true);
  });

  it('a mapping onto a unit still takes one copy onto another', async () => {
    const doc = copies(3, true);
    const { set } = await compile(doc, reg, fixtureLoaders());
    const ps = placed(resolveView(doc, set, 'v', reg, all).view, set, 'v.per');
    // Each copy's own, unpathed; and through the mapping, u0's data on u1 and u1's on u0.
    expect(ps.filter((p) => p[3] === -1).map((p) => p[1])).toEqual([0, 1, 2].map((k) => id(set, `u${k}`)));
    expect(ps.filter((p) => p[3] >= 0).map((p) => p[1])).toEqual([id(set, 'u1'), id(set, 'u0')]);
  });
});
