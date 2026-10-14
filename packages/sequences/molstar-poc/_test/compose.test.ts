/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Structures as data, and 3D views composed from the viewer (seqspec:libraries#the-mol-adapter): the PDB coverage
// provider's rows and fields, Seq* runs onto a structure's residues as MVS color tables, and MVS states composed from a
// headless viewer's state (a track's colors as drawn, the structure an item names, its default coloring and a switch of
// it, the selection), checked as data: no canvas.

import { readFileSync } from 'node:fs';
import { compile } from '@molstar/seq-core/format';
import { Extents, KeyKind, OrderedSet, type SeqSet } from '@molstar/seq-core/model';
import {
  ChannelMarks,
  createSeqViewer,
  DefaultRegistry,
  DefaultSeqViewerSpec,
  manualScheduler,
} from '@molstar/seq-viewer';
import { MVSData } from '@molstar/mvs-builder/mvs-data';
import { describe, expect, it } from 'vitest';
import type { SplitFeatures } from '../../viewer/src/layout.js';
import { fixtureMeasurer } from '../../viewer/src/testing/index.js';
import {
  composeStructureMvs,
  coverageRows,
  MolstarSvsEntry,
  type MvsColorTable,
  readStructure,
  StructureLink,
  type StructureRef,
  structureOf,
} from '../src/index.js';
import {
  colorByItemAction,
  colorTracks,
  composeMvs,
  composeStructureView,
  coordinatesOf,
  itemStructure,
  ownTracks,
  type Shown,
  showInStructureAction,
  splitRow,
  structureTrackActions,
} from '../src/viewer.js';
import { bytesOf, demoLoaders, P53, registries } from './fixtures.js';
import {
  coverageTrack,
  P53_STRUCTURES,
  P53_TRACKS,
  PTM_COLORS,
  p53CoverageOptions,
  p53DemoDoc,
  p53ShownPatch,
} from './p53-demo.js';

const unit = (set: SeqSet, ref: string) => {
  const t = { kind: 0, id: 0, row: -1 };
  if (!set.byName(ref, t) || t.kind !== KeyKind.unit) throw new Error(ref);
  return { kind: 1, id: t.id, row: -1 };
};
const tableOf = (uri: string) => JSON.parse(atob(uri.slice(uri.indexOf(',') + 1))) as MvsColorTable;
type Node = { kind: string; params?: Record<string, unknown>; children?: Node[] };
const nodes = (n: Node, kind: string): Node[] => [
  ...(n.kind === kind ? [n] : []),
  ...(n.children ?? []).flatMap((c) => nodes(c, kind)),
];
const colorsOf = (t: MvsColorTable) => [...new Set(t.color)].sort();
const ref = (label: string): StructureRef => {
  const s = P53_STRUCTURES.find((x) => x.label === label);
  if (!s) throw new Error(label);
  return { source: s.source, ...('assembly' in s ? { assembly: s.assembly } : {}) };
};

describe('molstar:sifts_coverage', () => {
  it('a row per resolved stretch of each chain, in UniProt numbering, each naming its structure', async () => {
    const { set, issues } = await compile(p53DemoDoc(P53.sequence), registries(), demoLoaders());
    expect(issues.filter((i) => i.severity === 'error')).toEqual([]);
    const a = set.annotations.find((x) => x && set.strings.get(x.key) === 'coverage');
    if (!a) throw new Error('no coverage');
    // The demo's rows sit on the chains' units (instance space); their UniProt numbers as the provider reads them.
    expect(set.strings.get(set.units.key[a.target.unit[0]])).toBe('1tup/C');
    const c = coverageRows(set, p53CoverageOptions);
    const rows = c.label.map((label, r) => [label, c.beg[r], c.end[r], c.structure[r]]);
    // 1TUP's chains A, B, C (label_asym_id C, D, E) resolve 94–289, 96–289 and 95–289 of P04637; the model all 393.
    expect(rows.slice(0, 3).map((r) => r.slice(0, 3))).toEqual([
      ['1TUP chain A (label_asym_id C)', 94, 289],
      ['1TUP chain B (label_asym_id D)', 96, 289],
      ['1TUP chain C (label_asym_id E)', 95, 289],
    ]);
    expect(rows[0][3]).toEqual({ source: { service: 'pdbe', id: '1tup' }, chains: ['C'], authChains: ['A'] });
    expect(rows.find((r) => r[0] === '1C26 chain A')?.[3]).toEqual({
      source: { service: 'pdbe', id: '1c26' },
      chains: ['A'],
      assembly: '1',
    });
    expect(rows.at(-1)).toEqual([
      'AlphaFold chain A',
      1,
      393,
      { source: { service: 'afdb', id: 'AF-P04637-F1' }, chains: ['A'] },
    ]);
    expect(a.custom.doc).toEqual({ molstar: { representation: 'cartoon', base: '#c8c8c8' } });
  });
});

describe('composeStructureMvs', () => {
  it('Seq* runs on P04637 become residue_range tables on every chain they reach, later tables winning', async () => {
    const { set } = await compile(p53DemoDoc(P53.sequence), registries(), demoLoaders());
    const p53 = unit(set, 'p53/unit');
    const runs = Extents.create(2);
    Extents.add(runs, p53, 93, 100);
    Extents.add(runs, p53, 174, 175);
    const r = composeStructureMvs({
      set,
      structure: { source: { service: 'pdbe', id: '1tup' }, chains: ['C'] },
      layers: [{ runs, colors: ['#ff0000', '#0000ff'] }],
    });
    if (!r) throw new Error('no state');
    expect(r.tables).toEqual([
      {
        label_asym_id: ['C', 'D', 'E', 'C', 'D', 'E'],
        beg_label_seq_id: [1, 1, 1, 82, 82, 82],
        end_label_seq_id: [7, 7, 7, 82, 82, 82],
        color: ['#ff0000', '#ff0000', '#ff0000', '#0000ff', '#0000ff', '#0000ff'],
      },
    ]);
    expect(MVSData.validationIssues(r.state)).toBeUndefined();
    const colorFromUri = nodes(r.state.root as Node, 'color_from_uri');
    expect(colorFromUri.map((n) => [n.params?.schema, n.params?.format])).toEqual([['residue_range', 'json']]);
    expect(tableOf(colorFromUri[0].params?.uri as string)).toEqual(r.tables[0]);
    expect(nodes(r.state.root as Node, 'download')[0].params?.url).toBe(
      'https://www.ebi.ac.uk/pdbe/entry-files/download/1tup.bcif',
    );
  });
});

async function viewer() {
  const spec = DefaultSeqViewerSpec();
  const ctx = createSeqViewer(
    { ...spec, registry: [DefaultRegistry, { svs: MolstarSvsEntry }] },
    { scheduler: manualScheduler(), loaders: demoLoaders(), text: fixtureMeasurer() },
  );
  const r = await ctx.setDocument(p53DemoDoc(P53.sequence));
  expect(r.issues.filter((i) => i.severity === 'error')).toEqual([]);
  return ctx;
}
type Ctx = Awaited<ReturnType<typeof viewer>>;

describe("the coverage's chain rows in instance space (QM1b, HR1)", () => {
  it("1TUP's rows on its chains' units; a chain's residue lights its own row alone, P04637's every row", async () => {
    const ctx = await viewer();
    const a = ctx.graph.assembled.value;
    const set = a?.set as SeqSet;
    const view = a?.view;
    const k = view?.resolved.tracks.split[trackOf(ctx, coverageTrack(P53_STRUCTURES[0]))] ?? -1;
    const split = (view?.prepared[k]?.state as SplitFeatures | undefined)?.split;
    // 1TUP's p53 chains by their author IDs (A, B, C; label_asym_id C, D, E), each on its unit.
    expect(split?.labels).toEqual(['A', 'B', 'C']);
    expect(Array.from(split?.unit ?? [])).toEqual(['1tup/C', '1tup/D', '1tup/E'].map((r) => unit(set, r).id));
    const ch = new ChannelMarks(ctx);
    const rows = (bit: 'highlight' | 'hover') =>
      ['1tup/C', '1tup/D', '1tup/E'].map((r) => {
        const m = ch.marks(0)?.ownOf(unit(set, r).id)?.[bit];
        return Array.from({ length: (m?.to(0) ?? 0) - (m?.from(0) ?? 0) }, (_, i) => [
          m?.start[(m?.from(0) ?? 0) + i],
          m?.stop[(m?.from(0) ?? 0) + i],
        ]);
      });
    const seq = (ref: string, i: number) =>
      ({
        kind: 'seq-loci',
        set,
        elements: [{ unit: unit(set, ref).id, indices: OrderedSet.ofSingleton(i) }],
      }) as never;
    ctx.highlight.set('t', seq('1tup/D', 81));
    expect(rows('highlight')).toEqual([[], [[174, 175]], []]);
    ctx.highlight.set('t', seq('p53/unit', 174));
    expect(rows('highlight')).toEqual([[[174, 175]], [[174, 175]], [[174, 175]]]);
    // Another structure's chains are not 1TUP's: a run on 1TUP's chain B lights none of 2OCJ's rows.
    ctx.highlight.set('t', seq('1tup/D', 81));
    expect(ch.marks(0)?.ownOf(unit(set, '2ocj/B').id)?.highlight.count).toBe(0);
    // A 3D hover on chain B's R175 (its chain system's coordinates, as linkViews sends them): that chain's row alone.
    const structure = await structureOf(await readStructure(bytesOf('1tup.bcif'), 'bcif'), 0, null);
    const link = StructureLink.attach(ctx.systems, '3d', structure, set, '1tup');
    const runs = Extents.create(1);
    Extents.add(runs, unit(set, '1tup/D'), 81, 82);
    const hover3d = coordinatesOf(link, link.fromLoci(link.toLoci(runs), Extents.create(1)));
    expect(hover3d.systems).toEqual(['molstar:3d/D@1_555']);
    ctx.highlight.set('t', hover3d);
    expect(rows('highlight')).toEqual([[], [[174, 175]], []]);
  });

  it('2OCJ shown: its entry marked, its four rows each on its chain', async () => {
    const ctx = await viewer();
    const r = await ctx.applyPatch(p53ShownPatch(P53_STRUCTURES[1]));
    expect(r.issues.filter((i) => i.severity === 'error')).toEqual([]);
    const view = ctx.graph.assembled.value?.view;
    const t = trackOf(ctx, coverageTrack(P53_STRUCTURES[1]));
    const split = (view?.prepared[view.resolved.tracks.split[t]]?.state as SplitFeatures | undefined)?.split;
    expect(split?.labels).toEqual(['A', 'B', 'C', 'D']);
    expect(view?.set.strings.get(view.resolved.tracks.label[t])).toBe('2OCJ · in 3D');
    // A residue of chain B lights its own row alone.
    const set = view?.set as SeqSet;
    const b = unit(set, '2ocj/B').id;
    ctx.highlight.set('t', {
      kind: 'seq-loci',
      set,
      elements: [{ unit: b, indices: OrderedSet.ofSingleton(81) }],
    } as never);
    const ch = new ChannelMarks(ctx);
    const units = Array.from(split?.unit ?? []);
    expect(units.map((u) => ch.marks(0)?.ownOf(u)?.highlight.count)).toEqual(units.map((u) => (u === b ? 1 : 0)));
  });
});
const trackOf = (ctx: Ctx, ref: string) => {
  const v = ctx.graph.assembled.value?.view;
  for (let t = 0; v && t < v.resolved.tracks.count; t++)
    if (v.set.strings.get(v.resolved.tracks.ref[t]) === ref) return t;
  throw new Error(ref);
};
/** The coverage item whose label starts with `label`. */
function itemOf(ctx: Ctx, label: string, split = -1) {
  const view = ctx.graph.assembled.value?.view;
  const set = view?.set as SeqSet;
  const a = set.annotations.find((x) => x && set.strings.get(x.key) === 'coverage');
  const r = Array.from({ length: a?.count ?? 0 }, (_, k) => k).find((k) =>
    set.strings.get(a?.rowLabel?.[k] ?? -1).startsWith(label),
  );
  return { annotation: 'coverage', item: a?.item[r ?? -1] ?? -1, track: P53_TRACKS.coverage, split };
}

describe('3D views composed from the viewer', () => {
  it('an item names its whole structure, or one chain when hovered on its own row', async () => {
    const ctx = await viewer();
    const view = ctx.graph.assembled.value?.view;
    if (!view) throw new Error('no view');
    expect(itemStructure(view, itemOf(ctx, '1TUP chain B'))?.structure).toEqual({
      source: { service: 'pdbe', id: '1tup' },
    });
    expect(itemStructure(view, itemOf(ctx, '1TUP chain B', 1))?.structure).toMatchObject({ chains: ['D'] });
    expect(itemStructure(view, itemOf(ctx, '1TUP chain B'))?.defaults).toEqual({
      representation: 'cartoon',
      base: '#c8c8c8',
    });
    ctx.dispose();
  });

  it("each structure's default coloring is its own data: the model's pLDDT; a PDB entry none", async () => {
    const ctx = await viewer();
    const view = ctx.graph.assembled.value?.view;
    if (!view) throw new Error('no view');
    expect(ownTracks(view, ref('AlphaFold'))).toEqual([P53_TRACKS.plddt]);
    expect(ownTracks(view, ref('1TUP'))).toEqual([]);
    const af = composeMvs(ctx, itemOf(ctx, 'AlphaFold'));
    expect(af?.tables.map(colorsOf)).toEqual([['#0053d6', '#65cbf3', '#ff7d45', '#ffdb13']]);
    expect(af && MVSData.validationIssues(af.state)).toBeUndefined();
    expect(composeMvs(ctx, itemOf(ctx, '1TUP chain A'))?.tables).toEqual([]);
    ctx.dispose();
  });

  it("the tracks that can color a structure: a sequence's through SIFTS on any, pLDDT only on its model, first", async () => {
    const ctx = await viewer();
    const af = colorTracks(ctx, ref('AlphaFold')).map((x) => x.track);
    const tup = colorTracks(ctx, ref('1TUP')).map((x) => x.track);
    expect(af[0]).toBe(P53_TRACKS.plddt);
    for (const t of [
      P53_TRACKS.domains,
      P53_TRACKS.regions,
      P53_TRACKS.letters,
      P53_TRACKS.ptms,
      P53_TRACKS.sites,
      P53_TRACKS.variants,
      P53_TRACKS.hydropathy,
    ]) {
      expect(af).toContain(t);
      expect(tup).toContain(t);
    }
    expect(tup).not.toContain(P53_TRACKS.plddt);
    expect(tup.some((t) => t.startsWith(P53_TRACKS.coverage))).toBe(false);
    ctx.dispose();
  });

  it('1TUP colored by its domains: the DNA-binding domain and the NLS through SIFTS, each in its Seq* color', async () => {
    const ctx = await viewer();
    const r = composeStructureView(ctx, ref('1TUP'), { tracks: [P53_TRACKS.domains] });
    if (!r) throw new Error('nothing composed');
    const t = r.tables[0];
    // DNA binding 102–292 is label_seq_id 9–199 of chains C to E; the NLS 305–321 reaches 212–219 (to 312).
    expect(t.label_asym_id.slice(0, 3)).toEqual(['C', 'D', 'E']);
    const at = t.beg_label_seq_id.indexOf(9);
    expect([t.label_asym_id[at], t.end_label_seq_id[at]]).toEqual(['C', 199]);
    expect(t.beg_label_seq_id).toContain(212);
    expect(colorsOf(t).length).toBe(2);
    ctx.dispose();
  });

  it('switching the coloring recomposes: the model by pLDDT, by its domains, by pLDDT again (the same MVS)', async () => {
    const ctx = await viewer();
    let shown: Shown | null = null;
    const states: string[] = [];
    const [show] = [
      (s: Shown) => {
        shown = s;
        const r = composeStructureView(ctx, s.structure, s.colorBy ? { tracks: [s.colorBy] } : {});
        states.push(JSON.stringify(r?.state.root));
      },
    ];
    const actions = structureTrackActions({
      structures: P53_STRUCTURES.map((s) => ({ label: s.label, structure: ref(s.label) })),
      shown: () => shown,
      show,
    });
    const colorBy = actions.find((a) => a.id === 'molstar:color-by');
    if (!colorBy) throw new Error('no action');
    // "Show AlphaFold in 3D" is offered on the coverage summary and on the model's own row, not on 1TUP's.
    const showAf = actions.find((a) => a.label === 'Show AlphaFold in 3D');
    expect(showAf?.when?.(ctx, 'main', P53_TRACKS.coverage, null)).toBe(true);
    expect(showAf?.when?.(ctx, 'main', `${P53_TRACKS.coverage}-af-p04637-f1`, null)).toBe(true);
    expect(showAf?.when?.(ctx, 'main', `${P53_TRACKS.coverage}-1tup`, null)).toBe(false);
    showAf?.run(ctx, 'main', P53_TRACKS.coverage, null);
    expect(colorBy.checked?.(ctx, 'main', P53_TRACKS.plddt, null)).toBe(true);
    colorBy.run(ctx, 'main', P53_TRACKS.domains, null);
    expect(colorBy.checked?.(ctx, 'main', P53_TRACKS.domains, null)).toBe(true);
    expect(colorBy.checked?.(ctx, 'main', P53_TRACKS.plddt, null)).toBe(false);
    colorBy.run(ctx, 'main', P53_TRACKS.plddt, null);
    expect(states.length).toBe(3);
    expect(states[1]).not.toBe(states[0]);
    expect(states[2]).toBe(states[0]);
    // On 1TUP, pLDDT is not offered; the domains are.
    show({ structure: ref('1TUP'), colorBy: null });
    expect(colorBy.when?.(ctx, 'main', P53_TRACKS.plddt, null)).toBe(false);
    expect(colorBy.when?.(ctx, 'main', P53_TRACKS.domains, null)).toBe(true);
    ctx.dispose();
  });

  it("the selection becomes a labeled component; the hover card's action offers itself on items with a structure", async () => {
    const ctx = await viewer();
    ctx.selection.set({ systems: ['seqstar:p53/unit'], runs: Int32Array.of(0, 174, 175) }, 'api');
    const r = composeMvs(ctx, itemOf(ctx, '1TUP chain A'), { selectionLabel: 'R175' });
    expect(nodes(r?.state.root as Node, 'label').map((n) => n.params?.text)).toEqual(['R175']);
    const component = nodes(r?.state.root as Node, 'component').find((n) =>
      JSON.stringify(n.params?.selector).includes('beg_label_seq_id'),
    );
    expect(component?.params?.selector).toContainEqual({
      label_asym_id: 'C',
      beg_label_seq_id: 82,
      end_label_seq_id: 82,
    });
    const shown: Shown[] = [];
    const action = showInStructureAction((s) => shown.push(s));
    expect(action.when(ctx, itemOf(ctx, '2OCJ chain A'))).toBe(true);
    expect(action.when(ctx, { annotation: 'features', item: 0 })).toBe(false);
    action.run(ctx, itemOf(ctx, '2OCJ chain A'), 'api');
    expect(shown[0]?.structure.source).toEqual({ service: 'pdbe', id: '2ocj' });
    ctx.dispose();
  });

  it("a sites track (PTMs) shows its residues as ball and stick, in its kinds' colors; a ranges track does not", async () => {
    const ctx = await viewer();
    const r = composeStructureView(ctx, ref('AlphaFold'), { tracks: [P53_TRACKS.ptms] });
    if (!r) throw new Error('nothing composed');
    // The PTM sites of the trimmed UniProt entry (modified residues and cross-links), one residue each.
    const ptms = (
      JSON.parse(readFileSync(new URL('./fixtures/p53/P04637.json', import.meta.url), 'utf8')) as {
        features: { type: string; location: { start: { value: number } } }[];
      }
    ).features
      .filter((f) => f.type === 'Modified residue' || f.type === 'Cross-link')
      .map((f) => f.location.start.value);
    expect(r.sticks.length).toBe(1);
    const residues = new Set<number>();
    const t = r.sticks[0];
    for (let k = 0; k < t.color.length; k++)
      for (let i = t.beg_label_seq_id[k]; i <= t.end_label_seq_id[k]; i++) residues.add(i);
    expect([...residues].sort((a, b) => a - b)).toEqual([...new Set(ptms)].sort((a, b) => a - b));
    // Each residue in its modification's color (the fixture's `kind`; where several share a residue, the last row's).
    const kinds = new Map<number, string>();
    for (const line of readFileSync(new URL('./fixtures/p53/P04637-ptms.tsv', import.meta.url), 'utf8')
      .trim()
      .split('\n')
      .slice(1)) {
      const [pos, kind] = line.split('\t');
      kinds.set(Number(pos), kind);
    }
    const drawn = new Map<number, string>();
    for (let k = 0; k < t.color.length; k++)
      for (let i = t.beg_label_seq_id[k]; i <= t.end_label_seq_id[k]; i++) drawn.set(i, t.color[k]);
    for (const [pos, kind] of kinds) expect([pos, drawn.get(pos)]).toEqual([pos, PTM_COLORS[kind]]);
    expect(drawn.get(15)).toBe(PTM_COLORS.phosphorylation);
    expect(drawn.get(386)).toBe(PTM_COLORS.sumoylation);
    const fromUri = nodes(r.state.root as Node, 'component_from_uri');
    expect(fromUri.length).toBe(1);
    expect(nodes(fromUri[0], 'representation')[0].params?.type).toBe('ball_and_stick');
    expect(tableOf(nodes(fromUri[0], 'color_from_uri')[0].params?.uri as string)).toEqual(t);
    expect(MVSData.validationIssues(r.state)).toBeUndefined();
    const domains = composeStructureView(ctx, ref('AlphaFold'), { tracks: [P53_TRACKS.domains] });
    expect(domains?.sticks).toEqual([]);
    ctx.dispose();
  });
});

describe('a split row applied to 3D', () => {
  it('the Regions track by one of its rows: Oligomerization alone colors 1C26, titled by the row', async () => {
    const ctx = await viewer();
    ctx.ui.toggleGroup(P53_TRACKS.regions);
    await ctx.graph.update();
    const view = ctx.graph.assembled.value?.view;
    if (!view) throw new Error('no view');
    const k = view.resolved.tracks.split[view.resolved.tracks.ref.indexOf(view.set.strings.find(P53_TRACKS.regions))];
    const labels = (view.prepared[k]?.state as { split?: { labels: string[] } } | undefined)?.split?.labels ?? [];
    const row = labels.indexOf('Oligomerization');
    expect(row).toBeGreaterThanOrEqual(0);
    expect(splitRow(view, P53_TRACKS.regions, row)).toEqual({ field: 'label', value: 'Oligomerization' });
    const r = composeStructureView(ctx, ref('1C26'), {
      tracks: [P53_TRACKS.regions],
      where: { field: 'label', value: 'Oligomerization' },
    });
    // P04637 325–356 is all of 1C26's chain A (label_seq_id 1–32), in one color.
    expect(r?.tables).toEqual([
      { label_asym_id: ['A'], beg_label_seq_id: [1], end_label_seq_id: [32], color: [expect.any(String)] },
    ]);
    const whole = composeStructureView(ctx, ref('1C26'), { tracks: [P53_TRACKS.regions] });
    expect(new Set(whole?.tables[0]?.color).size).toBeGreaterThan(1);
    // The track menu on that row's line: that row alone, checked there and not on the track's own line.
    let menuShown: Shown | null = { structure: ref('1C26'), colorBy: null };
    const applyTo3d = structureTrackActions({
      structures: [],
      shown: () => menuShown,
      show: (x) => {
        menuShown = x;
      },
    }).find((x) => x.id === 'molstar:color-by');
    if (!applyTo3d) throw new Error('no action');
    expect(applyTo3d.labelFor?.(ctx, 'main', P53_TRACKS.regions, 'Oligomerization')).toBe(
      'Apply Oligomerization to 3D',
    );
    expect(applyTo3d.when?.(ctx, 'main', P53_TRACKS.regions, 'Oligomerization')).toBe(true);
    applyTo3d.run(ctx, 'main', P53_TRACKS.regions, 'Oligomerization');
    expect(menuShown).toMatchObject({ where: { field: 'label', value: 'Oligomerization' } });
    expect(applyTo3d.checked?.(ctx, 'main', P53_TRACKS.regions, 'Oligomerization')).toBe(true);
    expect(applyTo3d.checked?.(ctx, 'main', P53_TRACKS.regions, null)).toBe(false);
    const applied = menuShown as Shown | null;
    const fromMenu = composeStructureView(ctx, ref('1C26'), {
      tracks: [applied?.colorBy ?? ''],
      where: applied?.where ?? null,
    });
    expect(fromMenu?.tables).toEqual(r?.tables);
    applyTo3d.run(ctx, 'main', P53_TRACKS.regions, null);
    expect(applyTo3d.checked?.(ctx, 'main', P53_TRACKS.regions, null)).toBe(true);
    // The hover card's action on a region item, hovered on that row.
    let shown: Shown | null = { structure: ref('1C26'), colorBy: null };
    const action = colorByItemAction(
      () => shown,
      (s) => {
        shown = s;
      },
    );
    const a = view.set.annotations.find((x) => x && view.set.strings.get(x.key) === 'features');
    const at = Array.from({ length: a?.count ?? 0 }, (_, i) => i).find(
      (i) => view.set.strings.get(a?.rowLabel?.[i] ?? -1) === 'Oligomerization',
    );
    const target = { annotation: 'features', item: a?.item[at ?? -1] ?? -1, track: P53_TRACKS.regions, split: row };
    expect(action.when(ctx, target)).toBe(true);
    expect(action.titleFor?.(ctx, target)).toBe('Color 3D by Oligomerization');
    action.run(ctx, target, 'pointer');
    expect(shown).toMatchObject({ colorBy: P53_TRACKS.regions, where: { field: 'label', value: 'Oligomerization' } });
    ctx.dispose();
  });
});

describe('molstar:sifts_depth', () => {
  it("the summary's depth: runs of positions with the structures covering them, counted and named", async () => {
    const { set } = await compile(p53DemoDoc(P53.sequence), registries(), demoLoaders());
    const byKey = (k: string) => set.annotations.find((a) => a && set.strings.get(a.key) === k);
    const depth = byKey('coverage-depth');
    if (!depth) throw new Error('no coverage');
    // The coverage's rows in P04637's numbering (the demo's own sit on the chains' units).
    const cov = coverageRows(set, p53CoverageOptions);
    const counts = new Map<number, Set<string>>();
    for (let r = 0; r < cov.beg.length; r++) {
      for (let i = cov.beg[r] - 1; i < cov.end[r]; i++) {
        const names = counts.get(i) ?? new Set<string>();
        names.add(cov.entry[r]);
        counts.set(i, names);
      }
    }
    // The depth's runs, position by position.
    const drawn = new Map<number, number>();
    for (let r = 0; r < depth.count; r++)
      for (let i = depth.start[r]; i < depth.stop[r]; i++) drawn.set(i, depth.sparse?.[r] ?? Number.NaN);
    expect(drawn.size).toBe(counts.size);
    for (const [i, names] of counts) expect(drawn.get(i)).toBe(names.size);
    // P04637 200 (index 199): 1TUP, 2OCJ, 2XWR and the model; 340: 1C26 and the model; 10: the model alone.
    expect([drawn.get(199), drawn.get(339), drawn.get(9)]).toEqual([4, 2, 1]);
    const at = Array.from({ length: depth.count }, (_, r) => r).find(
      (r) => depth.start[r] <= 339 && 339 < depth.stop[r],
    );
    const tip = depth.tooltip && at !== undefined ? set.strings.get(depth.tooltip[at]) : null;
    expect(tip).toBe('covered by 2 structures: 1C26, AlphaFold');
  });
});
