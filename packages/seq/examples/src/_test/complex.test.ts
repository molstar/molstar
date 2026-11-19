/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { readFileSync } from 'node:fs';
import { join } from 'node:path';
import type { Node } from '@molstar/core/repr/tree';
import { createSeqContext } from '@molstar/seqx-core/context';
import { DataStore, decodeSeq, type RangesList, SeqRuntime } from '@molstar/seqx-core/data';
import { RecordsRuntime } from '@molstar/seqx-core/data/records';
import { location } from '@molstar/seqx-core/query';
import {
  attachSources,
  type BlockPlacement,
  bundledTextMetrics,
  type GridLay,
  GridPlacement,
  gridSummary,
  hitView,
  layoutView,
  LightTokens,
  recordDrawLists,
} from '@molstar/seqx-core/render';
import { resolveSpec } from '@molstar/seqx-core/resolve';
import { resolveViews } from '@molstar/seqx-core/view';
import { describe, expect, it } from 'vitest';
import * as blocks from '../complex/blocks.js';
import * as complexGrid from '../complex/complex-grid.js';
import {
  Cdk2Features,
  complexPae,
  complexPaeDocument,
  constructOffset,
  Data,
  Files,
  Lengths,
} from '../complex/data.js';
import { cellReadout, rectReadout } from '../complex/readout.js';
import { nodeFixtureAssets, PackageFolder } from '../node-assets.js';
import { ComplexViews } from '../views.js';

// The Complex dataset against its data (its provenance is `sources.test.ts`'s): the derived data held by the views equal to the files, the anchors of the entries (CDK2's K33, E51,
// D145; T160 in the activation loop; the PSTAIRE helix at the interface of both copies; ATP at the hinge), the PAE's
// shape and the grid over the complex (a strip a pair of chains, a cell's hit on both chains' systems, the chain
// boundaries drawn, a rect's summary).

const Folder = join(PackageFolder, 'fixtures', 'complex');

interface SourceRow {
  readonly file: string;
}

const sources = JSON.parse(readFileSync(join(Folder, 'sources.json'), 'utf8')) as SourceRow[];
const read = (name: string) => readFileSync(join(Folder, name), 'utf8');

/** A `seq` file's ranges table on a system, its fields by name. */
function rangesOf(text: string, system: string) {
  const d = decodeSeq(text).data!;
  const pair = d.pairs.find((p) => p.request.system === system)!;
  const t = (pair.answer as { tables: { ranges: RangesList; letters?: { letters: string } } }).tables;
  const field = (name: string) => {
    const c = t.ranges.fields!.find((f) => f.name === name)!.column as unknown as {
      length: number;
      get?: (i: number) => unknown;
    } & ArrayLike<unknown>;
    return Array.from({ length: c.length }, (_, i) => (c.get ? c.get(i) : c[i]));
  };
  return { ranges: t.ranges, letters: t.letters?.letters, field };
}

describe('the Complex fixtures', () => {
  it('the letters the views hold are the files’: 1FIN’s entities, UniProt’s entries', () => {
    const fin = read('1FIN.seq');
    expect(rangesOf(fin, 'cdk2').letters).toBe(Data.letters.cdk2);
    expect(rangesOf(fin, 'cyclin-a').letters).toBe(Data.letters.construct);
    expect(JSON.parse(read('P24941.json')).sequence.value).toBe(Data.letters.P24941);
    expect(JSON.parse(read('P20248.json')).sequence.value).toBe(Data.letters.P20248);
  });

  it('the contacts the views hold are `contacts.seq`’s, row for row', () => {
    const c = rangesOf(read('contacts.seq'), 'cdk2');
    const rows: string[] = [];
    for (let i = 0; i < c.ranges.start.length; i++)
      rows.push(
        `${c.field('with')[i]}|${c.field('entry')[i]}|${c.field('copy')[i]}|${c.ranges.start[i]}|${c.field('partner')[i]}`,
      );
    const held: string[] = [];
    for (const [copy, a, b] of Data.contacts) held.push(`cyclin-a|1FIN|${copy}|${a}|${b}`);
    for (const l of Data.ligands) {
      const copy = l.entry === '1FIN' ? Data.chains.find((x) => x.id === l.chain)!.copy : 1;
      for (const [r, atom] of l.contacts) held.push(`${l.comp}|${l.entry}|${copy}|${r}|${atom}`);
    }
    expect(rows.sort()).toEqual(held.sort());
  });
});

describe('the Complex anchors', () => {
  const cdk2 = Data.letters.cdk2;

  it('CDK2: K33, E51 and D145 (DFG); PSTAIRE at 45–51; T160 in the activation loop; 1FIN’s CDK2 is UniProt’s', () => {
    expect([cdk2[32], cdk2[50], cdk2[144]]).toEqual(['K', 'E', 'D']);
    expect(cdk2.slice(44, 51)).toBe('PSTAIRE');
    expect(cdk2.slice(144, 147)).toBe('DFG');
    const loop = Cdk2Features.id.indexOf('aloop');
    expect(Cdk2Features.start[loop]).toBeLessThanOrEqual(159);
    expect(Cdk2Features.end[loop]).toBeGreaterThan(159);
    expect(cdk2[159]).toBe('T');
    expect(Data.letters.P24941).toBe(cdk2);
  });

  it('cyclin A: the construct is P20248 from 173 (its auth numbering), 260 of 432 residues', () => {
    expect(Data.cyclinFirst).toBe(173);
    expect(Data.letters.P20248.length).toBe(Lengths.cyclin);
    expect(Data.letters.construct.length).toBe(Lengths.construct);
    expect(Data.letters.P20248.slice(constructOffset(), constructOffset() + Lengths.construct)).toBe(
      Data.letters.construct,
    );
  });

  it('two copies (A with B, C with D); the PSTAIRE helix is at the interface in both', () => {
    expect(Data.chains.map((c) => `${c.id}${c.copy}`)).toEqual(['A1', 'B1', 'C2', 'D2']);
    for (const copy of [1, 2]) {
      const at = new Set(Data.contacts.filter((c) => c[0] === copy).map((c) => c[1]));
      // the helix 45–56 (0-based 44–55): most of its turns touch cyclin A
      expect([...at].filter((p) => p >= 44 && p < 56).length).toBeGreaterThanOrEqual(5);
    }
  });

  it('each copy’s ATP touches the hinge (E81, L83) and the glycine loop; K33 is just past 4 Å (the record’s fact)', () => {
    const atps = Data.ligands.filter((l) => l.entry === '1FIN');
    expect(atps.map((l) => [l.comp, l.chain, l.atoms.length])).toEqual([
      ['ATP', 'A', 31],
      ['ATP', 'C', 31],
    ]);
    for (const l of atps) {
      const at = new Set(l.contacts.map(([r]) => r + 1));
      expect(at.has(81) && at.has(83)).toBe(true);
      expect([10, 11, 12, 13].some((p) => at.has(p))).toBe(true);
      expect(at.has(33)).toBe(false);
    }
    // the other entries' ligands all sit at the hinge (L83)
    for (const l of Data.ligands.filter((x) => x.entry !== '1FIN'))
      expect(new Set(l.contacts.map(([r]) => r + 1)).has(83)).toBe(true);
  });
});

describe('the Complex’s synthetic PAE', () => {
  const n = Lengths.cdk2 + Lengths.cyclin;

  it('is AlphaFold-shaped, read through the matrix read: N × N over the chains in order, whole Å', () => {
    const doc = complexPaeDocument();
    expect(doc[0].predicted_aligned_error.length).toBe(n);
    const { matrix } = complexPae(doc);
    expect([matrix.rows, matrix.cols]).toEqual([n, n]);
    expect(matrix.values.every((v) => Number.isInteger(v) && v >= 0 && v <= 31)).toBe(true);
    // a wrong size is refused by the read (its axes' total checked)
    expect(() => complexPae([{ predicted_aligned_error: doc[0].predicted_aligned_error.slice(1) }])).toThrow();
  });

  it('low across 1FIN’s interface, high where cyclin A2 is disordered, low within a folded chain', () => {
    const { block } = complexPae();
    const across = block('cdk2', 'cyclin');
    const ny = Lengths.cyclin;
    const shift = constructOffset();
    let atInterface = 0;
    for (const [, a, b] of Data.contacts) atInterface += across[a * ny + b + shift];
    atInterface /= Data.contacts.length;
    let disordered = 0;
    for (let i = 0; i < Lengths.cdk2; i++) for (let j = 0; j < 150; j++) disordered += across[i * ny + j];
    disordered /= Lengths.cdk2 * 150;
    const within = block('cdk2', 'cdk2');
    let folded = 0;
    for (let i = 20; i < 280; i++) folded += within[i * Lengths.cdk2 + i + 3];
    folded /= 260;
    expect(atInterface).toBeLessThan(7);
    expect(disordered).toBeGreaterThan(25);
    expect(folded).toBeLessThan(6);
  });
});

/** A view laid out headless over the fixtures, its sources opened (the records runtime for its JSON). */
async function laid(spec: Node, size = { width: 720, height: 900 }) {
  const model = resolveSpec(spec, { runtimes: [SeqRuntime, RecordsRuntime] });
  const v = resolveViews(model)[0];
  const store = new DataStore();
  await attachSources(model, store, { assets: nodeFixtureAssets().reader });
  const input = { model, view: v, size, viewports: {}, tokens: LightTokens, text: bundledTextMetrics(), data: store };
  for (let i = 0; i < 3; i++) {
    store.want(layoutView(input).wants);
    await store.whenSettled();
  }
  return { layout: layoutView(input), model, view: v };
}

describe('1FIN’s blocks', () => {
  it('copies merged ×2, the ATPs after the polymers; split, a block a chain and a ligand', async () => {
    const { layout } = await laid(blocks.build());
    expect(layout.sections[0].section.blocks!.map((b) => b.label)).toEqual([
      'CDK2 (A) ×2',
      'cyclin A (B) ×2',
      'ATP ×2',
    ]);
    expect(layout.sections[0].section.axis.ranges).toEqual([0, Lengths.cdk2 + Lengths.construct + 1]);
    const split = await laid(blocks.split());
    expect(split.layout.sections[0].section.blocks!.length).toBe(6);
  });

  it('each feature on its own chain: the PSTAIRE helix under CDK2’s block only', async () => {
    const { layout } = await laid(blocks.build());
    const p = layout.sections[0].placement as BlockPlacement;
    const st = p.strips;
    const keysAt = (strip: number, pos: number) => {
      const a = st.axes[strip];
      const keys = new Set<string>();
      for (let y = st.bandY[strip * st.bands]; y < st.bandY[strip * st.bands] + st.bandH[strip * st.bands]; y += 2) {
        const h = hitView(layout, a.xOfUnit(pos - st.first[strip] + 0.5), y);
        if (h?.element) keys.add(String(h.element.key));
      }
      return [...keys];
    };
    expect(keysAt(0, 50)).toEqual(['pstaire']);
    // cyclin A's construct position 40 is P20248's 213: its first cyclin box
    expect(keysAt(1, Lengths.cdk2 + 40)).toEqual(['box-n']);
  });
});

describe('the grid over the complex', () => {
  const view = ComplexViews.find((v) => v.id === 'complex-grid')!;

  it('a strip a pair of chains, each pair’s grid drawn on it alone; the plot’s x is the blocks’ x', async () => {
    const { layout } = await laid(view.build());
    expect(layout.issues).toEqual([]);
    const top = layout.sections[0].placement.strips;
    const g = layout.sections[1];
    const st = g.placement.strips;
    expect(g.placement).toBeInstanceOf(GridPlacement);
    expect(st.n).toBe(4);
    const systems = (i: number) => `${st.axes[i].system} × ${st.ys[i]!.system}`;
    expect([0, 1, 2, 3].map(systems)).toEqual([
      '/data/cdk2-protein × /data/cdk2-protein',
      '/data/cdk2-protein × /data/cyclin-protein',
      '/data/cyclin-protein × /data/cdk2-protein',
      '/data/cyclin-protein × /data/cyclin-protein',
    ]);
    // each strip's one grid drawn, the other three not; x as the blocks above it
    for (let i = 0; i < 4; i++) {
      const drawn = g.lays[i].filter((l) => (l.lay as GridLay).grid);
      expect(drawn.length, systems(i)).toBe(1);
      expect([st.axes[i].x0, st.axes[i].width]).toEqual([top.axes[i >> 1].x0, top.axes[i >> 1].width]);
    }
  });

  it('a cell’s hit is on both chains’ systems; the boundaries are hairlines across the plot', async () => {
    const { layout, model, view: v } = await laid(view.build());
    const st = layout.sections[1].placement.strips;
    // CDK2's K33 against cyclin A2's 263 (strip 1: CDK2 across, cyclin A2 down)
    const ax = st.axes[1];
    const ay = st.ys[1]!;
    ax.cell(32);
    const x = (ax.sx0 + ax.sx1) / 2;
    ay.cell(262);
    const y = (ay.sx0 + ay.sx1) / 2;
    const h = hitView(layout, x, y)!;
    expect([h.system, h.position, h.ySystem, h.yPosition, h.target]).toEqual([
      '/data/cdk2-protein',
      32,
      '/data/cyclin-protein',
      262,
      'cell',
    ]);
    expect(cellReadout({ model, view: v }, h as never)).toMatch(/^CDK2 K33 × cyclin A2 [A-Z]263: \d+(\.\d)? Å$/);
  });

  it('the chain boundaries: a hairline down the plot and one across it, in canvas lists and the SVG', async () => {
    const ctx = createSeqContext({
      registry: [
        { kind: 'format-runtime', runtime: SeqRuntime },
        { kind: 'format-runtime', runtime: RecordsRuntime },
      ],
      ports: { assets: () => nodeFixtureAssets() },
    });
    const r = await ctx.renderSvg(view.build(), { size: { width: 720, height: 900 } });
    expect(r.issues).toEqual([]);
    const fills = recordDrawLists(r.lists, ['labels'])
      .split('\n')
      .filter((l) => l.startsWith('rect '));
    const parsed = fills.map((l) => l.split(' ').slice(1, 5).map(Number));
    const down = parsed.filter(([, , w, h]) => w === 1 && h > 100);
    const across = parsed.filter(([, , w, h]) => h === 1 && w > 100);
    expect([down.length, across.length]).toEqual([1, 1]);
    // square at aspect 1: the plot as tall as wide, each hairline the plot's whole extent
    expect(down[0][3]).toBe(across[0][2]);
    expect(r.svg).toContain(`width="1" height="${down[0][3]}"`);
    ctx.retire();
  });

  it('gridSummary over a rect across the boundary: the cells of each pair’s grid inside it, counted once', async () => {
    const { model, view: v } = await laid(view.build());
    const a = '/data/cdk2-protein';
    const c = '/data/cyclin-protein';
    const rect = location({
      rects: [
        { a: { system: a, ranges: [290, 298] }, b: { system: a, ranges: [0, 10] } },
        { a: { system: c, ranges: [0, 5] }, b: { system: a, ranges: [0, 10] } },
      ],
    });
    const s = gridSummary({ model, view: v }, rect);
    expect(s.count).toBe(8 * 10 + 5 * 10);
    const { block } = complexPae();
    const aa = block('cdk2', 'cdk2');
    const ca = block('cyclin', 'cdk2');
    let sum = 0;
    let min = Infinity;
    let max = -Infinity;
    const add = (x: number) => {
      sum += x;
      min = Math.min(min, x);
      max = Math.max(max, x);
    };
    for (let i = 290; i < 298; i++) for (let j = 0; j < 10; j++) add(aa[i * Lengths.cdk2 + j]);
    for (let i = 0; i < 5; i++) for (let j = 0; j < 10; j++) add(ca[i * Lengths.cdk2 + j]);
    expect(s.mean).toBeCloseTo(sum / 130, 5);
    expect([s.min, s.max]).toEqual([min, max]);
    expect(rectReadout({ model, view: v }, rect)).toBe(
      `Aligned CDK2 291–298, cyclin A2 1–5, Scored CDK2 1–10, mean ${(sum / 130).toFixed(1)} Å`,
    );
    expect(gridSummary({ model, view: v }, location({})).count).toBe(0);
  });
});

describe('the Complex views over the grid builder', () => {
  it('the grid view builds the four blocks of chain pairs on their systems', () => {
    const spec = complexGrid.build();
    const grids = (spec.children ?? [])
      .find((c) => c.id === 'data')!
      .children!.filter((c) => c.kind === 'grid')
      .map((c) => c.id);
    expect(grids).toEqual(['pae-cdk2-cdk2', 'pae-cdk2-cyclin', 'pae-cyclin-cdk2', 'pae-cyclin-cyclin']);
  });

  it('the fixtures’ locations are the files', () => {
    for (const f of Object.values(Files)) expect(sources.some((r) => `fixtures/${r.file}` === f)).toBe(true);
  });
});
