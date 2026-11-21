/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { DataStore } from '@molstar/seqx-core/data';
import {
  attachSources,
  bundledTextMetrics,
  type Frame,
  hitView,
  LightTokens,
  Renderer,
  renderConfig,
} from '@molstar/seqx-core/render';
import { BuiltInLayerTypes, BuiltInLayouts } from '@molstar/seqx-core/render/built-ins';
import { resolveSpec } from '@molstar/seqx-core/resolve';
import { projectOntoView, resolveViews } from '@molstar/seqx-core/view';
import { describe, expect, it } from 'vitest';
import { build, Cds, Domain, Exons, Gene, protein } from '../coding-dna.js';

// The headless twin of Features/Coordinates/Protein Over Its Coding DNA: the drawn marks against the codon map, every
// residue over its three nucleotides (a codon split by the intron in two parts), no protein letter over the intron, the
// domain's pieces ending where its residues' codons end, and a hovered residue lighting its codon. It fails when the
// codon alignment is off by one nucleotide (a nucleotide is 6.6 px at this width: the tolerance is a pixel).

const Size = { width: 900, height: 200 };
const Tracks = { ruler: 0, dna: 1, protein: 2, exons: 3, domain: 4, site: 5 } as const;

/** The genomic runs of coding indices `[a, b)` through the map (the exons' coding pieces, in order). */
function genomic(a: number, b: number): [number, number][] {
  const out: [number, number][] = [];
  let at = 0;
  for (let i = 0; i < Cds.start.length; i++) {
    const len = Cds.end[i] - Cds.start[i];
    const lo = Math.max(a, at);
    const hi = Math.min(b, at + len);
    if (hi > lo) out.push([Cds.start[i] + lo - at, Cds.start[i] + hi - at]);
    at += len;
  }
  return out;
}

async function frameOf(emphasis?: (f: Frame, r: Renderer, model: ReturnType<typeof resolveSpec>) => object) {
  const model = resolveSpec(build(), {});
  const view = resolveViews(model)[0];
  const store = new DataStore();
  const renderer = new Renderer({
    types: BuiltInLayerTypes,
    layouts: BuiltInLayouts,
    model,
    view,
    data: store,
    tokens: LightTokens,
    text: bundledTextMetrics(),
    config: renderConfig(),
  });
  renderer.frame({ size: Size, settled: false });
  await attachSources(model, store);
  let frame = renderer.frame({ size: Size, settled: true });
  for (let i = 0; i < 4; i++) {
    store.want(frame.wants);
    await store.whenSettled();
    frame = renderer.frame({ size: Size, settled: true });
    if (!frame.rebuilt) break;
  }
  const first = frame;
  if (emphasis) frame = renderer.frame({ size: Size, settled: true, emphasis: emphasis(first, renderer, model) });
  store.retire();
  return { f: frame, model, view };
}

/** The rects of a pass within a track's band (an outline starts a few px outside its mark), by x. */
function rectsIn(f: Frame, pass: 'data' | 'hover', track: number) {
  const lay = f.layout.sections[0].tracks[track].layers[0].lay as { box: { y: number; height: number } };
  const l = f.lists.passes[pass].rects;
  const out: { x0: number; x1: number; y: number; h: number; fill: number }[] = [];
  for (let i = 0; i < l.n; i++)
    if (l.y[i] >= lay.box.y - 4 && l.y[i] < lay.box.y + lay.box.height)
      out.push({ x0: l.x[i], x1: l.x[i] + l.w[i], y: l.y[i], h: l.h[i], fill: l.fill[i] });
  return out.sort((a, b) => a.x0 - b.x0);
}

const xOf = (f: Frame, p: number) =>
  (f.layout.sections[0].tracks[0].layers[0].lay as { axis: { x(p: number): number } }).axis.x(p);

/** Whether a mark `[x0, x1]` is a run `[a, b)` of the gene within a pixel. */
const at = (f: Frame, m: { x0: number; x1: number }, a: number, b: number) =>
  Math.abs(m.x0 - xOf(f, a)) <= 1 && Math.abs(m.x1 - xOf(f, b)) <= 1;

describe('A protein over its coding DNA, drawn against the map', () => {
  it('the protein is what the coding pieces translate to, the stop excluded', () => {
    expect(protein()).toBe('MASKELTPRDGHWYFNSQVE');
    expect(Gene.slice(Cds.end[1] - 3, Cds.end[1])).toBe('TAA');
  });

  it('every residue is drawn over its three nucleotides, a codon split by the intron in two parts, none over the intron', async () => {
    const { f } = await frameOf();
    const cells = rectsIn(f, 'data', Tracks.protein);
    // 20 residues: 19 one part, the 12th (index 11) in two: 21 cells
    const expected: [number, number][] = [];
    for (let r = 0; r < protein().length; r++) expected.push(...genomic(3 * r, 3 * r + 3));
    expect(expected.length).toBe(21);
    expect(cells.length).toBe(expected.length);
    expected.forEach(([a, b], i) => expect(at(f, cells[i], a, b), `cell ${i}: ${a}-${b}`).toBe(true));
    // the split codon: one nucleotide before the intron, two after it
    expect(genomic(33, 36)).toEqual([
      [45, 46],
      [76, 78],
    ]);
    // the intron (46–76) has no protein letters
    expect(cells.filter((c) => c.x0 < xOf(f, 76) - 1 && c.x1 > xOf(f, 46) + 1)).toEqual([]);
  });

  it('the domain’s two pieces end where its residues’ codons end, and the connector joins them', async () => {
    const { f } = await frameOf();
    const marks = rectsIn(f, 'data', Tracks.domain).filter((m) => m.h > 2);
    const want = genomic(3 * Domain.start[0], 3 * Domain.end[0]);
    expect(want).toEqual([
      [24, 46],
      [76, 84],
    ]);
    expect(marks.length).toBe(2);
    want.forEach(([a, b], i) => expect(at(f, marks[i], a, b), `piece ${i}`).toBe(true));
    const line = rectsIn(f, 'data', Tracks.domain).filter((m) => m.h <= 1.5);
    expect(line.length).toBe(1);
    expect(Math.abs(line[0].x0 - marks[0].x1)).toBeLessThanOrEqual(1);
    expect(Math.abs(line[0].x1 - marks[1].x0)).toBeLessThanOrEqual(1);
  });

  it('the exons are ranges over their stretches, coloured by exon', async () => {
    const { f } = await frameOf();
    const marks = rectsIn(f, 'data', Tracks.exons);
    expect(marks.length).toBe(2);
    expect(at(f, marks[0], Exons.start[0], Exons.end[0])).toBe(true);
    expect(at(f, marks[1], Exons.start[1], Exons.end[1])).toBe(true);
    expect(marks[0].fill).not.toBe(marks[1].fill);
  });

  it('hovering a residue lights it and its codon, on both tracks (the split one in both parts)', async () => {
    const { f } = await frameOf((first, _r, model) => {
      const cell = (a: number, b: number) => (xOf(first, a) + xOf(first, b)) / 2;
      const track = first.layout.sections[0].tracks[Tracks.protein].layers[0].lay as {
        box: { y: number; height: number };
      };
      // the first nucleotide's part of residue 11, before the intron
      const hit = hitView(first.layout, cell(45, 46), track.box.y + track.box.height / 2)!;
      // the hit is on the gene's axis, its extent the residue on the protein (what the viewer lifts a hover by)
      expect(hit.position).toBe(45);
      expect(hit.extent?.parts[0].system).toBe('/data/protein');
      const projected = projectOntoView(model, resolveViews(model)[0], hit.extent!, {
        own: true,
        close: true,
      });
      return { hover: hit, projected: { hover: projected, hoverSection: '/views/main/gene' } };
    });
    // the section's emphasis is a column across its tracks (`emphasis.selectionShape`): the codon's two parts, each a
    // rect from the ruler's top down through the DNA's and the protein's bands
    const all = rectsIn(f, 'hover', Tracks.dna).concat(rectsIn(f, 'hover', Tracks.protein), rectsIn(f, 'hover', 0));
    const lay = (t: number) => f.layout.sections[0].tracks[t].layers[0].lay as { box: { y: number; height: number } };
    const dna = lay(Tracks.dna).box;
    const prot = lay(Tracks.protein).box;
    const over = (a: number, b: number, box: { y: number; height: number }) =>
      all.some(
        (m) =>
          Math.abs(m.x0 - xOf(f, a)) <= 4 &&
          Math.abs(m.x1 - xOf(f, b)) <= 4 &&
          m.y <= box.y &&
          m.y + m.h >= box.y + box.height,
      );
    for (const [a, b] of [
      [45, 46],
      [76, 78],
    ])
      expect(over(a, b, dna) && over(a, b, prot), `codon part ${a}-${b} lit over both tracks`).toBe(true);
  });
});
