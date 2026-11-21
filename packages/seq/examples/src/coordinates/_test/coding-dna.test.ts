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

async function frameOf(
  emphasis?: (f: Frame, r: Renderer, model: ReturnType<typeof resolveSpec>) => object,
  config = {},
) {
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
    config: renderConfig(config),
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

  /** A frame with the pointer on the middle of a track at a gene position, the hover lifted as the viewer lifts it. */
  const hovered = (track: number, position: number, config = {}) =>
    frameOf((first, _r, model) => {
      const lay = first.layout.sections[0].tracks[track].layers[0].lay as { box: { y: number; height: number } };
      const x = (xOf(first, position) + xOf(first, position + 1)) / 2;
      const hit = hitView(first.layout, x, lay.box.y + lay.box.height / 2)!;
      const projected = hit.extent
        ? projectOntoView(model, resolveViews(model)[0], hit.extent, { own: true, close: true })
        : undefined;
      return {
        hover: hit,
        ...(projected ? { projected: { hover: projected, hoverSection: '/views/main/gene' } } : {}),
      };
    }, config);

  /** Whether the hover pass has a rect over a codon part `[a, b)` inside a track's band (an outline's top is outside). */
  const litOn = (f: Frame, track: number, a: number, b: number) =>
    rectsIn(f, 'hover', track).some((m) => Math.abs(m.x0 - xOf(f, a)) <= 4 && Math.abs(m.x1 - xOf(f, b)) <= 4);

  it.each([
    ['column (the default)', {}],
    ['tracks', { emphasis: { selectionShape: 'tracks' as const } }],
  ])('hovering a residue lights its codon on the DNA track only, both parts of a split one (%s)', async (_, config) => {
    // the first nucleotide of residue 11's codon, before the intron
    const { f } = await hovered(Tracks.protein, 45, config);
    expect(litOn(f, Tracks.protein, 45, 46)).toBe(true);
    expect(litOn(f, Tracks.dna, 45, 46)).toBe(true);
    expect(litOn(f, Tracks.dna, 76, 78)).toBe(true);
    // never a column: the ruler, the exons, the domain and the site stay as they were
    for (const t of [Tracks.ruler, Tracks.exons, Tracks.domain, Tracks.site])
      expect(rectsIn(f, 'hover', t).filter((m) => m.h > 0 && m.x0 < xOf(f, 80) && m.x1 > xOf(f, 40))).toEqual([]);
  });

  it('hovering an item through the path (the domain) lights its outline and its connector, nothing else', async () => {
    const { f } = await hovered(Tracks.domain, 60);
    const lit = [Tracks.ruler, Tracks.dna, Tracks.protein, Tracks.exons, Tracks.site].filter(
      (t) => rectsIn(f, 'hover', t).length > 0,
    );
    expect(lit).toEqual([]);
    // its two pieces outlined and the connector between them in the hover ink
    expect(rectsIn(f, 'hover', Tracks.domain).length).toBe(3);
  });
});
