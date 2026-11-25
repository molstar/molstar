/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The ruler over declared numberings (_spec/README.md "Render", Layer types, `ruler`; Numberings from a source): a
// heavy chain's Kabat and IMGT rulers draw their own labels, each the label of the position under its tick (insertion
// codes and IMGT's out-of-order CDR3 among them), in the draw lists and the SVG alike; a numbering whose labels are a
// records source's column draws nothing while they come, then its labels, and the query model reads them once held.

import { RecordsRuntime } from '@molstar/seqx-core/data/records';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { overlayQueryModel, EmptyOverlay } from '@molstar/seqx-core/resolve';
import { resolveAddressIn } from '@molstar/seqx-core/query';
import { describe, expect, it } from 'vitest';
import { H, Imgt, Kabat } from '../../resolve/_test/antibody.js';
import { type RulerLay, withSourceNumberings } from '../index.js';
import { renderSvg } from '../render-svg.js';
import { settle } from './helpers.js';
import { layerLay } from './small-helpers.js';

const Size = { width: 1400, height: 120 };

/** One section a numbering, a ruler track each; `source`: Kabat's labels are a records source's column. */
function spec(source = false) {
  const b = createSeqSpec({ title: 'rulers' });
  const records = source
    ? b.data.source('kabat-table', {
        format: 'json',
        system: 'h',
        value: Kabat.map((label, position) => ({ position, label })),
      })
    : undefined;
  const h = b.data.sequence('h', { alphabet: 'protein', letters: H }, (s) => {
    if (records) s.numbering('kabat', { method: 'labels', source: records, scheme: 'kabat' });
    else s.numbering('kabat', { method: 'labels', labels: Kabat, scheme: 'kabat' });
    s.numbering('imgt', { method: 'labels', labels: Imgt, scheme: 'imgt' });
  });
  b.views.view('main', (v) => {
    for (const n of ['kabat', 'imgt'])
      v.section(n, { of: h }).track('r', { label: n }).layer('r', { type: 'ruler', numbering: n });
  });
  return b.build();
}

/** The labels a ruler draws, each with the position under its tick. */
function shown(lay: RulerLay): { label: string; position: number }[] {
  const out: { label: string; position: number; x: number }[] = [];
  for (let j = 0; j < lay.n; j++)
    if (lay.shown[j]) out.push({ label: lay.labels[j], position: lay.axis.positionAt(lay.x[j]), x: lay.x[j] });
  return out.sort((a, b) => a.x - b.x).map(({ label, position }) => ({ label, position }));
}

describe('the ruler in a declared numbering', () => {
  it('draws Kabat’s and IMGT’s labels, each the label of the position under its tick', async () => {
    const s = await settle(spec());
    const f = await s.frame({ size: Size, viewport: { start: 0, end: H.length } });
    expect(f.layout.issues).toEqual([]);
    const kabat = shown(layerLay<RulerLay>(f, 0, 0, 0));
    const imgt = shown(layerLay<RulerLay>(f, 0, 0, 1));
    for (const { label, position } of kabat) expect(label).toBe(Kabat[position]);
    for (const { label, position } of imgt) expect(label).toBe(Imgt[position]);
    // over 48 px a position (the label gap): a step of 1, every label with room; insertion codes as minor ticks, the unnumbered tag none
    expect(kabat.map((l) => l.label)).toEqual(Kabat.filter((l) => l !== null));
    expect(imgt.map((l) => l.label)).toEqual(Imgt);
  });

  it('an insertion’s label where it clears what is drawn, though wider than its position: the nearest a major first', async () => {
    // about 24 px a position: `111.1` and `112.1` are wider than theirs
    const s = await settle(spec());
    const f = await s.frame({ size: { width: 640, height: 120 }, viewport: { start: 0, end: H.length } });
    const lay = layerLay<RulerLay>(f, 0, 0, 1);
    expect(lay.axis.scale).toBeLessThan(30);
    const imgt = shown(lay);
    for (const { label, position } of imgt) expect(label).toBe(Imgt[position]);
    const minors = imgt.filter((l) => l.label.includes('.')).map((l) => l.label);
    expect(minors.length).toBeGreaterThan(0);
    // wider than its position, drawn all the same where it clears its neighbours
    for (let j = 0; j < lay.n; j++)
      if (lay.shown[j] && lay.labels[j].includes('.')) expect(lay.widths[j]).toBeGreaterThan(lay.axis.scale);
    // never two drawn labels within the gap of each other
    const spans = Array.from({ length: lay.n }, (_, j) => j)
      .filter((j) => lay.shown[j])
      .map((j) => [lay.x[j] - lay.widths[j] / 2, lay.x[j] + lay.widths[j] / 2])
      .sort((a, b) => a[0] - b[0]);
    for (let k = 1; k < spans.length; k++) expect(spans[k][0]).toBeGreaterThanOrEqual(spans[k - 1][1] + 4 - 1e-6);
  });

  it('the SVG prints the same labels', async () => {
    const r = await renderSvg(spec(), { size: Size, viewport: { start: 0, end: H.length }, legend: 'none' });
    expect(r.issues).toEqual([]);
    const texts = [...r.svg.matchAll(/<text[^>]*>([^<]*)<\/text>/g)].map((m) => m[1]);
    for (const l of ['52A', '82C', '100B']) expect(texts).toContain(l);
    for (const l of ['111.1', '112.1', '113']) expect(texts).toContain(l);
  });
});

describe('a numbering from a records source', () => {
  it('draws nothing while its labels come, then its labels; the query model reads them once held', async () => {
    const s = await settle(spec(true), [RecordsRuntime]);
    const first = s.renderer.frame({ size: Size, viewport: { start: 0, end: H.length } });
    expect(shown(layerLay<RulerLay>(first, 0, 0, 0))).toEqual([]);
    // its tiles are wanted (the whole system at level 1)
    expect(first.wants.some((w) => w.source === '/data/kabat-table' && w.end === H.length)).toBe(true);
    const q = withSourceNumberings(overlayQueryModel(s.model, EmptyOverlay), s.model, s.store);
    // a form naming it waits for its labels, never says it is unknown
    const waiting = resolveAddressIn(s.model.spec, q, '/data/h/(at :label 82B :numbering kabat)');
    expect(waiting.status).toBe('waiting');
    expect(waiting.issues.map((i) => i.code)).toEqual(['seq.query.numbering.pending', 'seq.address.waiting']);

    const f = await s.frame({ size: Size, viewport: { start: 0, end: H.length } });
    expect(f.layout.issues).toEqual([]);
    const kabat = shown(layerLay<RulerLay>(f, 0, 0, 0));
    expect(kabat.map((l) => l.label)).toEqual(Kabat.filter((l) => l !== null));
    for (const { label, position } of kabat) expect(label).toBe(Kabat[position]);

    const held = withSourceNumberings(overlayQueryModel(s.model, EmptyOverlay), s.model, s.store);
    const r = resolveAddressIn(s.model.spec, held, '/data/h/(range :from 82A :to 83 :numbering kabat)');
    expect(r.result?.type === 'location' && Array.from(r.result.location.parts[0].ranges)).toEqual([11, 15]);
  });
});
