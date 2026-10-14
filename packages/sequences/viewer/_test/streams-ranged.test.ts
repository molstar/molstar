/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A ranged piece of a streamed chromosome (GM-01's regression test; seqspec:sources#the-viewers-side): a block
// `range: [A, E]` of chr2, forward and reversed, asks for and loads every region over the positions it shows, at the
// whole piece and after a zoom and a pan. The region math reads block indices, which for a unit's piece are the
// system's positions: no offset, and the margin clamped to the piece.

import { createSeqViewSpec } from '@molstar/seq-core/builder';
import type { Loaders } from '@molstar/seq-core/format';
import { describe, expect, it } from 'vitest';
import { CODING_CDS, hg38Streams } from '../../core/_test/format/stream-provider.js';
import { createSeqViewer, DefaultRegistry, DefaultSeqViewerSpec, manualScheduler, RegionState } from '../src/index.js';
import type { Streams } from '../src/streams.js';
import { fixtureMeasurer } from '../src/testing/index.js';

const A = CODING_CDS[0][0] - 5000;
const E = CODING_CDS[2][1] + 5000;
const flush = () => new Promise((ok) => setTimeout(ok, 0));

function doc(track: string, direction: 1 | -1) {
  const b = createSeqViewSpec({ title: 'ranged' });
  const genome = b.data.source({ ref: 'genome', stream: 'hg38_genome', data_kind: 'letters' });
  const hg38 = b.data.complex({ ref: 'hg38', source_ref: genome, label: 'GRCh38' });
  for (const k of ['genes', 'variants']) {
    b.data.source({ ref: `${k}_src`, stream: `hg38_${k}`, data_kind: 'items', coordinate_system: 'hg38' });
    b.data.annotation({ ref: k, geometry: 'ranges', source_ref: `${k}_src` });
  }
  b.data.source({ ref: 'phylop_src', stream: 'hg38_phylop', data_kind: 'values', coordinate_system: 'hg38' });
  b.data.annotation({ ref: 'phylop', geometry: 'values', source_ref: 'phylop_src' });
  const g = b.view({ ref: 'v', label: 'v' }).section({ ref: 's', label: 'chr2' });
  g.block({ unit: hg38.unit('chr2'), range: [A, E], ...(direction < 0 ? { direction } : {}) });
  if (track === 'genes') g.track({ ref: 't', label: 'Genes', height: 3 }).ranges({ data: 'genes', labels: true });
  if (track === 'variants') g.track({ ref: 't', label: 'Variants' }).variants({ data: 'variants' });
  if (track === 'phylop') g.track({ ref: 't', label: 'phyloP' }).values({ data: 'phylop' });
  return b.getDocument();
}

async function open(track: string, direction: 1 | -1) {
  const s = hg38Streams({ letters: true });
  const by: Record<string, { provider: unknown }> = {
    hg38_genome: s.genome,
    hg38_genes: s.genes,
    hg38_variants: s.variants,
    hg38_phylop: s.phylop,
  };
  const loaders: Loaders = {
    load: (r) => Promise.reject(new Error(`no ${r.uri}`)),
    stream: (src) => (src.stream ? ((by[src.stream]?.provider as never) ?? null) : null),
  };
  const scheduler = manualScheduler();
  const ctx = createSeqViewer(
    { ...DefaultSeqViewerSpec(), registry: [DefaultRegistry] },
    { scheduler, loaders, text: fixtureMeasurer() },
  );
  const done = ctx.setDocument(doc(track, direction));
  for (const ref of ctx.viewports.sections.value) ctx.viewports.get(ref)?.setSize(1000, 400, 1);
  const r = await done;
  if (!r.ok) throw new Error(JSON.stringify(r.issues));
  const m = ctx.streams as Streams;
  const settle = async () => {
    for (let k = 0; k < 12; k++) {
      m.pass();
      scheduler.advance(16);
      await flush();
    }
    await m.settled();
  };
  const view = ctx.graph.assembled.value?.view;
  let layer = '';
  for (let k = 0; view && k < view.resolved.layers.count; k++)
    if (view.set.strings.get(view.resolved.tracks.ref[view.resolved.layers.track[k]]) === 't')
      layer = view.set.strings.get(view.resolved.layers.ref[k]);
  /** The layer's state over the positions the section shows now, every 1/40 of them. */
  const shown = () => {
    const ax = ctx.layout.axis('s');
    if (!ax) throw new Error('no axis');
    const a = ax.indexOf(0, ax.visLo[0] as never);
    const b = ax.indexOf(0, (ax.visHi[0] - 1) as never);
    const [lo, hi] = a < b ? [a, b] : [b, a];
    const out: number[] = [];
    for (let x = lo; x <= hi; x += Math.max(1, Math.floor((hi - lo) / 40))) out.push(m.state(layer, 0, x));
    return { lo, hi, states: out };
  };
  return { ctx, m, settle, shown };
}

describe('a ranged piece of a streamed chromosome (GM-01)', () => {
  for (const direction of [1, -1] as const)
    for (const track of ['genes', 'variants', 'phylop'])
      it(`${track}, direction ${direction}: every region it shows loads, whole, zoomed and panned`, async () => {
        const t = await open(track, direction);
        await t.settle();
        const whole = t.shown();
        // The whole piece is in view, and nothing outside it.
        expect([whole.lo, whole.hi]).toEqual([A, E - 1]);
        expect(new Set(whole.states)).toEqual(new Set([RegionState.loaded]));
        const vp = t.ctx.viewports.get('s');
        vp?.zoomTo(0, A + 70_000, A + 90_000);
        await t.settle();
        const zoomed = t.shown();
        expect(zoomed.hi - zoomed.lo).toBeLessThan(25_000);
        expect(new Set(zoomed.states)).toEqual(new Set([RegionState.loaded]));
        vp?.panBy(300, 0);
        await t.settle();
        const panned = t.shown();
        expect(panned.lo).not.toBe(zoomed.lo);
        expect(new Set(panned.states)).toEqual(new Set([RegionState.loaded]));
        expect(t.m.stats.current.inFlight).toBe(0);
        t.ctx.dispose();
      });
});
