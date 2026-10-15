/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The performance pass's quick cases (.v6/seq/plans/lanes-input-perf.md §4.2; perf-pass.quick.ts): the host options'
// hot paths and this PR's new ones, each a section on the real composer and dispatcher with a step:
// - `unobservedLinear`: letters' cells with `unobservedGray`, zoomed in on a linear section over a unit with many
//   `observed` ranges, panned a step each (the per-unit mask, `LettersState.observed`);
// - `unobservedText`: wrapped text with `unobservedGray` (tinted by `wrappedColors`, or plain), letters over underline
//   and strip lanes, scrolled a step each;
// - `centeredText`: a centered reading column (`wrappedMaxWidth`, `wrappedAlign: 'center'`): a select drag's moves
//   with an overlay frame each, or the hit test alone;
// - `chips`: short ligand units as chips in the text: the pointer over them with an overlay frame each, or the rows
//   scrolled;
// - `textTracks`: a track drawn in full in the text (`wrappedTracks: 'full'`) per kind (values line, area and bars,
//   points, ranges), scrolled a step each.
// Each call loads the viewer afresh (`vi.resetModules`), as the guards do.

import type { SeqViewSpecDoc } from '@molstar/seq-core/format';
import { vi } from 'vitest';
import type { AllocOptions } from '../../core/_test/model/alloc.js';

const S = 'v.main';

/** A full-frame step's quick counts (a frame redraws every line or cell: about a millisecond). */
export const FRAME: AllocOptions = { warmSteps: 1_000, steps: 250, windows: 3 };

/** A unit of `n` residues observed over many ranges: 40 of them, uneven, with one-position islands. */
function observed(n: number): [number, number][] {
  const out: [number, number][] = [];
  for (let k = 0, at = 3; k < 40 && at < n - 2; k++) {
    const len = k % 5 === 4 ? 1 : 7 + ((k * 13) % 29);
    out.push([at, Math.min(n, at + len)]);
    at += len + 2 + ((k * 7) % 11);
  }
  return out;
}

/** The viewer over `doc` on the real composer (recording backend), `text` on the wrapped preset. */
async function bench(
  doc: SeqViewSpecDoc,
  o: { size: [number, number]; text?: boolean; render?: Record<string, unknown>; copyBack?: boolean },
) {
  const { driven } = await import('./input-fixtures.js');
  const { recording, COPY_BACK } = await import('./composer-fixtures.js');
  const { DefaultPresets, Hit, Inputs, KernelRenderers, PointerType } = await import('../src/index.js');
  const { Dirty, SectionComposer } = await import('../src/draw/index.js');
  const text = DefaultPresets.find((p) => p.name === 'wrapped');
  const t = await driven(doc, {
    size: { [S]: o.size },
    registry: o.text ? [{ presets: [{ ...(text as NonNullable<typeof text>), name: 'lab:text' }] }] : [],
    renderers: KernelRenderers as never,
    labels: false,
  });
  if (o.text) t.ctx.ui.setStyle(S, 'lab:text');
  if (o.render) t.ctx.setRenderProps(o.render as never);
  await t.ctx.graph.update();
  // The frames alone: the hover label's strings are one allocation per hovered element by design (IN2's guard).
  t.ctx.behaviors.get('hoverLabel')?.remove();
  const composer = new SectionComposer(t.ctx, S);
  const be = o.copyBack ? recording(COPY_BACK) : recording();
  composer.frame(be, Dirty.all);
  const l = t.layout(S);
  const h = Hit.create();
  /**
   * Up to 32 points on lines 1.. where `keep(hit)` holds, 7 px apart, 2.5 to 32.5 px under each line's top. The `hit`
   * case calls `hitTestAt` with the point in the record, as the dispatcher does: `hitTest(section, x, y)` from this
   * closure boxed a fractional coordinate, 16 bytes a call, where V8 did not inline it (the caller's cost).
   */
  const points = (keep: (hit: typeof h) => boolean) => {
    const pts: number[] = [];
    for (let line = 1; line < l.lines.count && pts.length < 64; line++)
      for (let dy = 2.5; dy < 40 && pts.length < 64; dy += 10)
        for (let x = l.dataLeft + 2; x < l.dataRight - 2 && pts.length < 64; x += 7) {
          const y = l.scrolledTop + l.lines.top[line] + dy;
          if (t.ctx.layout.hitTest(S, x, y, h) && keep(h)) pts.push(x, y);
        }
    return Float64Array.from(pts);
  };
  const p = Inputs.pointer();
  p.section = t.ctx.layout.sectionIndex(S);
  p.type = PointerType.move;
  p.button = -1;
  p.buttons = 0;
  return { t, composer, be, Dirty, vp: t.vp(S), l, h, p, PointerType, points };
}

/** Letters' cells with `unobservedGray` on a linear section zoomed in, panned left and right a step each. */
export async function unobservedLinear() {
  vi.resetModules();
  const { createSeqViewSpec } = await import('@molstar/seq-core/builder');
  const n = 2000;
  const b = createSeqViewSpec();
  b.data.sequence({ ref: 'chain', value: 'MKVLATGHEK'.repeat(n / 10) });
  b.data.unit({ ref: 'A', sequence: 'chain', observed: observed(n) } as never);
  const v = b.view({ ref: 'v', label: 'V' });
  v.block({ unit: 'A' } as never);
  v.track({ ref: 'v.letters', label: 'Letters' }).letters({
    color: { theme: 'letter', params: { colors: 'Hydrophobicity' } },
  } as never);
  const x = await bench(b.getDocument(), { size: [1000, 200], render: { unobservedGray: true } });
  x.vp?.setView(14, 300, 'api');
  await x.t.ctx.graph.update();
  x.composer.frame(x.be, x.Dirty.all);
  const step = (i: number) => {
    x.be.reset();
    x.vp?.panBy(i & 1 ? -11 : 11, 0, 'pointer');
    x.composer.frame(x.be, x.Dirty.tracks);
  };
  return { step, be: x.be };
}

/** Wrapped text with `unobservedGray` (`tinted`: and `wrappedColors`) over underline and strip lanes, scrolled. */
export async function unobservedText(tinted: boolean) {
  vi.resetModules();
  const { createSeqViewSpec } = await import('@molstar/seq-core/builder');
  const { BANDS, BAND_SCORES } = await import('./prepare-fixtures.js');
  const n = 2000;
  const b = createSeqViewSpec();
  const a = b.data.sequence({ ref: 'chain', value: 'MKVLATGHEK'.repeat(n / 10) });
  b.data.unit({ ref: 'A', sequence: 'chain', observed: observed(n) } as never);
  const f = b.data.annotation({
    ref: 'f',
    geometry: 'ranges',
    target: a,
    value: Array.from({ length: 60 }, (_, i) => ({ start: i * 30, stop: i * 30 + 45, label: `f${i}` })),
  });
  const value = Array.from({ length: n }, (_, i) => (i % 23 === 7 ? null : BAND_SCORES[Math.floor(i / 7) % 4]));
  const d = b.data.annotation({ ref: 'd', geometry: 'values', target: a, value: { value } } as never);
  const v = b.view({ ref: 'v', label: 'V', layout: 'text', height: 'box' });
  v.block({ unit: 'A' } as never);
  v.track({ ref: 'v.letters' }).letters({ color: { theme: 'letter', params: { colors: 'Hydrophobicity' } } } as never);
  v.track({ ref: 'v.features' }).ranges({ data: f, display: 'underline' } as never);
  v.track({ ref: 'v.plddt' }).values({
    data: d,
    domain: [0, 100],
    color: { theme: 'field_palette', params: { field: 'value', palette: BANDS } },
  } as never);
  const x = await bench(b.getDocument(), {
    size: [500, 1200],
    text: true,
    render: { unobservedGray: true, wrappedColors: tinted },
  });
  const step = (i: number) => {
    x.be.reset();
    x.vp?.panBy(0, i & 1 ? -11 : 11, 'pointer');
    x.composer.frame(x.be, x.Dirty.tracks);
  };
  return { step, be: x.be };
}

/** A centered reading column: `select` moves a press's drag over its lines, an overlay frame a move; `hit` tests. */
export async function centeredText(mode: 'select' | 'hit') {
  vi.resetModules();
  const { createSeqViewSpec } = await import('@molstar/seq-core/builder');
  const b = createSeqViewSpec();
  const a = b.data.sequence({ ref: 'a', value: 'ACDEFGHIKL'.repeat(200) });
  const f = b.data.annotation({
    ref: 'f',
    geometry: 'ranges',
    target: a,
    value: Array.from({ length: 60 }, (_, i) => ({ start: i * 30, stop: i * 30 + 45, label: `f${i}` })),
  });
  const v = b.view({ ref: 'v', label: 'V', layout: 'text', height: 'box' });
  v.block({ sequence: a });
  v.track({ ref: 'v.letters' }).letters();
  v.track({ ref: 'v.features' }).ranges({ data: f, display: 'underline' } as never);
  const x = await bench(b.getDocument(), {
    size: [1200, 1200],
    text: true,
    render: { wrappedMaxWidth: 420, wrappedAlign: 'center' },
    copyBack: true,
  });
  const pts = x.points((h) => h.index >= 0);
  const count = pts.length / 2;
  if (count < 16) {
    const y = x.l.scrolledTop + x.l.lines.top[1] + 10;
    const xm = (x.l.dataLeft + x.l.dataRight) / 2;
    const ok = x.t.ctx.layout.hitTest(S, xm, y, x.h);
    throw new Error(
      `only ${count} points: lines ${x.l.lines.count}, data ${x.l.dataLeft}..${x.l.dataRight}, top1 ${y}, ` +
        `hit ${ok} kind ${x.h.kind} block ${x.h.block} index ${x.h.index} track ${x.h.track}`,
    );
  }
  // Centered: the column's lines start well inside the section (420 px of 1,200).
  if (x.l.dataLeft < 300) throw new Error(`the column is not centered (it starts at ${x.l.dataLeft})`);
  const { t, p, PointerType } = x;
  if (mode === 'select') {
    p.type = PointerType.down;
    p.button = 0;
    p.buttons = 1;
    p.x = pts[0];
    p.y = pts[1];
    t.ctx.input.pointer(p);
    p.type = PointerType.move;
    p.button = -1;
    p.x += 30.5;
    t.ctx.input.pointer(p);
  }
  const pd = t.ctx.selection.pending.current;
  let hits = 0;
  let ends = 0;
  let last = -1;
  const step = (i: number) => {
    const k = i % count;
    if (mode === 'hit') {
      x.h.x = pts[2 * k];
      x.h.y = pts[2 * k + 1];
      if (t.ctx.layout.hitTestAt(p.section, x.h)) hits++;
      return;
    }
    x.be.reset();
    p.x = pts[2 * k];
    p.y = pts[2 * k + 1];
    t.ctx.input.pointer(p);
    if (pd.endShown !== last) ends++;
    last = pd.endShown;
    x.composer.frame(x.be, x.Dirty.overlay);
  };
  const state = () => ({ owner: t.ctx.ui.drag.value?.owner, active: pd.active, ends, hits });
  return { step, state };
}

/** Short ligand units as chips in the text: `hover` the pointer over them, an overlay frame each; `scroll` the rows. */
export async function chips(mode: 'hover' | 'scroll') {
  vi.resetModules();
  const { createSeqViewSpec } = await import('@molstar/seq-core/builder');
  const LIGANDS = ['HEM', 'ERY', 'NI', 'ZN', 'ATP', 'NAG', 'SO4', 'GOL', 'MG', 'CA', 'NA', 'CL'];
  const b = createSeqViewSpec();
  const a = b.data.sequence({ ref: 'a', alphabet: 'protein', value: 'MSEEKLKELAEKAIRLAKEGNPEEVLRALQ'.repeat(20) });
  const v = b.view({ ref: 'v', label: 'V', layout: 'text', height: 'box' });
  v.block({ sequence: a });
  for (let r = 0; r < 4; r++)
    LIGANDS.forEach((id, k) => {
      b.data.ligand({ ref: `l${r}.${k}`, value: ['C', 'N', 'O'], chem_comp_ids: [id] });
      v.block({ ligand: `l${r}.${k}` });
    });
  const r = b.data.annotation({
    ref: 'r',
    geometry: 'ranges',
    target: 'a',
    value: { start: [3, 200], stop: [20, 420] },
  });
  v.track({ ref: 'v.letters', label: 'Letters' }).letters();
  v.track({ ref: 'v.r', label: 'Region' }).ranges({ data: r, display: 'underline' } as never);
  const x = await bench(b.getDocument(), {
    size: [400, mode === 'scroll' ? 300 : 2400],
    text: true,
    copyBack: mode === 'hover',
  });
  const pts = x.points((h) => h.block >= 1 && h.index >= 0);
  const count = pts.length / 2;
  if (mode === 'hover' && count < 8) throw new Error(`only ${count} points on chips`);
  const { t, p } = x;
  let hovered = 0;
  const step = (i: number) => {
    x.be.reset();
    if (mode === 'scroll') {
      x.vp?.panBy(0, i & 1 ? -11 : 11, 'pointer');
      x.composer.frame(x.be, x.Dirty.tracks);
      return;
    }
    const k = i % count;
    p.x = pts[2 * k];
    p.y = pts[2 * k + 1];
    t.ctx.input.pointer(p);
    if (t.ctx.hover.live.current.index >= 0) hovered++;
    x.composer.frame(x.be, x.Dirty.overlay);
  };
  return { step, hovered: () => hovered, be: x.be };
}

/** The kinds `textTracks` draws in full. */
export const TEXT_TRACKS = ['line', 'area', 'bars', 'points', 'ranges'] as const;

/** A track of `kind` drawn in full in the text (`wrappedTracks: 'full'`), beside the letters, scrolled a step each. */
export async function textTracks(kind: (typeof TEXT_TRACKS)[number]) {
  vi.resetModules();
  const { createSeqViewSpec } = await import('@molstar/seq-core/builder');
  const n = 1200;
  const b = createSeqViewSpec();
  const a = b.data.sequence({ ref: 'a', alphabet: 'protein', value: 'MSEEKLKELAEKAIRLAKEGNPEEVLRALQ'.repeat(n / 30) });
  const v = b.view({ ref: 'v', label: 'V', layout: 'text', height: 'box' });
  v.block({ sequence: a });
  v.track({ ref: 'v.letters', label: 'Letters' }).letters();
  if (kind === 'points' || kind === 'ranges') {
    const value =
      kind === 'points'
        ? { start: Array.from({ length: 120 }, (_, i) => i * 10 + 3) }
        : {
            start: Array.from({ length: 40 }, (_, i) => i * 30),
            stop: Array.from({ length: 40 }, (_, i) => i * 30 + 18),
            label: Array.from({ length: 40 }, (_, i) => `d${i}`),
          };
    const d = b.data.annotation({ ref: 'd', geometry: 'ranges', target: a, value } as never);
    const tr = v.track({ ref: 'v.t', label: 'T', height: kind === 'ranges' ? 1 : 1 });
    if (kind === 'points') tr.points({ data: d });
    else tr.ranges({ data: d, labels: true } as never);
  } else {
    const value = Array.from({ length: n }, (_, i) => (i % 97 === 5 ? null : Math.round(50 + 40 * Math.sin(i / 6))));
    const d = b.data.annotation({ ref: 'd', geometry: 'values', target: a, value: { value } } as never);
    v.track({ ref: 'v.t', label: 'T', height: 2 }).values({ data: d, encoding: kind, domain: [0, 100] } as never);
  }
  const x = await bench(b.getDocument(), {
    size: [600, 800],
    text: true,
    render: { wrappedTracks: 'full', wrappedTrackAxes: true, wrappedTrackLabels: true },
  });
  const step = (i: number) => {
    x.be.reset();
    x.vp?.panBy(0, i & 1 ? -11 : 11, 'pointer');
    x.composer.frame(x.be, x.Dirty.tracks);
  };
  return { step, be: x.be, entries: x.l.entryCount };
}
