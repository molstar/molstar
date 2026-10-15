/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The overlay paths' cases, for their heap guard and its quick twin (overlay-paths.heap.ts, overlay-paths.quick.ts):
// the vertical overview and scrollbar's overlay through a vertical scroll (`overview`), the keyboard cursor's ring
// through a pan (`ring`), a tree clade's hover (`clade`), a held select drag's pending band on an alignment's rows
// (`selectDrag`) and a dot plot's pending rectangle (`plotDrag`), each on the recording backend. Each call loads the
// viewer afresh (`vi.resetModules`) and returns the case's step (it gets a running index) and, for a drag, `held`: the
// check that the selection is still pending.

import { expect, vi } from 'vitest';

/** A case: its step, and the check to make after it has run. */
export interface PathCase {
  step: (i: number) => void;
  held?: () => void;
}

async function load() {
  vi.resetModules();
  const draw = await import('../src/draw/index.js');
  const viewer = await import('../src/index.js');
  const { recording, COPY_BACK } = await import('./composer-fixtures.js');
  const { driven } = await import('./input-fixtures.js');
  const { lit, twoBlocks } = await import('./overlay-fixtures.js');
  const { rowLoci } = await import('./overlay-docs.js');
  const families = await import('./families.js');
  const { createSeqViewSpec } = await import('@molstar/seq-core/builder');
  return { ...draw, ...viewer, recording, COPY_BACK, driven, lit, twoBlocks, rowLoci, ...families, createSeqViewSpec };
}

/**
 * A held drag's moves, as the engine sends them: one pointer record, its position set from typed arrays and dispatched
 * (`ctx.input.pointer`). The driven fixture's `t.pointer(…, { alt: true })` with positions computed per call was the
 * step's own allocation, about 50 and 90 bytes a step (an options object, the positions boxed into a call not inlined;
 * 2026-10-14, the closing pass's first run of this file), which the two drag cases measured.
 */
function moves(
  m: Awaited<ReturnType<typeof load>>,
  ctx: { input: { pointer(p: object): unknown }; layout: { sectionIndex(s: string): number } },
  section: string,
  xs: readonly number[],
  ys: readonly number[],
  modifiers: number,
) {
  const p = m.Inputs.pointer();
  p.type = m.PointerType.move;
  p.section = ctx.layout.sectionIndex(section);
  p.button = -1;
  p.buttons = 1;
  p.modifiers = modifiers;
  const [x, y] = [Float64Array.from(xs), Float64Array.from(ys)];
  const [nx, ny] = [x.length - 1, y.length - 1];
  return (i: number) => {
    p.x = x[i & nx];
    p.y = y[i & ny];
    ctx.input.pointer(p);
  };
}

/**
 * An alignment of 400 rows in the `msa` style (the vertical overview and the scrollbar), rows selected and highlighted
 * (the overview's row bars, LC-21), its overlay frames.
 */
export async function overview(): Promise<PathCase> {
  const m = await load();
  const t = await m.lit(m.familyWithTree(400, 30), { width: 1000, height: 400, dpr: 2, section: m.FAMILY_SECTION });
  t.ctx.ui.setStyle(t.section, 'msa');
  await t.ctx.graph.update();
  const set = t.ctx.graph.assembled.value?.view?.set;
  if (!set) throw new Error('no set');
  t.ctx.selection.set(m.rowLoci(set, 3, 0, 5), 'api');
  t.ctx.highlight.set('host', m.rowLoci(set, 250, 0, 5));
  const marks = new m.ChannelMarks({
    graph: t.ctx.graph,
    hover: t.ctx.hover,
    selection: t.ctx.selection,
    highlight: t.ctx.highlight,
    systems: t.ctx.systems,
    layout: t.ctx.layout,
  });
  const vnav = new m.VerticalNav(t.ctx, t.section, marks);
  const be = m.recording(m.COPY_BACK);
  vnav.frame(be, m.Dirty.all);
  expect(t.layout(t.section).vnavWidth).toBeGreaterThan(0);
  // The row bars: the selection's border and the highlight, in the overview.
  const colors = new Set(be.rectItems.color.subarray(0, be.rectItems.n));
  expect([colors.has(m.SkinToken.selectionBorder), colors.has(m.SkinToken.highlight)]).toEqual([true, true]);
  return {
    step: (i: number) => {
      t.vp?.panBy(0, (i & 1) === 0 ? 23 : -23, 'pointer');
      be.reset();
      vnav.frame(be, m.Dirty.overlay);
    },
  };
}

/** The keyboard cursor's ring through a pan. */
export async function ring(): Promise<PathCase> {
  const m = await load();
  const t = await m.driven(m.twoBlocks(), { renderers: m.KernelRenderers as never, size: { v: [800, 300] } });
  const section = t.ctx.viewports.sections.value[0];
  const vp = t.ctx.viewports.get(section);
  vp?.setView(30, 0, 'api');
  t.focus(section);
  t.key(section, 'ArrowRight');
  expect(t.ctx.ui.keyCursor.value).not.toBeNull();
  const c = new m.SectionComposer(t.ctx, section);
  const be = m.recording(m.COPY_BACK);
  c.frame(be, m.Dirty.all);
  return {
    step: (i) => {
      vp?.panBy((i & 1) === 0 ? 3 : -3, 0, 'pointer');
      be.reset();
      c.frame(be, m.Dirty.all);
    },
  };
}

/** A tree clade's hover, on and off. */
export async function clade(): Promise<PathCase> {
  const m = await load();
  const t = await m.lit(m.familyWithTree(40, 10), { width: 1000, height: 400, dpr: 2, section: m.FAMILY_SECTION });
  t.ctx.ui.toggleClade('family', 2);
  await t.ctx.graph.update();
  const l = t.layout(t.section);
  const top = l.scrolledTop + l.entryTop[l.rows.entry];
  // The wedge's row (a clade's hover) and the last row (a leaf's), alternately.
  let wedge = -1;
  for (let d = 0; d < l.rows.count && wedge < 0; d++) if (l.rows.node?.[d] === 2) wedge = d;
  expect(wedge).toBeGreaterThanOrEqual(0);
  const ys = [wedge, l.rows.count - 1].map((d) => top + (d + 0.5) * l.rows.height);
  const x = l.headerWidth - 4;
  t.hover(x, ys[0]);
  expect(t.ctx.hover.live.current.kind).toBe(m.HoverKind.header);
  const be = m.recording(m.COPY_BACK);
  t.frame(be, m.Dirty.all);
  return {
    step: (i) => {
      t.hover(x, ys[i & 1]);
      t.frame(be, m.Dirty.overlay);
    },
  };
}

/** A held select drag's pending band on an alignment's rows. */
export async function selectDrag(): Promise<PathCase> {
  const m = await load();
  const t = await m.driven(m.familyWithTree(60, 20), {
    renderers: m.KernelRenderers as never,
    size: { [m.FAMILY_SECTION]: [1000, 400] },
  });
  const s = m.FAMILY_SECTION;
  t.ctx.viewports.get(s)?.setView(12, 0, 'api');
  const l = t.layout(s);
  const y = (d: number) => l.scrolledTop + l.entryTop[l.rows.entry] + (d + 0.5) * l.rows.height;
  t.pointer('down', s, t.xOf(s, 0, 5), y(2), { alt: true });
  t.pointer('move', s, t.xOf(s, 0, 9), y(4), { alt: true });
  expect(t.ctx.selection.pending.current.active).toBe(true);
  const c = new m.SectionComposer(t.ctx, s);
  const be = m.recording(m.COPY_BACK);
  c.frame(be, m.Dirty.all);
  const xs = [0, 1, 2, 3, 4, 5, 6, 7].map((k) => t.xOf(s, 0, 9 + k));
  const move = moves(m, t.ctx, s, xs, [4, 5, 6, 7].map(y), m.Modifier.alt);
  return {
    step: (i) => {
      move(i);
      be.reset();
      c.frame(be, m.Dirty.overlay);
    },
    held: () => expect(t.ctx.selection.pending.current.active).toBe(true),
  };
}

/** A held drag's pending rectangle on a dot plot. */
export async function plotDrag(): Promise<PathCase> {
  const m = await load();
  const b = m.createSeqViewSpec();
  const a = b.data.sequence({ ref: 'a', value: 'MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRV' });
  const value = Array(40 * 40).fill(1);
  const grid = b.data.annotation({
    ref: 'g',
    geometry: 'grid',
    value: { target: ['a'], b_target: ['a'], value: [value] },
  });
  const sec = b.view({ ref: 'v', label: 'V' }).section({ ref: 'v.s' });
  sec.block({ sequence: a });
  sec.track({ ref: 't-plot', label: 'Plot' }).dot_plot({ data: grid });
  const t = await m.driven(b.getDocument(), { renderers: m.KernelRenderers as never, size: { 'v.s': [800, 600] } });
  const l = t.layout('v.s');
  const v = t.ctx.layout.axis('v.s');
  if (!v) throw new Error('no axis');
  const e = l.entryCount - 1;
  const top = l.scrolledTop + l.entryTop[e] - v.scrollY;
  const yOf = (k: number) => top + ((t.xOf('v.s', 0, k) - v.left) / (v.right - v.left)) * l.entryHeight[e];
  t.pointer('down', 'v.s', t.xOf('v.s', 0, 5), yOf(5));
  t.pointer('move', 'v.s', t.xOf('v.s', 0, 12), yOf(14));
  expect(t.ctx.selection.pending.current.active).toBe(true);
  const c = new m.SectionComposer(t.ctx, 'v.s');
  const be = m.recording(m.COPY_BACK);
  c.frame(be, m.Dirty.all);
  const xs = [0, 1, 2, 3, 4, 5, 6, 7].map((k) => t.xOf('v.s', 0, 12 + k));
  const move = moves(m, t.ctx, 'v.s', xs, [14, 15, 16, 17].map(yOf), 0);
  return {
    step: (i) => {
      move(i);
      be.reset();
      c.frame(be, m.Dirty.overlay);
    },
    held: () => expect(t.ctx.selection.pending.current.active).toBe(true),
  };
}
