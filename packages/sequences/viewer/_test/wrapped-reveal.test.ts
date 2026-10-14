/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Wrapped text scrolls to a focus (seqspec:interaction#focus; the `focusViewport` behavior's `reveal`): the lines of
// the focus come into view, as the tracks' zoom brings the focus into view, over the same transition; nothing moves
// where they already are shown (the pad is where a scroll puts them, not a test of being shown). A new section opens
// on its focus, and a focus before the text has lines waits for them (seqspec:interaction#viewports). Hidden ends,
// reversed pieces and a focus over several blocks reveal every line they light.

import { createSeqViewSpec } from '@molstar/seq-core/builder';
import type { SeqViewSpecDoc } from '@molstar/seq-core/format';
import { describe, expect, it } from 'vitest';
import { DefaultPresets, KernelRenderers, type LayerRenderer } from '../src/index.js';
import { opened } from './layout-fixtures.js';

const S = 'v.main';
const SEQ = 'MSEEKLKELAEKAIRLAKEGNPEEVLRALQ';
const residues = (n: number) => SEQ.repeat(Math.ceil(n / 30)).slice(0, n);

/** `n` positions as text; a block's `range` and `direction`, a section `focus`, or a second sequence as a block. */
function text(
  n: number,
  o: { range?: [number, number]; direction?: -1; focus?: [number, number]; second?: number } = {},
): SeqViewSpecDoc {
  const b = createSeqViewSpec();
  b.data.sequence({ ref: 'a', alphabet: 'protein', value: residues(n) });
  if (o.second) b.data.sequence({ ref: 'c', alphabet: 'protein', value: residues(o.second) });
  const v = b.view({ ref: 'v', label: 'V', layout: 'text' });
  if (o.range || o.direction || o.second)
    v.block({ sequence: 'a', ...(o.range ? { range: o.range } : {}), ...(o.direction ? { direction: -1 } : {}) });
  if (o.second) v.block({ sequence: 'c' });
  if (o.focus) v.focus({ target: 'a', start: o.focus[0], stop: o.focus[1] });
  v.track({ ref: 'v.letters' }).letters();
  return b.getDocument();
}

const BUILT_IN = {
  presets: false as const,
  renderers: [...KernelRenderers] as LayerRenderer[],
  registry: [{ presets: [...DefaultPresets] }],
};

/** 2,000 positions (or `doc`) in a 400 × 120 px strip (`width` 0: no lines yet): about 50 lines, a few visible. */
async function open(doc: SeqViewSpecDoc = text(2000), width = 400) {
  const t = await opened(doc, {
    ...BUILT_IN,
    size: { [S]: [width, 120] },
    render: { wrappedGroupLabels: 'none', wrappedGutter: 0 },
  });
  t.layout(S);
  return t;
}
type T = Awaited<ReturnType<typeof open>>;

const scrollY = (t: T) => t.ctx.viewports.get(S)?.live.current.scrollY ?? -1;
/** A focus from loci on sequence `ref`, then the animation run out (`settle` false: not yet). */
const focus = (t: T, start: number, stop: number, o: { transitionMs?: number; settle?: false; ref?: string } = {}) => {
  const loci = t.ctx.refs.toLoci(o.ref ?? 'a', { start, stop });
  expect(loci).not.toBeNull();
  const ok = t.ctx.focus.fromLoci(loci as never, {
    section: S,
    ...(o.transitionMs !== undefined ? { transitionMs: o.transitionMs } : {}),
  });
  if (o.settle !== false) t.scheduler.advance(1000);
  return ok;
};
/** The view's height under what is pinned, px. */
const viewOf = (t: T) => t.layout(S).height - t.layout(S).scrolledTop;
/** The line (top, bottom: its last entry's end) holding block `b`'s shown position k. */
function lineOf(t: T, k: number, b = 0) {
  const l = t.layout(S);
  for (let n = 0; n < l.lines.count; n++) {
    if (l.lines.block[n] !== b || k < l.lines.start[n] || k >= l.lines.stop[n]) continue;
    let bottom = l.lines.top[n];
    for (let e = 0; e < l.entryCount; e++)
      if (l.entryRow[e] === n) bottom = Math.max(bottom, l.entryTop[e] + l.entryHeight[e]);
    return { n, top: l.lines.top[n], bottom };
  }
  throw new Error(`no line holds ${k} of block ${b}`);
}
/** Whether the line holding block b's shown k is wholly in view. */
const shown = (t: T, k: number, b = 0) => {
  const line = lineOf(t, k, b);
  return line.top >= scrollY(t) && line.bottom <= scrollY(t) + viewOf(t) + 0.5;
};

describe('wrapped text scrolls to a focus', () => {
  it('a focus on a far residue brings its line into view', async () => {
    const t = await open();
    expect(t.ctx.layout.isWrapped(S)).toBe(true);
    expect(scrollY(t)).toBe(0);
    expect(focus(t, 1500, 1510)).toBe(true);
    expect(scrollY(t)).toBeGreaterThan(0);
    expect(shown(t, 1500)).toBe(true);
  });

  it("is animated over the focus's transition, as the tracks' zoom is; a transition of 0 and `animate: false` jump", async () => {
    const t = await open();
    focus(t, 1500, 1510, { settle: false });
    const vp = t.ctx.viewports.get(S);
    expect(vp?.live.current.animating).toBe(true);
    t.scheduler.advance(50);
    const mid = scrollY(t);
    t.scheduler.advance(1000);
    const end = scrollY(t);
    expect(mid).toBeGreaterThan(0);
    expect(mid).toBeLessThan(end);
    expect(vp?.live.current.animating).toBe(false);
    focus(t, 0, 5, { transitionMs: 0, settle: false });
    expect(scrollY(t)).toBe(0);
    t.ctx.behaviors.get('focusViewport')?.setParams({ animate: false });
    focus(t, 1500, 1510, { settle: false });
    expect(scrollY(t)).toBe(end);
  });

  it('no scroll where the lines already are shown, even within the pad of an edge', async () => {
    const t = await open();
    expect(focus(t, 3, 8)).toBe(true);
    expect(scrollY(t)).toBe(0);
    focus(t, 1500, 1510);
    const there = scrollY(t);
    // Again, and a residue on the same line: nothing moves.
    focus(t, 1500, 1510);
    focus(t, 1502, 1504);
    expect(scrollY(t)).toBe(there);
    // A line shown 2 px above the bottom edge (inside `revealPad`, 12 px): shown, so nothing moves.
    const line = lineOf(t, 400);
    t.ctx.viewports.get(S)?.scrollTo(0, line.bottom - viewOf(t) + 2);
    const at = scrollY(t);
    expect(shown(t, 400)).toBe(true);
    focus(t, 400, 401);
    expect(scrollY(t)).toBe(at);
  });

  it('scrolls back up to a focus above the view, `revealPad` px clear of the top (0 at the start)', async () => {
    const t = await open();
    focus(t, 1500, 1510);
    expect(scrollY(t)).toBeGreaterThan(0);
    focus(t, 800, 805);
    expect(scrollY(t)).toBe(lineOf(t, 800).top - 12);
    focus(t, 0, 5);
    expect(scrollY(t)).toBe(0);
  });

  it('is a behavior param: `reveal: false` leaves the text where it is', async () => {
    const t = await open();
    t.ctx.behaviors.get('focusViewport')?.setParams({ reveal: false });
    focus(t, 1500, 1510);
    expect(scrollY(t)).toBe(0);
  });

  it("a new section opens on its focus: the document's focus node's line is in view", async () => {
    const t = await open(text(2000, { focus: [1500, 1510] }));
    expect(scrollY(t)).toBeGreaterThan(0);
    expect(shown(t, 1500)).toBe(true);
  });

  it('a focus before the text has lines (no width yet) waits for them, then scrolls to it', async () => {
    const t = await open(text(2000), 0);
    expect(t.layout(S).lines.count).toBe(0);
    expect(focus(t, 1500, 1510)).toBe(true);
    expect(scrollY(t)).toBe(0);
    t.ctx.viewports.get(S)?.setSize(400, 120, 1);
    t.layout(S);
    t.scheduler.advance(1000);
    expect(t.layout(S).lines.count).toBeGreaterThan(0);
    expect(shown(t, 1500)).toBe(true);
  });

  it('a range whose end is past the piece shows its last shown line too; one whose start is before it, its first', async () => {
    // The piece is [100, 1505): shown k is index − 100.
    const t = await open(text(2000, { range: [100, 1505] }));
    // [1490, 1600): shown from 1390 to the piece's end, 1404, across two lines (40 positions a line).
    expect(lineOf(t, 1390).n).not.toBe(lineOf(t, 1404).n);
    focus(t, 1490, 1600);
    expect(shown(t, 1390)).toBe(true);
    expect(shown(t, 1404)).toBe(true);
    // [50, 130): its start is before the piece; shown from 0, the first line.
    focus(t, 50, 130);
    expect(scrollY(t)).toBe(0);
  });

  it("a reversed piece reveals the range's lines in shown order (its stop's line comes first)", async () => {
    const t = await open(text(2000, { direction: -1 }));
    // Shown k is 1999 − index: [1500, 1560) is shown 440..499, lines 11 and 12; from the end, the first comes in at
    // the top.
    const vp = t.ctx.viewports.get(S);
    vp?.scrollTo(0, 1e9);
    expect(lineOf(t, 440).n).not.toBe(lineOf(t, 499).n);
    focus(t, 1500, 1560);
    expect(shown(t, 440)).toBe(true);
    expect(shown(t, 499)).toBe(true);
  });

  it("a focus over several blocks (its `blocks`) reveals every block's lines", async () => {
    // Two blocks of 1,000: the end of the first and the start of the second, on adjacent lines.
    const t = await open(text(1000, { second: 1000 }));
    const last = lineOf(t, 985, 0);
    const first = lineOf(t, 5, 1);
    expect(first.n).toBe(last.n + 1);
    t.ctx.focus.set(
      S,
      { block: 0, start: 980, stop: 990, blocks: Int32Array.of(0, 980, 990, 1, 0, 10), loci: null },
      { transitionMs: 0 },
    );
    expect(shown(t, 985, 0)).toBe(true);
    expect(shown(t, 5, 1)).toBe(true);
  });
});
