/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A select drag in a wrapped section allocates nothing per pointer event or per frame after warm-up
// (seqspec:interaction#selection; seqspec:viewer#performance): on the text preset, through the real dispatcher and
// behaviors, a press on the letters, then moves across four lines, past their ends and over their labels bands, the
// overlay drawing the pending range each step, and over the features' underlines by lane (each hit an item); and the
// scroll past the rows' edges, a frame a step, the range following.
// The positive control: the same moves allocating an object each (`leaking`), warmed and measured alike, are caught.
// Each case loads the viewer afresh (`vi.resetModules`): cases in one file share V8's feedback.

import { describe, expect, it, vi } from 'vitest';
import {
  controlCaught,
  describeScaling,
  FRAME_SETTLE_MS,
  FRAME_WARM,
  FRAME_WARM_MS,
  heapScaling,
  isFlat,
  leaking,
  settleFrame,
} from './heap.js';

/** A deep guard's steps: 500 catch about 43 bytes a step, these about 22. */
const STEPS = 1_000;

/** The drag's viewer and step, and the file's warm-up and measurement of a step (the guards' and the control's). */
async function setup(mode: 'lines' | 'scroll') {
  vi.resetModules();
  const { driven } = await import('./input-fixtures.js');
  const { recording, COPY_BACK } = await import('./composer-fixtures.js');
  const { createSeqViewSpec } = await import('@molstar/seq-core/builder');
  const { DefaultPresets, Hit, HoverKind, Inputs, KernelRenderers, PointerType } = await import('../src/index.js');
  const { runX } = await import('../src/runs.js');
  const { Dirty, SectionComposer } = await import('../src/draw/index.js');
  const S = 'v.main';
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
  v.track({ ref: 'v.features' }).ranges({ data: f });
  const text = DefaultPresets.find((p) => p.name === 'wrapped');
  const t = await driven(b.getDocument(), {
    size: { [S]: [500, mode === 'scroll' ? 300 : 1200] },
    registry: [{ presets: [{ ...(text as NonNullable<typeof text>), name: 'lab:text' }] }],
    renderers: KernelRenderers as never,
    labels: false,
  });
  t.ctx.ui.setStyle(S, 'lab:text');
  await t.ctx.graph.update();
  // The frames alone: the hover label's strings are one allocation per hovered element by design (IN2's guard).
  t.ctx.behaviors.get('hoverLabel')?.remove();
  // The scroll's frame, run by the step: the manual scheduler's frame map is the test's, not the viewer's.
  let frame: ((t: number) => void) | null = null;
  const port = t.scheduler as unknown as {
    requestFrame(cb: (t: number) => void): number;
    cancelFrame(id: number): void;
  };
  port.requestFrame = (cb) => {
    frame = cb;
    return 1;
  };
  port.cancelFrame = () => {
    frame = null;
  };
  const l = t.layout(S);
  const ax = t.ctx.layout.axis(S);
  const run0 = ax?.runs[0];
  if (!ax || !run0 || l.lines.count < 6) throw new Error('no lines');
  const composer = new SectionComposer(t.ctx, S);
  const be = recording(COPY_BACK);
  composer.frame(be, Dirty.all);
  // Points over lines 1 to 4: on letters, past a line's end, on its labels band, between groups, and on the features'
  // underlines, lane 0 and lane 1 (each an item's hit by lane, `hitLineLanes`).
  const features = t.trackIndex('v.features');
  const band = (line: number) => {
    for (let e = 0; e < l.entryCount; e++)
      if (l.entryTrack[e] === features && l.entryRow[e] === line) return l.entryTop[e] + l.entryHeight[e];
    throw new Error(`no features band on line ${line}`);
  };
  const pts = new Float64Array(64);
  for (let i = 0; i < 32; i++) {
    const line = 1 + (i % 4);
    const k = l.lines.start[line] + ((i * 7) % (l.lines.stop[line] - l.lines.start[line]));
    pts[2 * i] = i % 8 === 3 ? l.dataRight - 1.5 : l.dataLeft + runX(run0, l.lines.start[line], k, ax.gap) + 3.25;
    const lane = i % 8 === 6 ? 0 : i % 8 === 7 ? 1 : -1;
    pts[2 * i + 1] =
      lane >= 0
        ? l.scrolledTop + band(line) - 5 - 5 * lane + 1.25
        : l.scrolledTop + l.lines.top[line] + (i % 8 === 5 ? 2.5 : 22.5);
  }
  // The lanes' points hit items (the hover path the guard covers).
  const h = Hit.create();
  let items = 0;
  for (let i = 0; i < 32; i++)
    if (t.ctx.layout.hitTest(S, pts[2 * i], pts[2 * i + 1], h) && h.kind === HoverKind.item) items++;
  if (items < 4) throw new Error(`the lanes' points hit ${items} items`);
  const p = Inputs.pointer();
  p.section = t.ctx.layout.sectionIndex(S);
  p.type = PointerType.down;
  p.button = 0;
  p.buttons = 1;
  p.x = pts[0];
  p.y = l.scrolledTop + l.lines.top[1] + 22.5;
  t.ctx.input.pointer(p);
  p.type = PointerType.move;
  p.button = -1;
  p.x += 30.5;
  t.ctx.input.pointer(p);
  const pd = t.ctx.selection.pending.current;
  const vp = t.vp(S);
  let ends = 0;
  let last = -1;
  let scrolls = 0;
  /** The last scroll, in a Float64Array: a double in a closure's variable is boxed on each write. */
  const lastY = new Float64Array(1);
  const step = (i: number) => {
    if (mode === 'lines') {
      p.x = pts[2 * (i & 31)];
      p.y = pts[2 * (i & 31) + 1];
      t.ctx.input.pointer(p);
    } else {
      // 40 steps below the rows, 40 above: the rows scroll down and back, a frame a step.
      if (i % 40 === 0) {
        p.y = (i / 40) % 2 === 0 ? l.height + 25.5 : -25.5;
        t.ctx.input.pointer(p);
      }
      (frame as ((t: number) => void) | null)?.(16 * i);
    }
    if (pd.endShown !== last) ends++;
    last = pd.endShown;
    if (vp.live.current.scrollY !== lastY[0]) scrolls++;
    lastY[0] = vp.live.current.scrollY;
    be.reset();
    composer.frame(be, Dirty.overlay);
  };
  const measure = (s: (i: number) => void) => {
    const t0 = performance.now();
    for (let i = 0; i < FRAME_WARM || performance.now() - t0 < FRAME_WARM_MS; i++) s(i);
    // Settled after the timed warm-up, as every timed warm-up is (lib/seqstar/AGENTS.md, Tests; LC-20).
    let settled = 0;
    settleFrame(() => s(settled++), FRAME_SETTLE_MS);
    ends = scrolls = 0;
    let i = 0;
    return heapScaling(STEPS, () => s(i++));
  };
  const state = () => ({ owner: t.ctx.ui.drag.value?.owner, ends, active: pd.active, scrolls });
  return { step, measure, state };
}
async function run(mode: 'lines' | 'scroll') {
  const d = await setup(mode);
  const r = d.measure(d.step);
  return { r, ...d.state() };
}

describe('a select drag in text allocates nothing after warm-up', () => {
  it('moves across four lines, past their ends and over their labels bands, the pending range drawn each step', async () => {
    const { r, owner, ends, active } = await run('lines');
    expect([owner, active]).toEqual(['select', true]);
    // The range followed the pointer: its end changed on most steps.
    expect(ends).toBeGreaterThan(STEPS);
    expect(isFlat(r), describeScaling(r, STEPS)).toBe(true);
  }, 120_000);

  it("the scroll past the rows' edges, a frame a step, the range following", async () => {
    const { r, owner, ends, scrolls } = await run('scroll');
    // The rows scrolled on every step but the turns, and the range's end followed them across lines.
    expect(owner).toBe('select');
    expect(scrolls).toBeGreaterThan(9 * STEPS);
    expect(ends).toBeGreaterThan(STEPS);
    expect(isFlat(r), describeScaling(r, STEPS)).toBe(true);
  }, 120_000);

  it('the same moves allocating an object each are caught (positive control)', async () => {
    const d = await setup('lines');
    expect(controlCaught(() => d.measure(leaking(d.step)))).toBe(true);
    expect(d.state().owner).toBe('select');
  }, 120_000);
});
