/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Heap guards for IN2's pointer-rate paths (seqspec:viewer#performance; IN2's hot-path design), on `heapScaling` (a
// step is one input; fractional coordinates, as an engine's are; every default behavior subscribed): moves within an
// element over each region with the labels on, element changes with them off, a points track's items, a wrapped
// section, split rows and matrix bands (`state`); the label churn counted; drags, flat; a full press, drag and release,
// bounded and not retained; the wheel, a million calls, and a rows sweep under a resting pointer, its hover hit again
// each frame. The row zoom's guards are row-zoom.heap.ts'. The cases are input-cases.ts', shared with the quick twin
// (input.quick.ts). The positive controls, one per criterion: moves allocating an object each (`leaking`), primed and
// measured as the guards' are, are caught as not flat; the churn's moves allocating two each by the per-change bound;
// dragPan's cycles allocating eight each by the per-cycle bound.
// The public streams (seqspec:interaction#streams) cost nothing unobserved: every guard runs with them created and no
// subscriber; a storm of context-menu requests unobserved is flat, and element changes with a hover subscriber (which
// allocate a payload each, by design: that run is its positive control) are flat again once it unsubscribed.

import { describe, expect, it } from 'vitest';
import { case3, case10 } from '../../core/_test/format/cases.js';
import { controlCaught, describeScaling, HEAP_DEEP, heapGrowth, heapScaling, isFlat, leak, leaking } from './heap.js';
import {
  CHURNS,
  case3Spots,
  churnCase,
  cycleDoc,
  cycling,
  DRAGS,
  dragDoc,
  elementChanges,
  familySpots,
  featureItems,
  hoverStream,
  matrixBands,
  menuStorm,
  railSpots,
  rowsSweep,
  type Spot,
  scrollGroupDrag,
  startDrag,
  wheelCases,
  within,
  wrappedMoves,
} from './input-cases.js';

const N = HEAP_DEEP ? 50_000 : 25_000;

/**
 * At least 300,000 steps, then chunks until three in a row allocate under 2 KB: V8 runs a large path in its lower tiers, boxing
 * doubles, for several hundred thousand steps before it tiers up, and that is not the path's allocation. A path that
 * never settles within the cap fails the guard that follows.
 */
function prime(step: (i: number) => void) {
  const gc = (globalThis as { gc?: () => void }).gc;
  let k = 0;
  for (let chunk = 0, quiet = 0; chunk < 100 && (quiet < 3 || chunk < 15); chunk++) {
    gc?.();
    const before = process.memoryUsage().heapUsed;
    for (let j = 0; j < 20_000; j++) step(k++);
    quiet = process.memoryUsage().heapUsed - before < 2048 ? quiet + 1 : 0;
  }
}

/**
 * Bytes per unit of `run` once they have settled: V8 compiles an edge-rate path (a press's, a release's) late, its code
 * counts as heap, and how late depends on what ran before, so one measurement after a fixed warm-up is not
 * deterministic (556 or 96 bytes a press-drag-release cycle, by test order). Measured again until two runs agree
 * within 5% (at most ten runs).
 */
function settled(run: () => void, units: number): number {
  let last = heapGrowth(run) / units;
  for (let k = 0; k < 9; k++) {
    const next = heapGrowth(run) / units;
    if (Math.abs(next - last) <= Math.max(8, 0.05 * last)) return next;
    last = next;
  }
  return last;
}

/** `heapScaling` of a step, after it has settled. */
function guard(step: (i: number) => void) {
  prime(step);
  return heapScaling(N, step);
}

/** The label churn's bound on the bytes per hover change at a position. */
const PER_CHANGE = 1280;
/**
 * The label churn's measure, primed: the bytes per hover change (of the step `wrap` makes, the control's), checking that
 * each change computed its label once.
 */
async function churn(...args: Parameters<typeof churnCase>) {
  const c = await churnCase(...args);
  return (wrap: (step: (k: number) => void) => (k: number) => void = (f) => f): number => {
    const step = wrap(c.moved);
    prime(step);
    c.reset();
    const growth = heapGrowth(() => {
      for (let k = 0; k < 8000; k++) step(k);
    });
    // heapGrowth runs it twice: the counters hold both.
    const { computed, changes } = c.counts();
    expect([computed, changes]).toEqual([changes, changes]);
    return growth / (changes / 2);
  };
}
const DOCS = { case3, case10 };

describe('moves and element changes', () => {
  /** Every spot flat, one after another: a failure names its spot and the growth. */
  const flat = (spots: Spot[]) => {
    const failures: string[] = [];
    for (const s of spots) {
      const r = guard(within(s));
      if (!isFlat(r)) failures.push(`${s.name}: ${describeScaling(r, N)}`);
    }
    return failures;
  };

  it('moves within an element over every region of a box section with a tree, labels on: flat, one label per element', async () => {
    const { spots } = await familySpots();
    const regions = new Set(spots.map((s) => s.name.replace(/^family /, '')));
    for (const r of ['region 1', 'region 3 node', 'region 7', 'region 8', 'region 9'])
      expect(regions.has(r), r).toBe(true);
    expect(flat(spots)).toEqual([]);
  }, 120_000);

  it("the same over case 3: a column entry, a data column, the navigator's body, edge and outside, and outside the surface", async () => {
    expect(flat(await case3Spots())).toEqual([]);
  }, 120_000);

  it('and over the rulerRail', async () => {
    const spots = await railSpots();
    expect(spots.some((s) => s.name.includes('region 7'))).toBe(true);
    expect(flat(spots)).toEqual([]);
  }, 120_000);

  it('an element change every 8 moves over the same regions, labels off: flat', async () => {
    const { spots, step, changes } = await elementChanges();
    expect(spots.length).toBeGreaterThan(5);
    const r = guard(step);
    expect(changes()).toBeGreaterThan(N / 8);
    expect(isFlat(r), describeScaling(r, N)).toBe(true);
  }, 120_000);

  it("a points track's items (case 13), a wrapped section, split rows and matrix bands: flat", async () => {
    const items = await featureItems();
    const r1 = guard(items.step);
    expect(items.seen()).toBeGreaterThan(0);
    expect(isFlat(r1), describeScaling(r1, N)).toBe(true);

    // Wrapped, with inline labels.
    const wrapped = await wrappedMoves();
    expect(wrapped.spots).toBeGreaterThan(0);
    expect(wrapped.lines).toBeGreaterThan(5);
    const r2 = guard(wrapped.step);
    expect(isFlat(r2), describeScaling(r2, N)).toBe(true);

    // Matrix bands and stripes: `state` and `part` change with the point.
    const bands = await matrixBands();
    const r3 = guard(bands.step);
    expect(bands.states()).toBeGreaterThan(N / 4);
    expect(isFlat(r3), describeScaling(r3, N)).toBe(true);
  }, 120_000);

  it('label churn: compute runs once per hover change, never on a move; the bytes per change are bounded', async () => {
    const bytes: number[] = [];
    for (const c of CHURNS) bytes.push((await churn(DOCS[c.doc](), c.section, c.track, c.row, 20))());
    console.log(
      `label churn, bytes per hover change: position ${bytes[0].toFixed(0)}, aligned cell ${bytes[1].toFixed(0)}, column ${bytes[2].toFixed(0)}`,
    );
    // Measured: 850, 1,470 and 1,430: the strings and the label's own objects (no set, no closure).
    expect(bytes[0]).toBeLessThan(PER_CHANGE);
    expect(Math.max(bytes[1], bytes[2])).toBeLessThan(2048);
  }, 120_000);
});

/** The bound on the bytes a dragPan's press, drag and release cycle may make. */
const PER_PAN_CYCLE = 512;
/** Bytes per cycle, settled, after 4,000 cycles of warm-up (then `warmed`). */
function perCycle(cycle: (k: number) => void, warmed = () => {}): number {
  for (let k = 0; k < 4000; k++) cycle(k);
  warmed();
  return settled(() => {
    for (let k = 0; k < 2000; k++) cycle(k);
  }, 2000);
}

describe('drags', () => {
  it("dragPan, a select drag, the navigator's body and edge, the overview's bottom grip and the scrollbar's thumb: flat", async () => {
    const t = await dragDoc();
    const failures: string[] = [];
    for (const d of DRAGS) {
      const { step, release } = startDrag(t, d);
      const r = guard(step);
      expect([d.name, release()]).toEqual([d.name, true]);
      if (!isFlat(r)) failures.push(`${d.name}: ${describeScaling(r, N)}`);
    }
    expect(failures).toEqual([]);
  }, 120_000);

  it('dragPan in a scroll group of three: the others follow each move, flat', async () => {
    const { step, release, scrollX } = await scrollGroupDrag();
    const r = guard(step);
    release();
    expect(isFlat(r), describeScaling(r, N)).toBe(true);
    const [one, two, three] = scrollX();
    expect([two, three]).toEqual([one, one]);
  }, 120_000);

  it('a full press, drag and release, repeated: bytes per cycle bounded, nothing retained', async () => {
    const t = await cycleDoc();
    const cycle = cycling(t);
    const gc = (globalThis as { gc?: () => void }).gc as () => void;
    const bytes: Record<string, number> = {};
    for (const [name, alt] of [
      ['dragPan', false],
      ['select', true],
    ] as const) {
      let base = 0;
      bytes[name] = perCycle(
        (k) => cycle(alt, k),
        () => {
          gc();
          gc();
          base = process.memoryUsage().heapUsed;
        },
      );
      gc();
      gc();
      // Retained: what is left after collecting, over the whole run (the measured cycles and the priming).
      expect([name, process.memoryUsage().heapUsed - base < 256 * 1024]).toEqual([name, true]);
    }
    console.log(
      `press, drag and release: ${Object.entries(bytes)
        .map(([k, v]) => `${k} ${v.toFixed(0)} bytes`)
        .join(', ')}`,
    );
    expect(bytes.dragPan).toBeLessThan(PER_PAN_CYCLE);
    expect(bytes.select).toBeLessThan(4096);
    expect(t.ctx.ui.drag.value).toBeNull();
  }, 120_000);
});

describe('the wheel', () => {
  const WN = HEAP_DEEP ? 100_000 : 50_000;
  const wheelGuard = (step: (i: number) => void) => {
    prime(step);
    return heapScaling(WN, step);
  };

  it("Ctrl zoom, a horizontal pan, Shift with only deltaX and with only deltaY, a box section's rows: flat", async () => {
    const { cases, version } = await wheelCases();
    const v0 = version();
    const failures: string[] = [];
    for (const [name, step] of cases) {
      const r = wheelGuard(step);
      if (!isFlat(r)) failures.push(`${name}: ${describeScaling(r, WN)}`);
    }
    expect(version()).toBeGreaterThan(v0 + 4 * WN);
    expect(failures).toEqual([]);
  }, 180_000);

  it("a rows sweep under a resting pointer: the hover re-hit once a frame (the draw's flush) and following, flat", async () => {
    const { t, step, changes } = await rowsSweep();
    const r = wheelGuard(step);
    // The positive control: the hover followed the rows (a change about every other step), the view moved each step.
    expect(changes()).toBeGreaterThan(WN / 4);
    expect(t.ctx.hover.live.current.owner).toBe('pointer');
    expect(isFlat(r), describeScaling(r, WN)).toBe(true);
  }, 180_000);
});

describe('the public streams, unobserved', () => {
  it("a storm of the browser's context-menu requests on a track and off it, no menu subscriber: flat", async () => {
    const { step, asked, observed } = await menuStorm();
    expect(observed()).toBe(false);
    const r = guard(step);
    expect(asked()).toBeGreaterThan(N);
    expect(isFlat(r), describeScaling(r, N)).toBe(true);
  }, 120_000);

  it('element changes with a hover subscriber allocate a payload each (caught); once it unsubscribed: flat', async () => {
    const { step, subscribe, observed } = await hoverStream();
    const sub = subscribe();
    expect(controlCaught(() => guard(step))).toBe(true);
    sub.unsubscribe();
    expect(observed()).toBe(false);
    const r = guard(step);
    expect(isFlat(r), describeScaling(r, N)).toBe(true);
  }, 180_000);
});

// One control per criterion: not flat, bytes per hover change, bytes per cycle.
describe('positive controls', () => {
  it('the same moves allocating an object each are caught', async () => {
    const [s] = await case3Spots();
    expect(controlCaught(() => guard(leaking(within(s))))).toBe(true);
  }, 120_000);

  it("the churn's moves allocating two objects each are caught by the per-change bound", async () => {
    const measure = await churn(case10(), 'chains.main', 'chains.residues', -1, 20);
    expect(controlCaught(() => measure((step) => leaking(leaking(step))), PER_CHANGE)).toBe(true);
  }, 120_000);

  it("dragPan's cycles allocating eight objects each are caught by the per-cycle bound", async () => {
    const t = await cycleDoc();
    const cycle = cycling(t);
    const leaky = (k: number) => {
      cycle(false, k);
      for (let j = 0; j < 8; j++) leak(k + j);
    };
    expect(controlCaught(() => perCycle(leaky), PER_PAN_CYCLE)).toBe(true);
  }, 120_000);
});
