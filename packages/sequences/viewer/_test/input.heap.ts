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
// each frame. The row zoom's guards are rowZoom.heap.ts'. The positive controls, one per criterion, last: moves
// allocating an object each (`leaking`), primed and measured as the guards' are, are caught as not flat; the churn's
// moves allocating two each by the per-change bound; dragPan's cycles allocating eight each by the per-cycle bound.

import { createSeqViewSpec } from "@molstar/seq-core/builder";
import type { SeqViewSpecDoc } from "@molstar/seq-core/format";
import { describe, expect, it } from "vitest";
import { case2, case3, case10, case13 } from "../../core/_test/format/cases.js";
import {
  Hit,
  HitRegion,
  HoverKind,
  Inputs,
  type LabelProvider,
  type LayerRenderer,
  Modifier,
  PointerType,
} from "../src/index.js";
import type { InputDispatcher } from "../src/input.js";
import { hitTree, TREE_INSET, type TreeRows, treeOf, treeSpan } from "../src/tree.js";
import { WindowKind } from "../src/windows.js";
import { FAMILY_SECTION as F, familyWithTree, syntheticFamily } from "./families.js";
import { controlCaught, describeScaling, HEAP_DEEP, heapGrowth, heapScaling, isFlat, leak, leaking } from "./heap.js";
import { type Driven, driven } from "./input-fixtures.js";
import { renderers } from "./layout-fixtures.js";

const N = HEAP_DEEP ? 50_000 : 25_000;
const S = "family.main";
let sink = 0;

const rowsWithTree: LayerRenderer<TreeRows> = {
  name: "rows",
  geometries: ["linear"],
  prepare: (input) => {
    const a = input.section.blocks[0]?.alignment ?? -1;
    const aln = input.set.alignments[a];
    const track = input.section.tracks[input.layer.track];
    const rowOf = Int32Array.from({ length: track.laneCount }, (_, k) => input.section.lanes.row[track.laneStart + k]);
    return {
      set: input.set,
      alignment: a,
      rowOf,
      order: null,
      source: null,
      tree: aln?.tree ? treeOf(aln.tree) : null,
    };
  },
  header: { inset: () => TREE_INSET, hit: (v, s, _p, out) => hitTree(v, s, out), span: treeSpan },
};
const withTree = () => [...renderers().filter((r) => r.name !== "rows"), rowsWithTree as LayerRenderer];

/** A point on a section, and what is there. */
interface Spot {
  name: string;
  t: Driven;
  section: string;
  x: number;
  y: number;
}

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

/** One `PointerInput` for all guards, as an engine's. */
const e = Inputs.pointer();
function move(s: Spot, dx: number, dy: number) {
  e.type = PointerType.move;
  e.section = s.t.ctx.layout.sectionIndex(s.section);
  e.x = s.x + dx;
  e.y = s.y + dy;
  s.t.ctx.input.pointer(e);
}

/** The point is robust: the same element 0.3 px every way. */
function robust(t: Driven, section: string, x: number, y: number) {
  const h = Hit.create();
  const g = Hit.create();
  t.ctx.layout.hitTest(section, x, y, h);
  for (const [dx, dy] of [
    [0.3, 0],
    [-0.3, 0],
    [0, 0.3],
    [0, -0.3],
  ]) {
    t.ctx.layout.hitTest(section, x + dx, y + dy, g);
    if (
      g.region !== h.region ||
      g.entry !== h.entry ||
      g.block !== h.block ||
      g.shown !== h.shown ||
      g.treeNode !== h.treeNode
    )
      return false;
  }
  return true;
}

/** The first robust point of each (region, tree node or not) over a grid of a section. */
function spotsOf(t: Driven, section: string, tag: string): Spot[] {
  const l = t.layout(section);
  const found = new Map<number, Spot>();
  const h = Hit.create();
  for (let y = 1.25; y < l.height; y += 3.5) {
    for (let x = 1.25; x < l.width; x += 4.75) {
      t.ctx.layout.hitTest(section, x, y, h);
      const key = h.region * 4 + (h.treeNode >= 0 ? 1 : 0) + (h.kind === HoverKind.item ? 2 : 0);
      if (h.region === HitRegion.none || found.has(key) || !robust(t, section, x, y)) continue;
      found.set(key, {
        name: `${tag} region ${h.region}${h.treeNode >= 0 ? " node" : ""}${h.kind === HoverKind.item ? " item" : ""}`,
        t,
        section,
        x,
        y,
      });
    }
  }
  return [...found.values()];
}

/** The label churn's bound on the bytes per hover change at a position. */
const PER_CHANGE = 1280;
let computed = 0;
const counting: LabelProvider = {
  name: "lab:count",
  priority: 1,
  applies: (h) => h.kind !== HoverKind.none,
  label: () => void computed++,
};
/**
 * Element changes every 8 moves across 12 cells of a track, primed: its measure gives the bytes per hover change (of
 * the step `wrap` makes, the control's), and checks that each change computed its label once.
 */
async function churn(doc: SeqViewSpecDoc, section: string, track: string, row: number, first: number) {
  const t = await driven(doc, { size: { [section]: [1000, 150] }, registry: [{ labels: [counting] }] });
  const y = row < 0 ? t.yOf(section, track) : t.yOf(section, track, row);
  const xs = Array.from({ length: 12 }, (_, i) => t.xOf(section, 0, first + i * 3));
  let changes = 0;
  // The changes that get a label: an element in a labelled region.
  t.ctx.hover.live.version.subscribe(() => {
    const h = t.ctx.hover.live.current;
    if (h.kind !== HoverKind.none && h.region <= HitRegion.header) changes++;
  });
  const spot: Spot = { name: "", t, section, x: 0, y };
  const moved = (k: number) => {
    spot.x = xs[(k >> 3) % xs.length];
    move(spot, (k & 7) * 0.02, 0);
  };
  return (wrap: (step: (k: number) => void) => (k: number) => void = (f) => f): number => {
    const step = wrap(moved);
    prime(step);
    [changes, computed] = [0, 0];
    const growth = heapGrowth(() => {
      for (let k = 0; k < 8000; k++) step(k);
    });
    // heapGrowth runs it twice: the counters hold both.
    expect([computed, changes]).toEqual([changes, changes]);
    return growth / (changes / 2);
  };
}

/** Moves within spot s's element. */
const within = (s: Spot) => (i: number) => {
  move(s, (i & 7) * 0.03, (i & 3) * 0.03);
  sink = (sink + s.t.ctx.pointer.current.section) & 0xffff;
};

describe("moves and element changes", () => {
  /** Every spot flat, one after another: a failure names its spot and the growth. */
  const flat = (spots: Spot[]) => {
    const failures: string[] = [];
    for (const s of spots) {
      const r = guard(within(s));
      if (!isFlat(r)) failures.push(`${s.name}: ${describeScaling(r, N)}`);
    }
    return failures;
  };

  it("moves within an element over every region of a box section with a tree, labels on: flat, one label per element", async () => {
    const fam = await driven(familyWithTree(160, 20, { columns: true, sites: true }), {
      renderers: withTree(),
      size: { [F]: [1000, 500] },
    });
    const spots = spotsOf(fam, F, "family");
    const regions = new Set(spots.map((s) => s.name.replace(/^family /, "")));
    for (const r of ["region 1", "region 3 node", "region 7", "region 8", "region 9"])
      expect(regions.has(r), r).toBe(true);
    expect(flat(spots)).toEqual([]);
  }, 120_000);

  it("the same over case 3: a column entry, a data column, the navigator's body, edge and outside, and outside the surface", async () => {
    const c3 = await driven(case3(), { size: { [S]: [1000, 150] } });
    const spots = spotsOf(c3, S, "case3");
    const g = c3.frame(S, WindowKind.navigator);
    const nav = c3.layout(S);
    spots.push(
      { name: "navigator body", t: c3, section: S, x: (g[8] + g[9]) / 2, y: nav.navHeight / 2 },
      { name: "navigator edge", t: c3, section: S, x: g[8] + 0.7, y: nav.navHeight / 2 },
      { name: "navigator outside", t: c3, section: S, x: g[8] - 60, y: nav.navHeight / 2 },
      { name: "outside the surface", t: c3, section: S, x: -12.5, y: 40.5 },
    );
    expect(flat(spots)).toEqual([]);
  }, 120_000);

  it("and over the rulerRail", async () => {
    const rail = await driven(case2(), { size: { "fv.main": [1000, 200] } });
    rail.ctx.ui.setStyle("fv.main", "railed");
    await rail.ctx.graph.update();
    const spots = spotsOf(rail, "fv.main", "rail");
    expect(spots.some((s) => s.name.includes("region 7"))).toBe(true);
    expect(flat(spots)).toEqual([]);
  }, 120_000);

  it("an element change every 8 moves over the same regions, labels off: flat", async () => {
    const fam = await driven(familyWithTree(160, 20, { columns: true, sites: true }), {
      renderers: withTree(),
      size: { [F]: [1000, 500] },
      labels: false,
    });
    const spots = spotsOf(fam, F, "family");
    expect(spots.length).toBeGreaterThan(5);
    let changes = 0;
    fam.ctx.hover.live.version.subscribe(() => changes++);
    const r = guard((i) => {
      const s = spots[(i >> 3) % spots.length];
      move(s, (i & 7) * 0.03, 0);
    });
    expect(changes).toBeGreaterThan(N / 8);
    expect(isFlat(r), describeScaling(r, N)).toBe(true);
  }, 120_000);

  it("a points track's items (case 13), a wrapped section, split rows and matrix bands: flat", async () => {
    const f = "features.main";
    const feat = await driven(case13(), { size: { [f]: [1000, 900] }, labels: false });
    const NV = "features.natural-variants";
    const y = feat.yOf(f, "features.binding-site") + 0.25;
    const h = Hit.create();
    let item = 0;
    while (item < 393 && !(feat.ctx.layout.hitTest(f, feat.xOf(f, 0, item), y, h) && h.kind === HoverKind.item)) item++;
    expect(item < 393).toBe(true);
    const items = [
      { x: feat.xOf(f, 0, item), y },
      { x: feat.xOf(f, 0, item + 40), y },
      { x: 40.5, y },
    ];
    // Split rows of one residue column, one after another.
    const x174 = feat.xOf(f, 0, 174);
    const rows = Array.from({ length: 8 }, (_, r) => ({ x: x174, y: feat.yOf(f, NV, r + 3) }));
    let seen = 0;
    feat.ctx.interaction.hover.subscribe((ev) => {
      if (ev.hit.split >= 0) seen++;
    });
    const spot = (p: { x: number; y: number }): Spot => ({ name: "", t: feat, section: f, x: p.x, y: p.y });
    const list = [...items, ...rows].map(spot);
    const r1 = guard((i) => move(list[(i >> 3) % list.length], (i & 7) * 0.02, 0));
    expect(seen).toBeGreaterThan(0);
    expect(isFlat(r1), describeScaling(r1, N)).toBe(true);

    // Wrapped, with inline labels.
    const b = createSeqViewSpec();
    b.data.sequence({ ref: "a", value: "MKVLATGHEK".repeat(30) });
    b.view({ ref: "v", label: "V", layout: "text" }).track({ ref: "v.letters" }).letters();
    const w = await driven(b.getDocument(), { size: { "v.main": [400, 900] }, labels: false });
    const wl = w.layout("v.main");
    const ws = spotsOf(w, "v.main", "wrapped");
    expect(ws.length).toBeGreaterThan(0);
    expect(wl.lines.count).toBeGreaterThan(5);
    const r2 = guard((i) => move(ws[(i >> 3) % ws.length], (i & 7) * 0.03, 0));
    expect(isFlat(r2), describeScaling(r2, N)).toBe(true);

    // Matrix bands and stripes: `state` and `part` change with the point.
    const bands = await driven(case3(), {
      size: { [S]: [1000, 150] },
      labels: false,
      renderers: renderers().map((r) =>
        r.name === "points"
          ? {
              ...r,
              hit: (_v, _s, _p, out) => {
                out.state = Math.floor(out.point[1] / 3) % 4;
                out.part = Math.floor(out.point[0] / 2) % 3;
                return true;
              },
            }
          : r,
      ),
    });
    const by = bands.yOf(S, "family.flvr-arginine");
    const bx = bands.xOf(S, 0, 4);
    let states = 0;
    bands.ctx.hover.live.version.subscribe(() => states++);
    const r3 = guard((i) =>
      move({ name: "", t: bands, section: S, x: bx, y: by }, ((i & 7) - 4) * 0.7, ((i >> 3) & 7) - 4),
    );
    expect(states).toBeGreaterThan(N / 4);
    expect(isFlat(r3), describeScaling(r3, N)).toBe(true);
  }, 120_000);

  it("label churn: compute runs once per hover change, never on a move; the bytes per change are bounded", async () => {
    const position = (await churn(case10(), "chains.main", "chains.residues", -1, 20))();
    const cell = (await churn(case3(), S, "family.members", 0, 20))();
    const column = (await churn(case3(), S, "family.logo", -1, 20))();
    console.log(
      `label churn, bytes per hover change: position ${position.toFixed(0)}, aligned cell ${cell.toFixed(0)}, column ${column.toFixed(0)}`,
    );
    // Measured: 850, 1,470 and 1,430: the strings and the label's own objects (no set, no closure).
    expect(position).toBeLessThan(PER_CHANGE);
    expect(Math.max(cell, column)).toBeLessThan(2048);
  }, 120_000);
});

/** A pressed pointer's moves: buttons down, the press's modifiers, one `PointerInput` reused; the point in `p`. */
function dragged(t: Driven, section: string, mods: { alt?: boolean; shift?: boolean } = {}) {
  e.section = t.ctx.layout.sectionIndex(section);
  e.buttons = 1;
  e.button = -1;
  e.modifiers = (mods.alt ? Modifier.alt : 0) | (mods.shift ? Modifier.shift : 0);
  return (p: Float64Array) => {
    e.type = PointerType.move;
    e.x = p[0];
    e.y = p[1];
    t.ctx.input.pointer(e);
  };
}

/** The bound on the bytes a dragPan's press, drag and release cycle may make. */
const PER_PAN_CYCLE = 512;
/** A full press, four moves and a release at a cell of case 3 (Alt: a select drag), one `PointerInput` reused. */
function cycling(t: Driven) {
  t.vp(S).zoomTo(0, 5, 25);
  const y = t.layout(S).scrolledTop + 9.5;
  const cx = t.xOf(S, 0, 10);
  return (alt: boolean, _k: number) => {
    e.section = 0;
    e.pointerId = 0;
    e.modifiers = alt ? Modifier.alt : 0;
    e.button = 0;
    e.buttons = 1;
    e.type = PointerType.down;
    e.x = cx;
    e.y = y;
    t.ctx.input.pointer(e);
    e.button = -1;
    e.type = PointerType.move;
    for (let k = 1; k <= 4; k++) {
      e.x = cx + k * 9.5;
      t.ctx.input.pointer(e);
    }
    e.button = 0;
    e.buttons = 0;
    e.type = PointerType.up;
    t.ctx.input.pointer(e);
  };
}
/** Bytes per cycle, settled, after 4,000 cycles of warm-up (then `warmed`). */
function perCycle(cycle: (k: number) => void, warmed = () => {}): number {
  for (let k = 0; k < 4000; k++) cycle(k);
  warmed();
  return settled(() => {
    for (let k = 0; k < 2000; k++) cycle(k);
  }, 2000);
}

describe("drags", () => {
  /** Press, start the drag with a first move, check its owner, then the guard over `step`; released after. */
  async function dragFlat(
    name: string,
    t: Driven,
    section: string,
    press: [number, number],
    owner: string,
    wiggle: (i: number, p: Float64Array) => void,
    mods: { alt?: boolean; shift?: boolean } = {},
  ) {
    t.pointer("down", section, press[0], press[1], mods);
    const to = dragged(t, section, mods);
    const p = new Float64Array(2);
    for (const i of [0, 1]) {
      wiggle(i, p);
      to(p);
    }
    expect([name, t.ctx.ui.drag.value?.owner]).toEqual([name, owner]);
    const vp = t.vp(section);
    const before = vp.live.version.value;
    const r = guard((i) => {
      wiggle(i + 2, p);
      to(p);
    });
    t.pointer("up", section, press[0], press[1], mods);
    expect([name, vp.live.version.value > before || owner === "select"]).toEqual([name, true]);
    return isFlat(r) ? null : `${name}: ${describeScaling(r, N)}`;
  }

  it("dragPan, a select drag, the navigator's body and edge, the overview's bottom grip and the scrollbar's thumb: flat", async () => {
    const t = await driven(case3(), { size: { [S]: [1000, 150] } });
    const vp = t.vp(S);
    const y = t.layout(S).scrolledTop + 9.5;
    const failures: (string | null)[] = [];
    vp.zoomTo(0, 10, 20);
    failures.push(
      await dragFlat("dragPan", t, S, [500.5, y], "dragPan", (i, p) => {
        p[0] = 490.5 + (i & 15) * 1.25;
        p[1] = y + (i & 3) * 0.25;
      }),
    );
    // A select drag over cells (Alt): pending changes as it crosses them.
    vp.zoomTo(0, 5, 25);
    const cx = t.xOf(S, 0, 10);
    failures.push(
      await dragFlat(
        "select",
        t,
        S,
        [cx, y],
        "select",
        (i, p) => {
          p[0] = cx + (i % 40) * 3.7;
          p[1] = y + (i & 3) * 0.25;
        },
        { alt: true },
      ),
    );
    // The navigator: its body pans, its left edge zooms.
    vp.zoomTo(0, 5, 25);
    const g = t.frame(S, WindowKind.navigator);
    const ny = t.layout(S).navHeight / 2;
    const mid = (g[8] + g[9]) / 2;
    failures.push(
      await dragFlat("navigator body", t, S, [mid, ny], "navigator", (i, p) => {
        p[0] = mid + ((i & 31) - 16) * 0.9;
        p[1] = ny;
      }),
    );
    const left = t.frame(S, WindowKind.navigator)[8];
    failures.push(
      await dragFlat("navigator edge", t, S, [left + 1, ny], "navigator", (i, p) => {
        p[0] = left + 1 - ((i & 31) - 8) * 0.7;
        p[1] = ny;
      }),
    );
    // The scrollbar's thumb pans the rows (the overview's grips zoom them: rowZoom.heap.ts).
    vp.zoomTo(0, 0, 22);
    const b = t.frame(S, WindowKind.scrollbar);
    const bx = t.layout(S).width - t.layout(S).vnavWidth - 4;
    const by = (b[10] + b[11]) / 2;
    failures.push(
      await dragFlat("scrollbar thumb", t, S, [bx, by], "scrollbar", (i, p) => {
        p[0] = bx;
        p[1] = by + ((i & 31) - 16) * 0.3;
      }),
    );
    expect(failures.filter(Boolean)).toEqual([]);
  }, 120_000);

  it("dragPan in a scroll group of three: the others follow each move, flat", async () => {
    const b = createSeqViewSpec();
    const a = b.data.sequence({ ref: "a", value: "MKVLATGHEKLV".repeat(40) });
    const v = b.view({ ref: "g", label: "G" });
    for (const ref of ["g.one", "g.two", "g.three"])
      v.section({ ref, scroll_group: "s" }).block({ sequence: a }).track().letters();
    const t = await driven(b.getDocument(), { labels: false });
    const [one, two, three] = ["g.one", "g.two", "g.three"].map((r) => t.vp(r));
    one.zoomTo(0, 10, 30);
    expect(t.ctx.scrollGroups.groups.value.get("s")).toHaveLength(3);
    const y = t.layout("g.one").scrolledTop + 8.5;
    t.pointer("down", "g.one", 500.5, y);
    const to = dragged(t, "g.one");
    const p = Float64Array.of(490.5, y);
    to(p);
    p[0] = 480.5;
    to(p);
    expect(t.ctx.ui.drag.value?.owner).toBe("dragPan");
    const r = guard((i) => {
      p[0] = 490.5 + (i & 15) * 1.25;
      to(p);
    });
    t.pointer("up", "g.one", 500.5, y);
    expect(isFlat(r), describeScaling(r, N)).toBe(true);
    expect([two.live.current.scrollX, three.live.current.scrollX]).toEqual([
      one.live.current.scrollX,
      one.live.current.scrollX,
    ]);
  }, 120_000);

  it("a full press, drag and release, repeated: bytes per cycle bounded, nothing retained", async () => {
    const t = await driven(case3(), { size: { [S]: [1000, 150] }, labels: false });
    const cycle = cycling(t);
    const gc = (globalThis as { gc?: () => void }).gc as () => void;
    const bytes: Record<string, number> = {};
    for (const [name, alt] of [
      ["dragPan", false],
      ["select", true],
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
        .join(", ")}`,
    );
    expect(bytes.dragPan).toBeLessThan(PER_PAN_CYCLE);
    expect(bytes.select).toBeLessThan(4096);
    expect(t.ctx.ui.drag.value).toBeNull();
  }, 120_000);
});

describe("the wheel", () => {
  /** The wheel input, written inline by each step: a double passed to a helper is boxed, in the test too. */
  const w = Inputs.wheel();
  const send = (t: Driven) => {
    sink = (sink + (t.ctx.input.wheel(w) ? 1 : 0)) & 0xffff;
  };
  const WN = HEAP_DEEP ? 100_000 : 50_000;
  const wheelGuard = (step: (i: number) => void) => {
    prime(step);
    return heapScaling(WN, step);
  };

  it("Ctrl zoom, a horizontal pan, Shift with only deltaX and with only deltaY, a box section's rows: flat", async () => {
    const t = await driven(case3(), { size: { [S]: [1000, 150] } });
    const vp = t.vp(S);
    vp.zoomTo(0, 5, 20);
    const y = t.layout(S).scrolledTop + 9.5;
    const v0 = vp.live.version.value;
    const failures: string[] = [];
    w.section = t.ctx.layout.sectionIndex(S);
    w.y = y;
    const cases: [string, (i: number) => void][] = [
      [
        "ctrl zoom",
        (i) => {
          w.x = 400.5 + (i & 7);
          w.deltaX = 0.5;
          w.deltaY = i & 1 ? -3.25 : 3.25;
          w.modifiers = Modifier.ctrl;
          send(t);
        },
      ],
      [
        "horizontal pan",
        (i) => {
          w.x = 400.5;
          w.deltaX = i & 1 ? -5.5 : 5.5;
          w.deltaY = 0.25;
          w.modifiers = 0;
          send(t);
        },
      ],
      [
        "shift, deltaX only",
        (i) => {
          w.x = 400.5;
          w.deltaX = i & 1 ? -5.5 : 5.5;
          w.deltaY = 0;
          w.modifiers = Modifier.shift;
          send(t);
        },
      ],
      [
        "shift, deltaY only",
        (i) => {
          w.x = 400.5;
          w.deltaX = 0;
          w.deltaY = i & 1 ? -5.5 : 5.5;
          w.modifiers = Modifier.shift;
          send(t);
        },
      ],
      [
        "rows scroll",
        (i) => {
          w.x = 400.5;
          w.deltaX = 0;
          w.deltaY = i & 1 ? -2.5 : 2.5;
          w.modifiers = 0;
          send(t);
        },
      ],
    ];
    for (const [name, step] of cases) {
      const r = wheelGuard(step);
      if (!isFlat(r)) failures.push(`${name}: ${describeScaling(r, WN)}`);
    }
    expect(vp.live.version.value).toBeGreaterThan(v0 + 4 * WN);
    expect(failures).toEqual([]);
  }, 180_000);

  it("a rows sweep under a resting pointer: the hover re-hit once a frame (the draw's flush) and following, flat", async () => {
    const t = await driven(syntheticFamily(120, 8), { labels: false, size: { [F]: [1000, 260] } });
    const input = t.ctx.input as InputDispatcher;
    const x = t.xOf(F, 0, 4) + 0.25;
    const y = t.layout(F).scrolledTop + 2.5 * t.layout(F).rows.height;
    t.pointer("move", F, x, y);
    w.section = t.ctx.layout.sectionIndex(F);
    w.x = x;
    w.y = y;
    w.deltaX = 0;
    w.modifiers = 0;
    // 32 wheels down a row's half each, 32 up: the row under the pointer changes every other step.
    const step = (i: number) => {
      w.deltaY = (i & 63) < 32 ? 9.5 : -9.5;
      send(t);
      input.flush();
    };
    const v0 = t.ctx.hover.live.version.value;
    const r = wheelGuard(step);
    // The positive control: the hover followed the rows (a change about every other step), the view moved each step.
    const changes = t.ctx.hover.live.version.value - v0;
    expect(changes).toBeGreaterThan(WN / 4);
    expect(t.ctx.hover.live.current.owner).toBe("pointer");
    expect(isFlat(r), describeScaling(r, WN)).toBe(true);
  }, 180_000);
});

// The controls last: cases in one file share V8's feedback. One per criterion: not flat, bytes per hover change, bytes
// per cycle.
describe("positive controls", () => {
  it("the same moves allocating an object each are caught", async () => {
    const c3 = await driven(case3(), { size: { [S]: [1000, 150] } });
    const [s] = spotsOf(c3, S, "case3");
    expect(controlCaught(() => guard(leaking(within(s))))).toBe(true);
  }, 120_000);

  it("the churn's moves allocating two objects each are caught by the per-change bound", async () => {
    const measure = await churn(case10(), "chains.main", "chains.residues", -1, 20);
    expect(controlCaught(() => measure((step) => leaking(leaking(step))), PER_CHANGE)).toBe(true);
  }, 120_000);

  it("dragPan's cycles allocating eight objects each are caught by the per-cycle bound", async () => {
    const t = await driven(case3(), { size: { [S]: [1000, 150] }, labels: false });
    const cycle = cycling(t);
    const leaky = (k: number) => {
      cycle(false, k);
      for (let j = 0; j < 8; j++) leak(k + j);
    };
    expect(controlCaught(() => perCycle(leaky), PER_PAN_CYCLE)).toBe(true);
  }, 120_000);
});
