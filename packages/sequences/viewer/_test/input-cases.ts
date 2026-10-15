/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// IN2's pointer-rate cases, for its heap guards and their quick twin (input.heap.ts, input.quick.ts); every step is one
// input, at fractional coordinates as an engine's are, through one reused `PointerInput` (or `WheelInput`), with every
// default behavior subscribed and the public streams created with no subscriber: moves within an element at robust
// points of every region (`spotsOf`, `within`); element changes every 8 moves (`elementChanges`); a points track's
// items and split rows, a wrapped section, matrix bands (`featureItems`, `wrappedMoves`, `matrixBands`); the label
// churn's moves (`churnCase`); drags held at their owner (`DRAGS`, `startDrag`, `scrollGroupDrag`); a full press, drag
// and release (`cycling`); the wheel (`wheelCases`, `rowsSweep`); context-menu requests unobserved (`menuStorm`).

import { createSeqViewSpec } from '@molstar/seq-core/builder';
import type { SeqViewSpecDoc } from '@molstar/seq-core/format';
import { expect } from 'vitest';
import { case2, case3, case13 } from '../../core/_test/format/cases.js';
import {
  Hit,
  HitRegion,
  HoverKind,
  Inputs,
  type LabelProvider,
  type LayerRenderer,
  Modifier,
  PointerType,
} from '../src/index.js';
import type { InputDispatcher } from '../src/input.js';
import { hitTree, TREE_INSET, type TreeRows, treeOf, treeSpan } from '../src/tree.js';
import { WindowKind } from '../src/windows.js';
import { FAMILY_SECTION as F, familyWithTree, syntheticFamily } from './families.js';
import { type Driven, driven } from './input-fixtures.js';
import { renderers } from './layout-fixtures.js';

export const S = 'family.main';
let sink = 0;

const rowsWithTree: LayerRenderer<TreeRows> = {
  name: 'rows',
  geometries: ['linear'],
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
const withTree = () => [...renderers().filter((r) => r.name !== 'rows'), rowsWithTree as LayerRenderer];

/** A point on a section, and what is there. */
export interface Spot {
  name: string;
  t: Driven;
  section: string;
  x: number;
  y: number;
}

/** One `PointerInput` for all cases, as an engine's. */
const e = Inputs.pointer();
export function move(s: Spot, dx: number, dy: number) {
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
export function spotsOf(t: Driven, section: string, tag: string): Spot[] {
  const l = t.layout(section);
  const found = new Map<number, Spot>();
  const h = Hit.create();
  for (let y = 1.25; y < l.height; y += 3.5) {
    for (let x = 1.25; x < l.width; x += 4.75) {
      t.ctx.layout.hitTest(section, x, y, h);
      const key = h.region * 4 + (h.treeNode >= 0 ? 1 : 0) + (h.kind === HoverKind.item ? 2 : 0);
      if (h.region === HitRegion.none || found.has(key) || !robust(t, section, x, y)) continue;
      found.set(key, {
        name: `${tag} region ${h.region}${h.treeNode >= 0 ? ' node' : ''}${h.kind === HoverKind.item ? ' item' : ''}`,
        t,
        section,
        x,
        y,
      });
    }
  }
  return [...found.values()];
}

/** Moves within spot s's element. */
export const within = (s: Spot) => (i: number) => {
  move(s, (i & 7) * 0.03, (i & 3) * 0.03);
  sink = (sink + s.t.ctx.pointer.current.section) & 0xffff;
};

/** A family with a tree, columns and sites, in a box section; its spots over every region. */
export async function familySpots(labels = true) {
  const t = await driven(familyWithTree(160, 20, { columns: true, sites: true }), {
    renderers: withTree(),
    size: { [F]: [1000, 500] },
    ...(labels ? {} : { labels: false as const }),
  });
  return { t, spots: spotsOf(t, F, 'family') };
}

/** Case 3's spots, and its navigator's body, edge and outside, and a point outside the surface. */
export async function case3Spots() {
  const c3 = await driven(case3(), { size: { [S]: [1000, 150] } });
  const spots = spotsOf(c3, S, 'case3');
  const g = c3.frame(S, WindowKind.navigator);
  const nav = c3.layout(S);
  spots.push(
    { name: 'navigator body', t: c3, section: S, x: (g[8] + g[9]) / 2, y: nav.navHeight / 2 },
    { name: 'navigator edge', t: c3, section: S, x: g[8] + 0.7, y: nav.navHeight / 2 },
    { name: 'navigator outside', t: c3, section: S, x: g[8] - 60, y: nav.navHeight / 2 },
    { name: 'outside the surface', t: c3, section: S, x: -12.5, y: 40.5 },
  );
  return spots;
}

/** Case 2 on the `railed` style: the rulerRail's spots. */
export async function railSpots() {
  const rail = await driven(case2(), { size: { 'fv.main': [1000, 200] } });
  rail.ctx.ui.setStyle('fv.main', 'railed');
  await rail.ctx.graph.update();
  return spotsOf(rail, 'fv.main', 'rail');
}

/** An element change every 8 moves over the family's regions, labels off; `changes` counts the hover's versions. */
export async function elementChanges() {
  const { t, spots } = await familySpots(false);
  let changes = 0;
  t.ctx.hover.live.version.subscribe(() => changes++);
  const step = (i: number) => {
    const s = spots[(i >> 3) % spots.length];
    move(s, (i & 7) * 0.03, 0);
  };
  return { t, spots, step, changes: () => changes };
}

/** A points track's items (case 13) and split rows of one residue column; `seen`: hover events on a split row. */
export async function featureItems() {
  const f = 'features.main';
  const feat = await driven(case13(), { size: { [f]: [1000, 900] }, labels: false });
  const NV = 'features.natural-variants';
  const y = feat.yOf(f, 'features.binding-site') + 0.25;
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
  feat.ctx.input.events.hover.subscribe((ev) => {
    if (ev.hit.split >= 0) seen++;
  });
  const spot = (p: { x: number; y: number }): Spot => ({ name: '', t: feat, section: f, x: p.x, y: p.y });
  const list = [...items, ...rows].map(spot);
  return { step: (i: number) => move(list[(i >> 3) % list.length], (i & 7) * 0.02, 0), seen: () => seen };
}

/** A wrapped section with inline labels: moves over its spots; `spots`, `lines`: their counts. */
export async function wrappedMoves() {
  const b = createSeqViewSpec();
  b.data.sequence({ ref: 'a', value: 'MKVLATGHEK'.repeat(30) });
  b.view({ ref: 'v', label: 'V', layout: 'text' }).track({ ref: 'v.letters' }).letters();
  const w = await driven(b.getDocument(), { size: { 'v.main': [400, 900] }, labels: false });
  const ws = spotsOf(w, 'v.main', 'wrapped');
  return {
    step: (i: number) => move(ws[(i >> 3) % ws.length], (i & 7) * 0.03, 0),
    spots: ws.length,
    lines: w.layout('v.main').lines.count,
  };
}

/** Matrix bands and stripes: `state` and `part` change with the point; `states` counts the hover's versions. */
export async function matrixBands() {
  const bands = await driven(case3(), {
    size: { [S]: [1000, 150] },
    labels: false,
    renderers: renderers().map((r) =>
      r.name === 'points'
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
  const by = bands.yOf(S, 'family.flvr-arginine');
  const bx = bands.xOf(S, 0, 4);
  let states = 0;
  bands.ctx.hover.live.version.subscribe(() => states++);
  const at: Spot = { name: '', t: bands, section: S, x: bx, y: by };
  return { step: (i: number) => move(at, ((i & 7) - 4) * 0.7, ((i >> 3) & 7) - 4), states: () => states };
}

let computed = 0;
const counting: LabelProvider = {
  name: 'lab:count',
  priority: 1,
  applies: (h) => h.kind !== HoverKind.none,
  label: () => void computed++,
};

/** The label churn's tracks: a position, an aligned cell, a column (each with its document, section, track, row). */
export const CHURNS = [
  { name: 'position', doc: 'case10', section: 'chains.main', track: 'chains.residues', row: -1 },
  { name: 'aligned cell', doc: 'case3', section: S, track: 'family.members', row: 0 },
  { name: 'column', doc: 'case3', section: S, track: 'family.logo', row: -1 },
] as const;

/**
 * Element changes every 8 moves across 12 cells of a track from position `first`, labels on (a counting provider):
 * `moved` is one move; `counts` the labels computed and the changes that get one (an element in a labelled region)
 * since `reset`.
 */
export async function churnCase(doc: SeqViewSpecDoc, section: string, track: string, row: number, first: number) {
  const t = await driven(doc, { size: { [section]: [1000, 150] }, registry: [{ labels: [counting] }] });
  const y = row < 0 ? t.yOf(section, track) : t.yOf(section, track, row);
  const xs = Array.from({ length: 12 }, (_, i) => t.xOf(section, 0, first + i * 3));
  let changes = 0;
  t.ctx.hover.live.version.subscribe(() => {
    const h = t.ctx.hover.live.current;
    if (h.kind !== HoverKind.none && h.region <= HitRegion.header) changes++;
  });
  const spot: Spot = { name: '', t, section, x: 0, y };
  const moved = (k: number) => {
    spot.x = xs[(k >> 3) % xs.length];
    move(spot, (k & 7) * 0.02, 0);
  };
  return {
    moved,
    counts: () => ({ computed, changes }),
    reset: () => {
      [changes, computed] = [0, 0];
    },
  };
}

/** A pressed pointer's moves: buttons down, the press's modifiers, one `PointerInput` reused; the point in `p`. */
export function dragged(t: Driven, section: string, mods: { alt?: boolean; shift?: boolean } = {}) {
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

/** A drag over case 3: its owner, modifiers, and where it presses and moves (`at`, after the view it sets). */
export interface DragCase {
  name: string;
  owner: string;
  mods?: { alt?: boolean };
  at(t: Driven): { press: [number, number]; wiggle: (i: number, p: Float64Array) => void };
}

/** Case 3 for the drags. */
export const dragDoc = () => driven(case3(), { size: { [S]: [1000, 150] } });

/** dragPan, a select drag, the navigator's body and edge, the scrollbar's thumb, in order on one viewer. */
export const DRAGS: DragCase[] = [
  {
    name: 'dragPan',
    owner: 'dragPan',
    at: (t) => {
      t.vp(S).zoomTo(0, 10, 20);
      const y = t.layout(S).scrolledTop + 9.5;
      return {
        press: [500.5, y],
        wiggle: (i, p) => {
          p[0] = 490.5 + (i & 15) * 1.25;
          p[1] = y + (i & 3) * 0.25;
        },
      };
    },
  },
  {
    // A select drag over cells (Alt): pending changes as it crosses them.
    name: 'select',
    owner: 'select',
    mods: { alt: true },
    at: (t) => {
      t.vp(S).zoomTo(0, 5, 25);
      const y = t.layout(S).scrolledTop + 9.5;
      const cx = t.xOf(S, 0, 10);
      return {
        press: [cx, y],
        wiggle: (i, p) => {
          p[0] = cx + (i % 40) * 3.7;
          p[1] = y + (i & 3) * 0.25;
        },
      };
    },
  },
  {
    // The navigator: its body pans.
    name: 'navigator body',
    owner: 'navigator',
    at: (t) => {
      t.vp(S).zoomTo(0, 5, 25);
      const g = t.frame(S, WindowKind.navigator);
      const ny = t.layout(S).navHeight / 2;
      const mid = (g[8] + g[9]) / 2;
      return {
        press: [mid, ny],
        wiggle: (i, p) => {
          p[0] = mid + ((i & 31) - 16) * 0.9;
          p[1] = ny;
        },
      };
    },
  },
  {
    // Its left edge zooms.
    name: 'navigator edge',
    owner: 'navigator',
    at: (t) => {
      const left = t.frame(S, WindowKind.navigator)[8];
      const ny = t.layout(S).navHeight / 2;
      return {
        press: [left + 1, ny],
        wiggle: (i, p) => {
          p[0] = left + 1 - ((i & 31) - 8) * 0.7;
          p[1] = ny;
        },
      };
    },
  },
  {
    // The scrollbar's thumb pans the rows (the overview's grips zoom them: row-zoom.heap.ts).
    name: 'scrollbar thumb',
    owner: 'scrollbar',
    at: (t) => {
      t.vp(S).zoomTo(0, 0, 22);
      const b = t.frame(S, WindowKind.scrollbar);
      const bx = t.layout(S).width - t.layout(S).vnavWidth - 4;
      const by = (b[10] + b[11]) / 2;
      return {
        press: [bx, by],
        wiggle: (i, p) => {
          p[0] = bx;
          p[1] = by + ((i & 31) - 16) * 0.3;
        },
      };
    },
  },
];

/**
 * A drag pressed and started with two moves, its owner checked; `step` moves it on, `release` lets go and says whether
 * it moved the view (a select drag moves none).
 */
export function startDrag(t: Driven, d: DragCase) {
  const { press, wiggle } = d.at(t);
  t.pointer('down', S, press[0], press[1], d.mods);
  const to = dragged(t, S, d.mods);
  const p = new Float64Array(2);
  for (const i of [0, 1]) {
    wiggle(i, p);
    to(p);
  }
  expect([d.name, t.ctx.ui.drag.value?.owner]).toEqual([d.name, d.owner]);
  const vp = t.vp(S);
  const before = vp.live.version.value;
  return {
    step: (i: number) => {
      wiggle(i + 2, p);
      to(p);
    },
    release: () => {
      t.pointer('up', S, press[0], press[1], d.mods);
      return vp.live.version.value > before || d.owner === 'select';
    },
  };
}

/** dragPan in a scroll group of three, pressed and started; `scrollX` the three sections' after. */
export async function scrollGroupDrag() {
  const b = createSeqViewSpec();
  const a = b.data.sequence({ ref: 'a', value: 'MKVLATGHEKLV'.repeat(40) });
  const v = b.view({ ref: 'g', label: 'G' });
  for (const ref of ['g.one', 'g.two', 'g.three'])
    v.section({ ref, scroll_group: 's' }).block({ sequence: a }).track().letters();
  const t = await driven(b.getDocument(), { labels: false });
  const [one, two, three] = ['g.one', 'g.two', 'g.three'].map((r) => t.vp(r));
  one.zoomTo(0, 10, 30);
  expect(t.ctx.scrollGroups.groups.value.get('s')).toHaveLength(3);
  const y = t.layout('g.one').scrolledTop + 8.5;
  t.pointer('down', 'g.one', 500.5, y);
  const to = dragged(t, 'g.one');
  const p = Float64Array.of(490.5, y);
  to(p);
  p[0] = 480.5;
  to(p);
  expect(t.ctx.ui.drag.value?.owner).toBe('dragPan');
  return {
    step: (i: number) => {
      p[0] = 490.5 + (i & 15) * 1.25;
      to(p);
    },
    release: () => t.pointer('up', 'g.one', 500.5, y),
    scrollX: () => [one.live.current.scrollX, two.live.current.scrollX, three.live.current.scrollX],
  };
}

/** Case 3 for the press, drag and release cycles, labels off. */
export const cycleDoc = () => driven(case3(), { size: { [S]: [1000, 150] }, labels: false });

/** A full press, four moves and a release at a cell of case 3 (Alt: a select drag), one `PointerInput` reused. */
export function cycling(t: Driven) {
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

/** The wheel input, written inline by each step: a double passed to a helper is boxed, in the test too. */
const w = Inputs.wheel();
const send = (t: Driven) => {
  sink = (sink + (t.ctx.input.wheel(w) ? 1 : 0)) & 0xffff;
};

/**
 * Case 3 zoomed, the wheel over its rows: Ctrl zoom, a horizontal pan, Shift with only deltaX and with only deltaY, the
 * box section's rows; `version` is the viewport's.
 */
export async function wheelCases() {
  const t = await driven(case3(), { size: { [S]: [1000, 150] } });
  const vp = t.vp(S);
  vp.zoomTo(0, 5, 20);
  w.section = t.ctx.layout.sectionIndex(S);
  w.y = t.layout(S).scrolledTop + 9.5;
  const cases: [string, (i: number) => void][] = [
    [
      'ctrl zoom',
      (i) => {
        w.x = 400.5 + (i & 7);
        w.deltaX = 0.5;
        w.deltaY = i & 1 ? -3.25 : 3.25;
        w.modifiers = Modifier.ctrl;
        send(t);
      },
    ],
    [
      'horizontal pan',
      (i) => {
        w.x = 400.5;
        w.deltaX = i & 1 ? -5.5 : 5.5;
        w.deltaY = 0.25;
        w.modifiers = 0;
        send(t);
      },
    ],
    [
      'shift, deltaX only',
      (i) => {
        w.x = 400.5;
        w.deltaX = i & 1 ? -5.5 : 5.5;
        w.deltaY = 0;
        w.modifiers = Modifier.shift;
        send(t);
      },
    ],
    [
      'shift, deltaY only',
      (i) => {
        w.x = 400.5;
        w.deltaX = 0;
        w.deltaY = i & 1 ? -5.5 : 5.5;
        w.modifiers = Modifier.shift;
        send(t);
      },
    ],
    [
      'rows scroll',
      (i) => {
        w.x = 400.5;
        w.deltaX = 0;
        w.deltaY = i & 1 ? -2.5 : 2.5;
        w.modifiers = 0;
        send(t);
      },
    ],
  ];
  return { cases, version: () => vp.live.version.value };
}

/**
 * A rows sweep under a resting pointer: 32 wheels down a row's half each, 32 up, each flushed (the draw's flush hits
 * the hover again), so the row under the pointer changes every other step; `changes` counts the hover's versions.
 */
export async function rowsSweep() {
  const t = await driven(syntheticFamily(120, 8), { labels: false, size: { [F]: [1000, 260] } });
  const input = t.ctx.input as InputDispatcher;
  const x = t.xOf(F, 0, 4) + 0.25;
  const y = t.layout(F).scrolledTop + 2.5 * t.layout(F).rows.height;
  t.pointer('move', F, x, y);
  w.section = t.ctx.layout.sectionIndex(F);
  w.x = x;
  w.y = y;
  w.deltaX = 0;
  w.modifiers = 0;
  const step = (i: number) => {
    w.deltaY = (i & 63) < 32 ? 9.5 : -9.5;
    send(t);
    input.flush();
  };
  const v0 = t.ctx.hover.live.version.value;
  return { t, step, changes: () => t.ctx.hover.live.version.value - v0 };
}

/**
 * A storm of the browser's context-menu requests on case 3's spots, on a track and off it, with no `interaction.menu`
 * subscriber; `asked` counts the dispatcher's own menu events.
 */
export async function menuStorm() {
  const c3 = await driven(case3(), { size: { [S]: [1000, 150] }, labels: false });
  const spots = spotsOf(c3, S, 'case3');
  const m = Inputs.pointer();
  m.section = c3.ctx.layout.sectionIndex(S);
  m.button = 2;
  m.buttons = 2;
  m.pointerId = 1;
  let asked = 0;
  c3.ctx.input.events.menu.subscribe(() => asked++);
  const step = (i: number) => {
    const s = spots[(i >> 3) % spots.length];
    m.x = s.x + (i & 7) * 0.03;
    m.y = s.y;
    c3.ctx.input.menu(m);
  };
  return { step, asked: () => asked, observed: () => c3.ctx.interaction.observed('menu') };
}

/** Element changes for the public hover stream: `subscribe` adds a subscriber (a payload a change, by design). */
export async function hoverStream() {
  const { t, spots } = await familySpots(false);
  expect(spots.length).toBeGreaterThan(5);
  return {
    step: (i: number) => move(spots[(i >> 3) % spots.length], (i & 7) * 0.03, 0),
    subscribe: () =>
      t.ctx.interaction.hover.subscribe((ev) => {
        sink = (sink + ev.place.x) & 0xffff;
      }),
    observed: () => t.ctx.interaction.observed('hover'),
  };
}
