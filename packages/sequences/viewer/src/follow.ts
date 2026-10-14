/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The `follow` behavior (seqspec:viewer#behaviors; seqspec:format#sections; LK1): sections that follow each other set
// their windows through a coordinate map. When a section moves (a gesture, the keyboard, a host's or the api's move),
// the element at its center goes through the route of the pair from its block (`ResolvedView.follows`) to the other
// section, which centers that image, with the fraction of the cell, at the leader's cell times the route's
// scale (a codon step: a third, or three). Followers move with origin `behavior`, which nothing follows, so no change
// echoes; chains and scroll groups are walked from the leader, each section once. A change that is not a move of the
// section's own (origin `graph`: a new axis, its data loaded; or a resize) leads nothing: the section is set again from
// the last section that led in its links, so a follower whose data loads late opens where the leader is, and the
// leader stays; but a section whose own focus a patch changed leads (as a snapshot's focus does). A section on rings
// shows its whole axis: a pass-through, never moved and never leading from its own viewport (CG3; its window leads its
// linear section). Vertical scroll is each section's own. Per change nothing is allocated: the tables are built when
// the view, the scroll groups or the sections change, and doubles pass through Float64Arrays.

import type { ResolvedView } from "@molstar/seq-core/format";
import type { CoordinateMap } from "@molstar/seq-core/model";
import type { Axis } from "./axis.ts";
import { defineBehavior } from "./default-behaviors.ts";
import type { Unsubscribable } from "./state.ts";
import { type SectionViewport, Viewport } from "./viewport.ts";

/**
 * `io[0]`, a position in the map's `from` system (an element and the fraction of it, in index order), → its image in
 * `to`, likewise; where it maps nowhere (a gap, an intron, past the ends), the middle of the nearest mapped element.
 * False: an empty map. One search over the segments; nothing allocated.
 */
export function imageForward(m: CoordinateMap, io: Float64Array): boolean {
  const n = m.count;
  if (n === 0) return false;
  let q = io[0] * m.span;
  const e = Math.floor(q);
  let lo = 0;
  let hi = n;
  while (lo < hi) {
    const mid = (lo + hi) >> 1;
    if (m.at[mid] <= e) lo = mid + 1;
    else hi = mid;
  }
  let k = lo - 1;
  // From ranges overlap only in `order` 2 (a source used twice): the first segment back that holds e.
  for (let j = k; m.order === 2 && j >= 0; j--) {
    if (e < m.at[j] + m.length[j]) {
      k = j;
      break;
    }
  }
  if (k < 0 || e >= m.at[k] + m.length[k]) {
    const next = k + 1;
    if (k < 0 || (next < n && m.at[next] - e < e - (m.at[k] + m.length[k] - 1))) {
      k = next;
      q = m.at[k] + 0.5;
    } else q = m.at[k] + m.length[k] - 0.5;
  }
  const d = q - m.at[k];
  io[0] = m.dir[k] > 0 ? m.toLo[k] + d : m.toLo[k] + m.length[k] - d;
  return true;
}

/** The first segment in `byTo` order holding to position j (copies in from order); -1: none. `h` out: the last start at or below j. */
function holding(m: CoordinateMap, j: number, h: Int32Array): number {
  let lo = 0;
  let hi = m.count;
  while (lo < hi) {
    const mid = (lo + hi) >> 1;
    if (m.toLo[m.byTo[mid]] <= j) lo = mid + 1;
    else hi = mid;
  }
  h[0] = lo - 1;
  let s = -1;
  for (let x = lo - 1; x >= 0; x--) {
    const y = m.byTo[x];
    if (j < m.toLo[y] + m.length[y]) s = y;
    if (x > 0 && m.byToMax[x - 1] <= j) break;
  }
  return s;
}
const below = new Int32Array(1);

/** `imageForward` read back: `io[0]` in `to` → in `from` (a position two copies hold: the first in from order). */
export function imageBackward(m: CoordinateMap, io: Float64Array): boolean {
  const n = m.count;
  if (n === 0) return false;
  let t = io[0];
  const j = Math.floor(t);
  let s = holding(m, j, below);
  // A strand - segment's image is closed at its high edge: t on it is that segment's first element's low edge.
  if (s < 0 && t === j && j > 0) {
    const r = holding(m, j - 1, below);
    if (r >= 0 && m.dir[r] < 0 && j === m.toLo[r] + m.length[r]) s = r;
    else holding(m, j, below);
  }
  if (s < 0) {
    // The segment below j (by its start: the one ending last when to ranges are disjoint) or the next above it.
    const h = below[0];
    const lo = h >= 0 ? m.byTo[h] : -1;
    const up = h + 1 < n ? m.byTo[h + 1] : -1;
    if (lo < 0 || (up >= 0 && m.toLo[up] - j < j - (m.toLo[lo] + m.length[lo] - 1))) {
      s = up;
      t = m.toLo[s] + 0.5;
    } else {
      s = lo;
      t = m.toLo[s] + m.length[s] - 0.5;
    }
  }
  const d = m.dir[s] > 0 ? t - m.toLo[s] : m.toLo[s] + m.length[s] - t;
  io[0] = (m.at[s] + d) / m.span;
  return true;
}

/** Of the follower's blocks showing the route's target system, the one whose piece shows index `io[0]`; else `to`. */
function blockFor(ax: Axis | null, to: number, io: Float64Array): number {
  if (!ax || to < 0 || to >= ax.count) return to;
  const i = Math.floor(io[0]);
  for (let b = 0; b < ax.count; b++) {
    if (ax.unit[b] === ax.unit[to] && ax.alignment[b] === ax.alignment[to] && ax.shownOf(b, i) >= 0) return b;
  }
  return to;
}

/**
 * Block `b`'s index `io[0]` (an element and its fraction, as `centerElement` gives it) and px a position `io[1]`
 * through link `k` of the view's follows (`2 * pair + side`, side 1: from the pair's b) onto the other section: its
 * block there (of `to`'s axis, the one showing the image), `io` its index and px a position there; -1: no route from
 * that block (`io` then as it was). Also the ring's window (ring-window.ts). Allocation-free.
 */
export function followImage(v: ResolvedView, k: number, b: number, io: Float64Array, to: Axis | null): number {
  const f = v.follows;
  const p = k >> 1;
  if (b < 0 || b >= v.sections.blockCount[k & 1 ? f.b[p] : f.a[p]]) return -1;
  const r = f.routeStart[p] + (k & 1 ? v.sections.blockCount[f.a[p]] : 0) + b;
  const path = f.routes.path[r];
  if (path < 0) return -1;
  const m = v.paths[path].map;
  const forward = f.routes.forward[r] === 1;
  if (!(forward ? imageForward(m, io) : imageBackward(m, io))) return -1;
  io[1] = forward ? io[1] / m.span : io[1] * m.span;
  return blockFor(to, f.routes.to[r], io);
}

const GROUP = -1;

export const FollowBehavior = defineBehavior(
  "follow",
  "view",
  0,
  {},
  "Sections that follow each other (`follow`) set their windows through a coordinate map: the other centers the image " +
    "of the leader's center element at the leader's cell times the route's scale; each keeps its vertical scroll.",
  ({ ctx, add }) => {
    /** The view the tables were built from: a change while the graph runs ahead of them leads nothing. */
    let view: ResolvedView | null = null;
    let nodes: (SectionViewport | null)[] = [];
    /** CSR of each node's links: the other node, and the pair and side (`2 * pair + side`, side 1: from b) or GROUP. */
    let start = new Int32Array(1);
    let other = new Int32Array(0);
    let kind = new Int32Array(0);
    /** Per node: the cell, scrollX and width last seen (a bump with none changed is not a move). */
    let seen = new Float64Array(0);
    /** Per node: when it last led (`clock`), and 1 for a section that follows none (the authority among equals). */
    let led = new Float64Array(0);
    let root = new Uint8Array(0);
    let mark = new Int32Array(0);
    let queue = new Int32Array(0);
    let epoch = 0;
    let clock = 0;
    let busy = false;
    const io = new Float64Array(2);
    /**
     * Per node on rings: the image it was last reached at (its block, index and px a position; block -1: none). A ring
     * shows its whole axis, so it is a pass-through: reached, it keeps the image and is never moved; leading (as the
     * middle of a chain), it passes that image on, never its own viewport (CG3: its window leads instead).
     */
    let kept = new Float64Array(0);
    let subs: Unsubscribable[] = [];
    /** Sections whose focus a patch changed before the tables held them: they led last (`rebuild` marks them). */
    const focused = new Set<SectionViewport>();
    /** Whether node n's section is on rings (its preset's geometry: no layout is built for it). */
    const ring = (n: number): boolean => {
      const vp = nodes[n];
      return vp !== null && ctx.layout.presetOf(ctx.layout.sectionIndex(vp.ref))?.geometry === "circular";
    };

    /** Sets y from x along one link; false when x has no center or no route from there (y stays where it is). */
    const apply = (x: number, y: number, k: number): boolean => {
      const from = nodes[x];
      const to = nodes[y];
      if (!from || !to) return false;
      if (k === GROUP) {
        if (ring(x) || ring(y)) return false;
        to.follow(from, "behavior");
        return true;
      }
      const v = view;
      if (!v || !(from instanceof Viewport) || !(to instanceof Viewport)) return false;
      let a = -1;
      if (!ring(x)) a = from.centerElement(io);
      else if (kept[3 * x] >= 0) {
        a = kept[3 * x];
        io[0] = kept[3 * x + 1];
        io[1] = kept[3 * x + 2];
      }
      const b = followImage(v, k, a, io, to.axis);
      if (b < 0) return false;
      if (!ring(y)) {
        to.centerOn(b, io, "behavior");
        return true;
      }
      kept[3 * y] = b;
      kept[3 * y + 1] = io[0];
      kept[3 * y + 2] = io[1];
      return true;
    };
    /** From x, every node its links reach, each set once (a node a route does not reach may be reached by another). */
    const propagate = (x: number) => {
      busy = true;
      try {
        epoch++;
        mark[x] = epoch;
        queue[0] = x;
        let tail = 1;
        for (let head = 0; head < tail; head++) {
          const u = queue[head];
          for (let e = start[u]; e < start[u + 1]; e++) {
            const y = other[e];
            if (mark[y] === epoch || !apply(u, y, kind[e])) continue;
            mark[y] = epoch;
            queue[tail++] = y;
          }
        }
      } finally {
        busy = false;
      }
    };
    /** A section that can lead a resync: shown with a data width, and not on rings (a ring leads only through its window). */
    const ready = (n: number) => {
      const vp = nodes[n];
      if (ring(n)) return false;
      const ax = vp instanceof Viewport ? vp.axis : null;
      return ax !== null && ax.count > 0 && (vp as Viewport).live.current.width - ax.insetX > 0;
    };
    /** n and the nodes its links reach, into `queue`: their count. */
    const component = (n: number): number => {
      epoch++;
      mark[n] = epoch;
      queue[0] = n;
      let tail = 1;
      for (let head = 0; head < tail; head++) {
        const u = queue[head];
        for (let e = start[u]; e < start[u + 1]; e++) {
          if (mark[other[e]] === epoch) continue;
          mark[other[e]] = epoch;
          queue[tail++] = other[e];
        }
      }
      return tail;
    };
    /** n's links set again from their authority: the ready node that led last (else one that follows none). */
    const resync = (n: number) => {
      const count = component(n);
      let best = -1;
      for (let h = 0; h < count; h++) {
        const u = queue[h];
        if (ready(u) && (best < 0 || led[u] > led[best] || (led[u] === led[best] && root[u] > root[best]))) best = u;
      }
      if (best >= 0) propagate(best);
    };
    const onVersion = (n: number) => () => {
      const vp = nodes[n];
      if (!vp) return;
      const s = vp.live.current;
      const o = 3 * n;
      const moved = s.cell !== seen[o] || s.scrollX !== seen[o + 1];
      const resized = s.width !== seen[o + 2];
      seen[o] = s.cell;
      seen[o + 1] = s.scrollX;
      seen[o + 2] = s.width;
      // A change the graph made before these tables were rebuilt is resynced when they are.
      if (busy || (!moved && !resized) || s.origin === "behavior" || ctx.graph.assembled.value?.view?.resolved !== view)
        return;
      if (s.origin === "graph" || resized) return resync(n);
      // A ring's own viewport leads nothing (its keys and gestures go to its window's section).
      if (ring(n)) return;
      led[n] = ++clock;
      propagate(n);
    };

    const release = () => {
      for (const s of subs) s.unsubscribe();
      subs = [];
    };
    /** The tables, at edge rate: the current view's sections and pairs, scroll groups (joins included). */
    const rebuild = () => {
      release();
      const before = new Map(nodes.map((vp, i) => [vp, i]));
      const [oldSeen, oldLed] = [seen, led];
      view = ctx.graph.assembled.value?.view?.resolved ?? null;
      const set = ctx.graph.assembled.value?.view?.set;
      const v = view;
      const sec = v?.sections;
      nodes = [];
      for (let i = 0; sec && set && i < sec.count; i++) nodes.push(ctx.viewports.get(set.strings.get(sec.ref[i])));
      const index = new Map(nodes.map((vp, i) => [vp, i]));
      const links: [number, number, number][] = [];
      for (let p = 0; v && p < v.follows.count; p++) {
        links.push([v.follows.a[p], v.follows.b[p], 2 * p], [v.follows.b[p], v.follows.a[p], 2 * p + 1]);
      }
      // A scroll group joins the walk where a member follows or is followed (the rest is `scrollGroups`' alone).
      const paired = new Set(links.map((l) => nodes[l[0]]));
      for (const members of ctx.scrollGroups.groups.value.values()) {
        if (!members.some((m) => paired.has(m))) continue;
        const at = members.map((m) => index.get(m) ?? index.set(m, nodes.push(m) - 1).get(m) ?? -1);
        for (const x of at) for (const y of at) if (x !== y) links.push([x, y, GROUP]);
      }
      const n = nodes.length;
      // Nodes without links are dropped from the walk (their bumps are not followed).
      links.sort((p, q) => p[0] - q[0]);
      start = new Int32Array(n + 1);
      for (const [x] of links) start[x + 1]++;
      for (let i = 0; i < n; i++) start[i + 1] += start[i];
      other = Int32Array.from(links, (l) => l[1]);
      kind = Int32Array.from(links, (l) => l[2]);
      seen = new Float64Array(3 * n).fill(Number.NaN);
      kept = new Float64Array(3 * n).fill(-1);
      led = new Float64Array(n);
      root = new Uint8Array(n);
      mark = new Int32Array(n);
      queue = new Int32Array(n);
      for (let i = 0; i < n; i++) {
        const j = before.get(nodes[i]);
        if (j !== undefined) {
          seen.set(oldSeen.subarray(3 * j, 3 * j + 3), 3 * i);
          led[i] = oldLed[j];
        }
        if (nodes[i] && focused.has(nodes[i] as SectionViewport)) led[i] = ++clock;
        root[i] = sec && i < sec.count && typeof sec.params[i]?.follow !== "string" ? 1 : 0;
      }
      focused.clear();
      busy = true;
      try {
        for (let i = 0; i < n; i++) {
          const vp = nodes[i];
          if (vp && start[i + 1] > start[i]) subs.push(vp.live.version.subscribe(onVersion(i)));
        }
      } finally {
        busy = false;
      }
      // Whatever the graph moved (a new axis, data loaded, a route that now exists) is set from its authority.
      const done = new Uint8Array(n);
      for (let i = 0; i < n; i++) {
        if (done[i] || start[i + 1] === start[i]) continue;
        const count = component(i);
        for (let h = 0; h < count; h++) done[queue[h]] = 1;
        resync(i);
      }
    };
    // A section's own focus change from a patch (origin `graph`) moves it and leads, as a snapshot's does
    // (seqspec:interaction#replays-and-carries): it is the last to lead before its viewport emits.
    add(
      ctx.focus.changed.subscribe((c) => {
        const vp = c.origin === "graph" && c.entry !== null ? ctx.viewports.get(c.section) : null;
        if (!vp) return;
        const n = nodes.indexOf(vp);
        if (n >= 0) led[n] = ++clock;
        else focused.add(vp);
      }),
    );
    add(ctx.graph.assembled.subscribe(rebuild));
    add(ctx.scrollGroups.groups.subscribe(rebuild));
    add(release);
  },
);
