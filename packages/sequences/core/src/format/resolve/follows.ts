/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Axis links (seqspec:format#sections; LK1): `ResolvedView.follows`, each `follow` pair once, with a route per block of
// either side: from that block's coordinate system through its section's anchor to a block of the other section. A
// blocks axis anchors on its own blocks (the block under the center leads); an alignment axis on `follow_anchor`, else
// its reference row (columns reach other systems only through rows). Checked here: a section of another view, itself, a
// loop (through follows and scroll groups), and two sections no path joins are `svs.section.follow`, and that pair is
// left out. A pair whose side has no blocks yet, or whose anchor row is not held (a lazy alignment), is kept without
// routes: the viewer leaves the follower until it loads. Resolve-time code (edges): it allocates.

import {
  type ComposedPath,
  type CoordinateSystem,
  CoordinateSystemKind,
  Coords,
  Issue,
  KeyKind,
  type KeyTarget,
  Mappings,
  type PathStep,
  type SeqSet,
} from "@molstar/seq-core/model";
import { SectionAxis } from "./codes.ts";
import { blockSystem } from "./reach.ts";
import type { ResolvedSection } from "./section.ts";

/** One row per pair (`a` followed by `b`), and per pair a route per block of `a`, then per block of `b`. */
export interface ResolvedFollows {
  readonly count: number;
  /** Sections, `b` the one whose `follow` names `a`; each leads when it moves. */
  readonly a: Int32Array;
  readonly b: Int32Array;
  /** b's cell over a's through the first route: 1, or 1/3 and 3 through a codon step. */
  readonly scale: Float64Array;
  /** Into `routes`: a's blocks at `routeStart`, then b's (`blockCount` of a later). */
  readonly routeStart: Int32Array;
  readonly routes: {
    readonly count: number;
    /** Into `paths`; -1: none from this block (yet). */
    readonly path: Int32Array;
    /** The block of the other section it reaches (an index within that section); -1: none. */
    readonly to: Int32Array;
    /** 1: the path's map reads from this block's system (`map.from`); 0: from the other's (to → from). */
    readonly forward: Uint8Array;
  };
}

const t: KeyTarget = { kind: -1, id: -1, row: -1 };
const part: PathStep[] = [];
const PENDING = null;
type Anchor = CoordinateSystem | typeof PENDING;

/** Appends the steps from → to; false when no single shortest path joins them. */
function append(set: SeqSet, from: CoordinateSystem, to: CoordinateSystem, out: PathStep[]): boolean {
  const n = Mappings.path(set, from, to, part);
  for (let i = 0; n > 0 && i < part.length; i++) out.push(part[i]);
  return n >= 0;
}

/** The pairs of one view (`composeView`), their paths through `pathOf`; issues for what cannot be followed. */
export function followsOf(
  set: SeqSet,
  view: string,
  sections: readonly ResolvedSection[],
  pathOf: (p: ComposedPath) => number,
  issues: Issue[],
): ResolvedFollows {
  const index = new Map(sections.map((s, i) => [s.ref, i]));
  // Sections already linked, by follows and scroll groups: a follow within one is a loop.
  const up = Int32Array.from(sections.keys());
  const find = (i: number): number => {
    while (up[i] !== i) i = up[i];
    return i;
  };
  const groups = new Map<number, number>();
  for (let i = 0; i < sections.length; i++) {
    const g = sections[i].row.scrollGroup;
    if (g < 0) continue;
    if (!groups.has(g)) groups.set(g, i);
    up[find(i)] = find(groups.get(g) ?? i);
  }
  const say = (s: ResolvedSection, param: string, message: string) =>
    issues.push(Issue.error("svs.section.follow", `${s.path}/params/${param}`, message));
  const anchors = sections.map((s) => anchorOf(set, s, say));
  const out = { a: [] as number[], b: [] as number[], scale: [] as number[], routeStart: [] as number[] };
  const routes = { path: [] as number[], to: [] as number[], forward: [] as number[] };
  for (let bi = 0; bi < sections.length; bi++) {
    const s = sections[bi];
    const ref = s.row.params.follow;
    if (typeof ref !== "string") continue;
    const ai = index.get(ref) ?? -1;
    if (ai < 0) {
      say(s, "follow", `"${ref}" is not a section of view "${view}": a section follows one of its own view`);
      continue;
    }
    if (ai === bi) {
      say(s, "follow", "a section cannot follow itself: name another section of this view");
      continue;
    }
    if (find(ai) === find(bi)) {
      const why = `"${s.ref}" and "${ref}" are already linked (by follow or a scroll group); a loop has no leader`;
      say(s, "follow", `following "${ref}" closes a loop: ${why}`);
      continue;
    }
    up[find(bi)] = find(ai);
    const start = routes.path.length;
    let found = 0;
    let scale = 1;
    for (const [x, y] of [
      [ai, bi],
      [bi, ai],
    ]) {
      for (let k = 0; k < sections[x].blocks.length; k++) {
        const r = route(set, sections[x], k, anchors[x], sections[y]);
        const forward = r !== null && Coords.sameSystem(r.path.map.from, r.path.from);
        routes.path.push(r ? pathOf(r.path) : -1);
        routes.to.push(r ? r.to : -1);
        routes.forward.push(forward ? 1 : 0);
        // b's cell over a's: a residue's cell holds its codon's three nucleotides.
        if (r && found++ === 0) scale = forward === (x === ai) ? 1 / r.path.map.span : r.path.map.span;
      }
    }
    const pending = [ai, bi].some((i) => sections[i].blocks.length === 0 || anchors[i] === PENDING);
    if (found === 0 && !pending) {
      routes.path.length = start;
      routes.to.length = start;
      routes.forward.length = start;
      say(s, "follow", `no mapping path joins "${s.ref}" and "${ref}": map what one shows to what the other does`);
      continue;
    }
    out.a.push(ai);
    out.b.push(bi);
    out.routeStart.push(start);
    out.scale.push(scale);
  }
  const I = Int32Array;
  return {
    count: out.a.length,
    a: I.from(out.a),
    b: I.from(out.b),
    scale: Float64Array.from(out.scale),
    routeStart: I.from(out.routeStart),
    routes: {
      count: routes.path.length,
      path: I.from(routes.path),
      to: I.from(routes.to),
      forward: Uint8Array.from(routes.forward),
    },
  };
}

/**
 * A section's anchor: on an alignment axis, `follow_anchor` (a row of its alignment, or a unit) or the reference row;
 * PENDING while the axis has no blocks or the row is not held; on a blocks axis the block itself (undefined).
 */
function anchorOf(
  set: SeqSet,
  s: ResolvedSection,
  say: (s: ResolvedSection, param: string, message: string) => void,
): Anchor | undefined {
  const ref = s.row.params.follow_anchor;
  const block = s.blocks[0];
  if (!block) return PENDING;
  if (s.row.axis !== SectionAxis.columns) {
    if (typeof ref === "string")
      say(
        s,
        "follow_anchor",
        "follow_anchor anchors an alignment axis; a blocks axis follows the block under its center",
      );
    return undefined;
  }
  const aln = set.alignments[block.alignment];
  if (!aln) return PENDING;
  let a: CoordinateSystem = { kind: CoordinateSystemKind.row, id: aln.id, row: aln.referenceRow };
  if (typeof ref === "string") {
    if (!set.byName(ref, t)) return PENDING;
    if (t.kind === KeyKind.unit) a = { kind: CoordinateSystemKind.unit, id: t.id, row: -1 };
    else if (t.kind === KeyKind.row && t.id === aln.id) a = { kind: CoordinateSystemKind.row, id: t.id, row: t.row };
    else {
      say(s, "follow_anchor", `"${ref}" is not a row of this section's alignment, nor a unit`);
      return PENDING;
    }
  }
  return Mappings.systemLength(set, a) < 0 ? PENDING : a;
}

/**
 * The first block of y that block k of x reaches through x's anchor (each side's anchor → the other's block system; the
 * other's anchor is its own when it leads), and the path; null: none.
 */
function route(
  set: SeqSet,
  x: ResolvedSection,
  k: number,
  ax: Anchor | undefined,
  y: ResolvedSection,
): { path: ComposedPath; to: number } | null {
  const from = blockSystem(x.blocks[k]);
  const a = ax ?? from;
  for (let l = 0; l < y.blocks.length; l++) {
    const to = blockSystem(y.blocks[l]);
    const steps: PathStep[] = [];
    // One coordinate system on both axes (an alignment's columns twice, a unit and its reversed piece): the identity.
    const same = Coords.sameSystem(from, to);
    if (!same && (ax === PENDING || !append(set, from, a, steps) || !append(set, a, to, steps))) continue;
    try {
      const path = Mappings.compose(set, from, steps);
      if (path.map.count > 0) return { path, to: l };
    } catch {
      // A route through a row not held, or two codon steps: none.
    }
  }
  return null;
}
