/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A use's `via` (_spec/README.md "Resolve", Paths): the maps a
// section (or a structure binding) names; for that use the search leaves out their rivals, the other maps joining the
// same two systems, and paths that never reach such a pair are unchanged.

import type { Value } from '@molstar/core/repr/tree';
import { exactAddress, type Issue, issueOf } from '@molstar/core/repr/tree';
import { ResolveIssues } from './issues.js';
import type { PathResult } from './model.js';
import type { Edge, PathGraph } from './paths.js';

/** The maps a `via` names, sorted and deduplicated (a key part); empty for none. */
export function viaList(via: Value | readonly string[] | undefined): readonly string[] {
  const list = Array.isArray(via) ? via : via === undefined ? [] : [via];
  const out = new Set<string>();
  for (const v of list) {
    const a = exactAddress(v as Value);
    if (a !== undefined) out.add(a);
  }
  return [...out].sort();
}

/** A path entry's key: its pair, and its `via` when it names any. */
export function pathKey(from: string, to: string, via: readonly string[] = []): string {
  return via.length ? `${from}\u0000${to}\u0000${via.join('\u0000')}` : `${from}\u0000${to}`;
}

const pairKey = (a: string, b: string) => (a < b ? `${a}\u0000${b}` : `${b}\u0000${a}`);

/**
 * The rivals of the named maps: every other map joining the same two systems (either way). A name that is no map
 * has none (`viaPath` names it unused).
 */
export function viaRivals(
  via: readonly string[],
  maps: Iterable<{ readonly address: string; readonly from: string; readonly to: string }>,
): ReadonlySet<string> {
  const named = new Set(via);
  const pairs = new Set<string>();
  const all = [...maps];
  for (const m of all) if (named.has(m.address)) pairs.add(pairKey(m.from, m.to));
  const out = new Set<string>();
  if (!pairs.size) return out;
  for (const m of all) if (!named.has(m.address) && pairs.has(pairKey(m.from, m.to))) out.add(m.address);
  return out;
}

/** The graph without the steps of the rivals (a unit's steps are never a rival's). */
export function withoutRivals(graph: PathGraph, rivals: ReadonlySet<string>): PathGraph {
  if (!rivals.size) return graph;
  const kept = (e: Edge) => e.unit !== undefined || !rivals.has(e.map);
  return { edgesFrom: (s) => graph.edgesFrom(s).filter(kept) };
}

/**
 * A path for a use with `via` (the author, 2026-11-20: "Warn"): the search without the named maps' rivals; when the
 * path it finds steps through none of the named maps, the via settled nothing on its pair, so the path is chosen as if
 * it were absent. `unused` lists the named maps the chosen path does not step through (a name that is no map too):
 * what the use reports as `seq.resolve.via.unused`. A search that finds no path keeps its issue and names none.
 */
export function viaPath(named: readonly string[], search: (via: readonly string[]) => PathResult): PathResult {
  if (!named.length) return search(named);
  const r = search(named);
  if (!r.path) return r;
  const on = new Set(r.path.edges.map((e) => e.map));
  const unused = named.filter((m) => !on.has(m));
  if (!unused.length) return r;
  const chosen = unused.length === named.length ? search([]) : r;
  return { ...chosen, unused };
}

/** A use's `via` map that joins no pair on its paths: a warning at the use (a section), naming the map. */
export function viaUnusedIssue(use: string, map: string): Issue {
  const code = 'seq.resolve.via.unused';
  return {
    ...issueOf(ResolveIssues, code, `The via names ${map}, which joins no pair on the paths`, {
      address: use,
      path: ['params', 'via'],
      hint: ResolveIssues.codes[code].hint,
    }),
    places: [{ address: map }],
  };
}
