/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Map builders from common identities (seqspec:libraries#translate): a host's system against Seq*'s, each one
// `CoordinateMap` (span 1), built at the edge.

import { type CoordinateSystem, Coords, type NameId, NO_LABEL, Numbering } from "@molstar/seq-core/model";

/** Pairs from element `f` with to position `t`, extending the last segment when both continue it. */
function pair(quads: number[], f: number, t: number) {
  const k = quads.length - 4;
  if (k >= 0 && quads[k] + quads[k + 1] === f && quads[k + 2] + quads[k + 1] === t) quads[k + 1]++;
  else quads.push(f, 1, t, 1);
}

/** Label → the first position with it; NO_LABEL is none. */
export function firstOf(labels: Int32Array): Map<number, number> {
  const first = new Map<number, number>();
  for (let i = labels.length - 1; i >= 0; i--) if (labels[i] !== NO_LABEL) first.set(labels[i], i);
  return first;
}

/** Positions with equal labels pair (`label_seq_id`): each from position with the first to position of its label. */
function byLabels(from: CoordinateSystem, fromLabels: Int32Array, to: CoordinateSystem, toLabels: Int32Array) {
  const first = firstOf(toLabels);
  const quads: number[] = [];
  for (let i = 0; i < fromLabels.length; i++) {
    const j = first.get(fromLabels[i]);
    if (j !== undefined) pair(quads, i, j);
  }
  return Coords.create(from, to, 1, quads);
}

/**
 * `to` positions labeled in `from`'s numbering (auth numbers): each pairs with the from index its label names. `codes`
 * gives insertion codes, matched against a `labels` numbering's suffixes; a code it does not have pairs nothing.
 */
function byNumbering(
  from: CoordinateSystem,
  n: Numbering,
  to: CoordinateSystem,
  labels: Int32Array,
  codes?: readonly (string | null)[],
) {
  const quads: number[] = [];
  for (let j = 0; j < labels.length; j++) {
    const code = codes?.[j];
    const suffix = code ? (n.method === "labels" ? n.strings.find(code) : -1) : -1;
    if (code && suffix < 0) continue;
    const i = Numbering.indexOf(n, labels[j], suffix as NameId);
    if (i >= 0) pair(quads, i, j);
  }
  return Coords.create(from, to, 1, quads);
}

/** A SIFTS-style table: (fromStart, toStart, length) triples, 0-based; segments may leave holes on either side. */
function bySegments(from: CoordinateSystem, to: CoordinateSystem, segments: Int32Array) {
  const quads: number[] = [];
  for (let k = 0; k + 2 < segments.length; k += 3) quads.push(segments[k], segments[k + 2], segments[k + 1], 1);
  return Coords.create(from, to, 1, quads);
}

export const TranslateMaps = { byLabels, byNumbering, bySegments };
