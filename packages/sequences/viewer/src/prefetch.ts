/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Prefetch hints (seqspec:viewer#snapshots-and-playback, seqspec:sources#the-viewers-side; ST6): the focus regions of
// the snapshots next to the one shown, which the stream manager asks for after the view's and the margin's, so entering
// one draws from the cache. A hint is a place (a target's name and a range of its positions) read from the snapshots'
// ops over the base's focus nodes; the manager matches it to the blocks that show the target, at the level that view
// would need. Edge-rate: made when the document or the snapshot changes, read as a stable list.

import type { DocNode, SeqViewSpecDoc } from '@molstar/seq-core/format';

export interface Hint {
  /** A block's unit by name: `hg38/chr17`, or its first segment(s): `hg38`'s members and `chr17/unit` both name it. */
  readonly target: string;
  readonly start: number;
  readonly stop: number;
}

/** The `focus` nodes of `doc` by ref. */
function foci(doc: SeqViewSpecDoc): Map<string, Record<string, unknown>> {
  const out = new Map<string, Record<string, unknown>>();
  const walk = (n: DocNode) => {
    if (n.kind === 'focus' && n.ref && n.params) out.set(n.ref, n.params);
    for (const c of n.children ?? []) walk(c);
  };
  walk(doc.root);
  return out;
}

/**
 * The hints of the snapshots before and after `current` (null: the base, whose next is the first): each snapshot's ops
 * that set a focus node's `target`, `start` or `stop` over the base's, and for the base itself its focus nodes. A focus
 * without a target and a range is no hint.
 */
export function snapshotHints(doc: SeqViewSpecDoc | null, current: string | null): Hint[] {
  const list = doc?.snapshots ?? [];
  if (!doc || list.length === 0) return [];
  const base = foci(doc);
  const at = current === null ? -1 : list.findIndex((s) => s.ref === current);
  const out: Hint[] = [];
  const add = (params: Record<string, unknown>) => {
    const { target, start, stop } = params;
    if (typeof target === 'string' && typeof start === 'number' && typeof stop === 'number' && stop > start)
      out.push({ target, start, stop });
  };
  for (const k of [at - 1, at + 1]) {
    if (k < -1 || k >= list.length || (k === -1 && at === -1)) continue;
    if (k === -1) {
      for (const p of base.values()) add(p);
      continue;
    }
    const merged = new Map<string, Record<string, unknown>>();
    for (const op of list[k].ops ?? []) {
      if (op.op !== 'set' || !op.ref || !base.has(op.ref)) continue;
      merged.set(op.ref, { ...(merged.get(op.ref) ?? base.get(op.ref)), ...op.params });
    }
    for (const p of merged.values()) add(p);
  }
  return out;
}
