/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A `seq` source's dynamic nodes found by key (_spec/README.md "Data", The seq format: its nodes found by key): a lazy
// chunk of the runtime, loaded on the first find (most views never ask).

import { abortError, type FindNamesOptions, type NameHit, SeqNodeKeyField } from './protocol.js';
import type { SeqData } from './seq-format.js';

/** A data's index: each key's hit, and the keys by their case-folded form (for case-insensitive and prefix finds). */
interface Index {
  readonly byKey: Map<string, { hit: NameHit; level: number }>;
  readonly folded: Map<string, string[]>;
}

const Indexes = new WeakMap<SeqData, Index>();

/**
 * The index: a key's hit is the first region whose answer brings it, at the finest level that brings it (reading that
 * one region brings the node; never the span between regions far apart). Keys are exact: two that differ by case are two.
 */
function indexOf(data: SeqData): Index {
  let index = Indexes.get(data);
  if (index) return index;
  const byKey = new Map<string, { hit: NameHit; level: number }>();
  const folded = new Map<string, string[]>();
  for (const { request: r, answer } of data.pairs) {
    if (!answer.ok) continue;
    for (const n of answer.nodes ?? []) {
      const had = byKey.get(n.key);
      const finer = !had || r.level < had.level || (r.level === had.level && r.start < had.hit.start);
      if (!had) {
        const f = n.key.toLowerCase();
        const keys = folded.get(f);
        if (keys) keys.push(n.key);
        else folded.set(f, [n.key]);
      }
      if (finer)
        byKey.set(n.key, {
          hit: { system: r.system, start: r.start, end: r.end, id: n.key, label: n.key, field: SeqNodeKeyField },
          level: r.level,
        });
    }
  }
  index = { byKey, folded };
  Indexes.set(data, index);
  return index;
}

/**
 * The nodes of a key: the exact key first, then the keys equal to it case-insensitively, then by prefix
 * (case-insensitive) when asked, at most `max` (20).
 */
export function findNodeKeys(
  data: SeqData,
  text: string,
  options: FindNamesOptions,
  signal: AbortSignal,
): Promise<readonly NameHit[]> {
  if (signal.aborted) return Promise.reject(abortError());
  const { byKey, folded } = indexOf(data);
  const t = text.trim();
  const q = t.toLowerCase();
  const max = options.max ?? 20;
  const keys: string[] = [];
  if (byKey.has(t)) keys.push(t);
  for (const k of folded.get(q) ?? []) if (k !== t) keys.push(k);
  if (options.prefix) for (const [f, ks] of folded) if (f !== q && f.startsWith(q)) keys.push(...ks);
  return Promise.resolve(keys.slice(0, max).map((k) => byKey.get(k)!.hit));
}
