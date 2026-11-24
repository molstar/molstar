/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A `seq` source's dynamic nodes found by key (_spec/README.md "Data", The seq format: its nodes found by key): a lazy
// chunk of the runtime, loaded on the first find (most views never ask).

import { abortError, type FindNamesOptions, type NameHit } from './protocol.js';
import type { SeqData } from './seq-format.js';
import { SeqNodeKeyField } from './seq-runtime.js';

/** Each data's index, by lowercased key: built on the first find. */
const Indexes = new WeakMap<SeqData, Map<string, NameHit>>();

/** The index: a key's hit is the extent of the pairs whose answers bring it, on the first system that does. */
function indexOf(data: SeqData): Map<string, NameHit> {
  let m = Indexes.get(data);
  if (m) return m;
  m = new Map();
  for (const { request: r, answer } of data.pairs) {
    if (!answer.ok) continue;
    for (const n of answer.nodes ?? []) {
      const k = n.key.toLowerCase();
      const h = m.get(k);
      if (!h)
        m.set(k, { system: r.system, start: r.start, end: r.end, id: n.key, label: n.key, field: SeqNodeKeyField });
      else if (h.system === r.system)
        m.set(k, { ...h, start: Math.min(h.start, r.start), end: Math.max(h.end, r.end) });
    }
  }
  Indexes.set(data, m);
  return m;
}

/** The nodes of a key (exact, case-insensitive), then by prefix when asked, at most `max` (20). */
export function findNodeKeys(
  data: SeqData,
  text: string,
  options: FindNamesOptions,
  signal: AbortSignal,
): Promise<readonly NameHit[]> {
  if (signal.aborted) return Promise.reject(abortError());
  const index = indexOf(data);
  const q = text.trim().toLowerCase();
  const max = options.max ?? 20;
  const exact = index.get(q);
  const out: NameHit[] = exact ? [exact] : [];
  if (options.prefix)
    for (const [k, h] of index) {
      if (out.length >= max) break;
      if (k !== q && k.startsWith(q)) out.push(h);
    }
  return Promise.resolve(out.slice(0, max));
}
