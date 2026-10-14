/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Item coordinate systems (ST7, seqspec:sources#data-patches): a streamed items patch whose rows carry a `map` (a
// `codon` or `segments` mapping's params from the item's own positions onto the source's system) gives, per item, two
// dynamic nodes of its source, keyed by the item's `id`: its coordinate system (`#item_system`: a letterless molecule
// and its own unit, `<source ref>/<id>`) and its mapping (`<source ref>/<id>.map`). They come and go with the regions
// that carry them, as a patch's own nodes do (dynamic.ts), but are read into the set only when something names the item
// (`#assemble.stream`): until then they are its row's `map`, held in the cache. Arrival-rate.

import { ItemSystemKind, type StreamNode, type StreamTable } from '@molstar/seq-core/format';

/**
 * An item's own length from its map's params, or why none can be read: a codon map's residues (its `from_range`'s
 * end, else its CDS, a multiple of 3, over 3; each range `[start, stop)` with start before stop, on either strand, as
 * the model's codon map reads it), a segments map's furthest from position (every segment's length given).
 */
function lengthOf(p: Readonly<Record<string, unknown>>): number | string {
  const rows = (x: unknown): unknown[][] => (Array.isArray(x) ? x.filter((r): r is unknown[] => Array.isArray(r)) : []);
  const int = (v: unknown) => typeof v === 'number' && Number.isInteger(v);
  if (p.method === 'codon') {
    const range = p.from_range;
    if (Array.isArray(range) && range.length === 2 && int(range[1])) return range[1] as number;
    if (!Array.isArray(p.cds)) return 'a codon map needs cds (or from_range)';
    let n = 0;
    for (const [a, e] of rows(p.cds)) {
      if (!int(a) || !int(e)) continue;
      if ((e as number) <= (a as number))
        return `its CDS range [${a}, ${e}) is empty or reversed (write [start, stop))`;
      n += (e as number) - (a as number);
    }
    return n % 3 === 0
      ? n / 3
      : `its CDS is ${n} bases, not a multiple of 3 (write a partial CDS from its first whole codon)`;
  }
  if (p.method !== 'segments') return `a map is a codon or segments mapping's params, not "${String(p.method)}"`;
  let n = 0;
  for (const [f, , len] of rows(p.segments)) {
    if (!int(f) || !int(len)) return "a segments map needs each segment's length";
    n = Math.max(n, (f as number) + (len as number));
  }
  return n > 0 ? n : 'a segments map needs at least one segment';
}

/**
 * The item coordinate systems and mappings of items table `t` from `source`, over its system `to` (a ref): one pair per
 * item `id` whose rows carry a `map` (the first row's), into `out`. An item whose map gives no length (or a CDS not a
 * multiple of 3) is reported (`svs.stream.item_map`) and left out; its mapping's other faults are its node's checks'.
 */
export function itemSystemNodes(
  source: string,
  to: string,
  t: StreamTable,
  out: StreamNode[],
  report: (code: string, message: string) => void,
) {
  const map = t.map;
  const id = t.id;
  if (!map || !id) return;
  const seen = new Set<string>();
  for (let j = 0; j < map.length && j < id.length; j++) {
    const m = map[j];
    const key = id[j] == null ? '' : String(id[j]);
    if (m == null || typeof m !== 'object' || Array.isArray(m) || key === '' || seen.has(key)) continue;
    seen.add(key);
    const p = m as Readonly<Record<string, unknown>>;
    const length = lengthOf(p);
    if (typeof length === 'string') {
      report('svs.stream.item_map', `${source}: item "${key}"'s map is left out: ${length}`);
      continue;
    }
    const label = t.label && j < t.label.length && t.label[j] != null ? String(t.label[j]) : key;
    const alphabet = p.method === 'codon' ? 'protein' : 'generic';
    out.push({ kind: ItemSystemKind as StreamNode['kind'], key, params: { length, alphabet, label } });
    out.push({ kind: 'mapping', key: `${key}.map`, params: { ...p, from: `${source}/${key}`, to } });
  }
}
