/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// What a search reads (_spec/README.md "Search"): every system a letters layer of the resolved view draws, whole, in
// the view's order (sections top to bottom, a section's blocks in order), or what `within` names (its parts the targets,
// their ranges the windows; a system no section shows only then). Letters through the store: a region asked by a want
// of the search's own and waited for, never read as empty while pending (EC-169); a failed read is the target's issue.

import type { Location } from '@molstar/seqx-core/query';
import type { DataEntry, DataStore, LettersTable, Want } from '@molstar/seqx-core/data';
import { LettersType } from '@molstar/seqx-core/render';
import type { StaticModel } from '@molstar/seqx-core/resolve';
import type { LetterReader, SearchTarget } from '@molstar/seqx-core/search';
import type { ResolvedLayer, ResolvedView } from '@molstar/seqx-core/view';

/** A target and where its letters are read from. */
export interface SearchTargetEntry extends SearchTarget {
  /** The system's words (its label, else its id). */
  readonly label: string;
}

/** A letters table's code at a position (a coded table through its alphabet), or 0. */
function codeAt(t: LettersTable, p: number): number {
  const i = p - t.start;
  if (i < 0 || i >= t.letters.length) return 0;
  if (typeof t.letters === 'string') return t.letters.charCodeAt(i);
  return t.alphabet ? t.alphabet.charCodeAt(t.letters[i]) || 0 : 0;
}

/** Thrown for a region whose read failed (its target skipped). */
class ReadFailed extends Error {}

/**
 * A system's letters read through the store by region: wants of the reader's own, as many tiles as the store asks for
 * at once, each waited for until settled; a pending tile is waited for, never read as empty (EC-169).
 */
export function storeLetters(store: DataStore, source: string, system: string): LetterReader {
  const step = Math.max(1024, store.settings.tileBins * Math.max(1, store.settings.maxTilesPerWant >> 1));
  const fill = (out: Uint8Array, start: number, want: Want): boolean => {
    if (store.unopened(source)) throw new ReadFailed('its source could not be opened');
    // a source not open yet answers nothing a search can wait on: its target is skipped, said, never read as empty
    if (store.coverage(source) === undefined) throw new ReadFailed('its source is not open yet');
    const keys = store.keysOf(want);
    if (!keys.length && store.tooMany(want)) throw new ReadFailed('past the store’s bound on a want');
    const tables: LettersTable[] = [];
    for (const k of keys) {
      const e: DataEntry | undefined = store.get(k);
      if (e?.status === 'failed') throw new ReadFailed('a region could not be read');
      if (e?.status === 'loaded') {
        if (e.tables.letters) tables.push(e.tables.letters);
      } else if (e?.status !== 'empty') return false;
    }
    for (let p = want.start; p < want.end; p++) {
      let c = 0;
      for (let k = 0; k < tables.length && !c; k++) c = codeAt(tables[k], p);
      out[p - start] = c;
    }
    return true;
  };
  const settled = (signal: AbortSignal, want: Want): Promise<void> =>
    new Promise((resolve) => {
      const w = store.wanter();
      w.want([want]);
      let done = false;
      const end = () => {
        if (done) return;
        done = true;
        signal.removeEventListener('abort', end);
        w.dispose();
        resolve();
      };
      signal.addEventListener('abort', end);
      void w.whenSettled().then(end);
    });
  return {
    async read(start, end, signal) {
      const out = new Uint8Array(end - start);
      for (let s = start; s < end; s += step) {
        const want: Want = { source, system, start: s, end: Math.min(end, s + step), level: 1, priority: 'view' };
        while (!fill(out, start, want)) {
          if (signal.aborted) return out;
          await settled(signal, want);
        }
      }
      return out;
    },
  };
}

/** The letters layers of a resolved view, in its order: sections top to bottom, a section's blocks in order. */
function lettersLayers(view: ResolvedView): ResolvedLayer[] {
  const out: ResolvedLayer[] = [];
  const of = (layers: readonly ResolvedLayer[]) => {
    for (const l of layers) if (l.type === 'letters' && l.placed) out.push(l);
  };
  for (const s of view.sections) {
    if (s.blocks) for (const b of s.blocks) for (const ls of b.layers) of(ls ?? []);
    else for (const t of s.tracks) of(t.layers);
  }
  return out;
}

/**
 * The targets of a search: every system a letters layer draws, once, in the view's order; with `within`, its parts'
 * systems (those shown first, in the view's order, then the others by address), each with its ranges.
 */
export function searchTargets(
  model: StaticModel,
  view: ResolvedView | undefined,
  store: DataStore,
  within?: Location,
): { readonly targets: SearchTargetEntry[]; readonly unknown: string[] } {
  const order: string[] = [];
  const read = new Map<string, { source: string; system: string }>();
  for (const l of view ? lettersLayers(view) : []) {
    const src = LettersType.prepare(l, model).src;
    if (!src || read.has(src.dataSystem)) continue;
    order.push(src.dataSystem);
    read.set(src.dataSystem, { source: src.source, system: src.system });
  }
  const ranges = new Map<string, Int32Array>();
  if (within) {
    for (const p of within.parts) ranges.set(p.system, Int32Array.from(p.ranges));
    const shown = order.filter((s) => ranges.has(s));
    const others = [...ranges.keys()].filter((s) => !read.has(s)).sort();
    order.length = 0;
    order.push(...shown, ...others);
  }
  const out: SearchTargetEntry[] = [];
  // a system whose length is not known yet (its source opening): said, never searched as empty
  const unknown: string[] = [];
  for (const address of order) {
    const sys = model.systems.get(address);
    const from = read.get(address) ?? (sys?.source ? { source: sys.source, system: sys.name } : undefined);
    if (!sys?.alphabet || !from) continue;
    if (sys.length === undefined) {
      unknown.push(address);
      continue;
    }
    const r = ranges.get(address);
    out.push({
      system: address,
      label: sys.label ?? address.slice(address.lastIndexOf('/') + 1),
      alphabet: sys.alphabet,
      length: sys.length,
      cyclic: sys.cyclic,
      ...(r ? { within: r } : {}),
      letters: storeLetters(store, from.source, from.system),
    });
  }
  return { targets: out, unknown };
}
