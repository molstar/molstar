/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */
// @module-tag heavy

// Against immutable.js (README §8): get, set, delete, iterate and build at 10, 1,000 and 100,000 entries, keys
// including 22-character string ids, as ratios (ours over immutable.js; the target is at most 1). A miss is printed
// with its case; the test fails only past `Catastrophic`, a regression no machine's noise explains. And the flat
// threshold's measurement (README §2): a lookup and a persistent set on a flat array and on a trie of 1 to 32 entries.

import { List, Map as IMap, OrderedMap, Set as ISet } from 'immutable';
import { describe, expect, it } from 'vitest';
import { ChampNode, FlatNode, type Root, rootGet, rootSet } from '../champ.js';
import { DefaultConfig } from '../config.js';
import { HashMap } from '../hash-map.js';
import { HashSet } from '../hash-set.js';
import { OrderedHashMap } from '../ordered-hash-map.js';
import { Vector } from '../vector.js';
import { time } from './bench.js';
import { id22, Rng } from './model.js';

const Catastrophic = 2.5;
const Sizes = [10, 1_000, 100_000];

interface Row {
  name: string;
  ours: number;
  theirs: number;
}

function report(title: string, rows: Row[]) {
  const lines = rows.map(
    (r) =>
      `${r.name.padEnd(44)} ${r.ours.toFixed(1).padStart(9)} ns ${r.theirs.toFixed(1).padStart(9)} ns  ${(r.ours / r.theirs).toFixed(2)}${r.ours > r.theirs ? '  MISS' : ''}`,
  );
  console.info(
    `${title}\n${'case'.padEnd(44)} ${'ours'.padStart(12)} ${'immutable'.padStart(12)}  ratio\n${lines.join('\n')}`,
  );
  for (const r of rows) expect(r.ours / r.theirs, r.name).toBeLessThan(Catastrophic);
}

function keysOf(kind: 'id22' | 'int', n: number, seed: number): unknown[] {
  const rng = new Rng(seed);
  if (kind === 'int') return rng.shuffle(Array.from({ length: n }, (_, i) => i * 7 + 1));
  const seen = new Set<string>();
  while (seen.size < n) seen.add(id22(rng));
  return [...seen];
}

type MapLike = {
  get(k: unknown): unknown;
  set(k: unknown, v: number): MapLike;
  delete(k: unknown): MapLike;
  forEach(fn: (v: number, k: unknown) => unknown): unknown;
};

function mapCases(
  label: string,
  build: (entries: [unknown, number][]) => MapLike,
  ibuild: (entries: [unknown, number][]) => MapLike,
): Row[] {
  const rows: Row[] = [];
  for (const kind of ['id22', 'int'] as const) {
    for (const n of Sizes) {
      const keys = keysOf(kind, n, n);
      const present = new Set(keys);
      const fresh = (kind === 'int' ? keysOf(kind, n, n).map((x) => (x as number) + 3) : keysOf(kind, n, n + 1)).filter(
        (k) => !present.has(k),
      );
      const entries = keys.map((k, i) => [k, i] as [unknown, number]);
      const ours = build(entries);
      const theirs = ibuild(entries);
      const k = Math.min(n, 1000);
      const name = (op: string) => `${label} ${op} ${kind} n=${n}`;
      const pair = (op: string, f: (m: MapLike) => unknown, ops: number, g: (m: MapLike) => unknown) =>
        rows.push({ name: name(op), ours: time(() => f(ours), ops), theirs: time(() => g(theirs), ops) });

      const get = (m: MapLike) => {
        let s = 0;
        for (let i = 0; i < keys.length; i++) s += m.get(keys[i]) as number;
        return s;
      };
      pair('get', get, n, get);
      const update = (m: MapLike) => {
        for (let i = 0; i < k; i++) m = m.set(keys[i], -i);
        return m;
      };
      pair('set (update)', update, k, update);
      const insert = (m: MapLike) => {
        for (let i = 0; i < k && i < fresh.length; i++) m = m.set(fresh[i], i);
        return m;
      };
      pair('set (insert)', insert, Math.min(k, fresh.length), insert);
      const del = (m: MapLike) => {
        for (let i = 0; i < k; i++) m = m.delete(keys[i]);
        return m;
      };
      pair('delete', del, k, del);
      const iterate = (m: MapLike) => {
        let s = 0;
        m.forEach((v) => {
          s += v;
        });
        return s;
      };
      pair('iterate', iterate, n, iterate);
      rows.push({ name: name('build'), ours: time(() => build(entries), n), theirs: time(() => ibuild(entries), n) });
    }
  }
  return rows;
}

describe('against immutable.js', () => {
  it('HashMap and Map', () => {
    const rows = mapCases(
      'HashMap',
      (e) => HashMap.of(e) as unknown as MapLike,
      (e) =>
        IMap<unknown, number>().withMutations((m) => {
          for (const [k, v] of e) m.set(k, v);
        }) as unknown as MapLike,
    );
    report('HashMap / immutable Map', rows);
  });

  it('OrderedHashMap and OrderedMap', () => {
    const rows = mapCases(
      'OrderedHashMap',
      (e) => OrderedHashMap.of(e) as unknown as MapLike,
      (e) =>
        OrderedMap<unknown, number>().withMutations((m) => {
          for (const [k, v] of e) m.set(k, v);
        }) as unknown as MapLike,
    );
    report('OrderedHashMap / immutable OrderedMap', rows);
  });

  it('HashSet and Set', () => {
    const rows: Row[] = [];
    for (const n of Sizes) {
      const keys = keysOf('id22', n, n);
      const ours = HashSet.of(keys);
      const theirs = ISet(keys);
      const k = Math.min(n, 1000);
      rows.push({
        name: `HashSet has id22 n=${n}`,
        ours: time(() => keys.reduce((s: number, x) => s + (ours.has(x) ? 1 : 0), 0), n),
        theirs: time(() => keys.reduce((s: number, x) => s + (theirs.has(x) ? 1 : 0), 0), n),
      });
      rows.push({
        name: `HashSet delete id22 n=${n}`,
        ours: time(() => {
          let s = ours;
          for (let i = 0; i < k; i++) s = s.delete(keys[i]);
          return s;
        }, k),
        theirs: time(() => {
          let s = theirs;
          for (let i = 0; i < k; i++) s = s.delete(keys[i]);
          return s;
        }, k),
      });
      rows.push({
        name: `HashSet build id22 n=${n}`,
        ours: time(() => HashSet.of(keys), n),
        theirs: time(() => ISet(keys), n),
      });
    }
    report('HashSet / immutable Set', rows);
  });

  it('Vector and List', () => {
    const rows: Row[] = [];
    for (const n of Sizes) {
      const xs = Array.from({ length: n }, (_, i) => i);
      const ours = Vector.of(xs);
      const theirs = List(xs);
      const k = Math.min(n, 1000);
      rows.push({
        name: `Vector get n=${n}`,
        ours: time(() => {
          let s = 0;
          for (let i = 0; i < n; i++) s += ours.get(i)!;
          return s;
        }, n),
        theirs: time(() => {
          let s = 0;
          for (let i = 0; i < n; i++) s += theirs.get(i)!;
          return s;
        }, n),
      });
      rows.push({
        name: `Vector set n=${n}`,
        ours: time(() => {
          let v = ours;
          for (let i = 0; i < k; i++) v = v.set(i, -i);
          return v;
        }, k),
        theirs: time(() => {
          let v = theirs;
          for (let i = 0; i < k; i++) v = v.set(i, -i);
          return v;
        }, k),
      });
      rows.push({
        name: `Vector push n=${n}`,
        ours: time(() => {
          let v = ours;
          for (let i = 0; i < k; i++) v = v.push(i);
          return v;
        }, k),
        theirs: time(() => {
          let v = theirs;
          for (let i = 0; i < k; i++) v = v.push(i);
          return v;
        }, k),
      });
      rows.push({
        name: `Vector pop n=${n}`,
        ours: time(() => {
          let v = ours;
          for (let i = 0; i < k; i++) v = v.pop();
          return v;
        }, k),
        theirs: time(() => {
          let v = theirs;
          for (let i = 0; i < k; i++) v = v.pop();
          return v;
        }, k),
      });
      rows.push({
        name: `Vector iterate n=${n}`,
        ours: time(() => {
          let s = 0;
          ours.forEach((x) => {
            s += x;
          });
          return s;
        }, n),
        theirs: time(() => {
          let s = 0;
          theirs.forEach((x) => {
            s += x;
          });
          return s;
        }, n),
      });
      rows.push({
        name: `Vector build n=${n}`,
        ours: time(() => Vector.of(xs), n),
        theirs: time(() => List(xs), n),
      });
    }
    report('Vector / immutable List', rows);
  });
});

describe('the flat threshold (README §2)', () => {
  it('a lookup and a persistent set, flat and as a trie, 1 to 32 entries', () => {
    const lines: string[] = [];
    for (const kind of ['id22', 'int'] as const) {
      for (let n = 1; n <= 32; n++) {
        const keys = keysOf(kind, n, 100 + n);
        const entries: unknown[] = [];
        for (let i = 0; i < n; i++) entries.push(keys[i], i);
        const flat: Root = new FlatNode(null, entries);
        let trie: Root = new ChampNode({ dead: false }, 0, 0, []);
        for (let i = 0; i < n; i++) trie = rootSet(trie, i, keys[i], i, DefaultConfig, null);
        const imap = IMap<unknown, number>(keys.map((k, i) => [k, i]));
        const get = (r: Root) => {
          let s = 0;
          for (let i = 0; i < n; i++) s += rootGet(r, keys[i], DefaultConfig);
          return s;
        };
        const set = (r: Root) => {
          for (let i = 0; i < n; i++) rootSet(r, n, keys[i], -1, DefaultConfig, null);
          return r;
        };
        const fg = time(() => get(flat), n, 10, 5);
        const tg = time(() => get(trie), n, 10, 5);
        const ig = time(
          () => {
            let s = 0;
            for (let i = 0; i < n; i++) s += imap.get(keys[i])!;
            return s;
          },
          n,
          10,
          5,
        );
        const fs = time(() => set(flat), n, 10, 5);
        const ts = time(() => set(trie), n, 10, 5);
        const is = time(
          () => {
            for (let i = 0; i < n; i++) imap.set(keys[i], -1);
            return imap;
          },
          n,
          10,
          5,
        );
        lines.push(
          `${kind.padEnd(5)} n=${String(n).padStart(2)}  get flat ${fg.toFixed(1).padStart(6)} trie ${tg.toFixed(1).padStart(6)} immutable ${ig.toFixed(1).padStart(6)}   set flat ${fs.toFixed(1).padStart(6)} trie ${ts.toFixed(1).padStart(6)} immutable ${is.toFixed(1).padStart(6)}  (ns)`,
        );
      }
    }
    console.info(`flat threshold\n${lines.join('\n')}`);
    expect(lines.length).toBe(64);
  });
});
