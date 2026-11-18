/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The hash trie behind `HashMap` and `HashSet` (README §2): one flat array up to `FlatThreshold` entries, a CHAMP trie
// above it (Steindorfer and Vinju, OOPSLA 2015), collision nodes at the bottom. Keys and values interleave in one
// array (`[k0, v0, k1, v1, …]`); a CHAMP node keeps its inline entries from the front in bit order and its sub-nodes
// from the back. Edits copy the path, or write in place a node made under the edit's live owner (README §3).

import { hashEntry, type Owner, type Resolved } from './config.js';

/**
 * The flat representation's largest size (README §2): measured against immutable.js at 1 to 16 entries of
 * 22-character string keys and small integers (`_test/threshold.bench.test.ts`).
 */
export const FlatThreshold = 8;

export const KindFlat = 0;
export const KindChamp = 1;
export const KindCollision = 2;

/** Absent, inside the trie: `undefined` is a value a map may store. */
export const NotFound: unique symbol = Symbol('NotFound') as any;

export class FlatNode {
  readonly kind = KindFlat;
  constructor(
    public owner: Owner | null,
    /** `[k0, v0, k1, v1, …]` in insertion order */
    public entries: any[],
  ) {}
}

export class ChampNode {
  readonly kind = KindChamp;
  constructor(
    public owner: Owner | null,
    public dataMap: number,
    public nodeMap: number,
    /** inline entries `[k, v]` from the front in bit order, sub-nodes from the back in bit order */
    public content: any[],
  ) {}
}

export class CollisionNode {
  readonly kind = KindCollision;
  constructor(
    public owner: Owner | null,
    public hash: number,
    /** `[k0, v0, …]`, unordered */
    public entries: any[],
  ) {}
}

export type TrieNode = ChampNode | CollisionNode;
export type Root = FlatNode | ChampNode;

export const EmptyFlat = new FlatNode(null, []);

/** What the last edit did: read right after the call (an edit's `hash` and `equals` must not edit collections). */
export const Edit = { added: false, removed: false, changed: false };

function resetEdit() {
  Edit.added = false;
  Edit.removed = false;
  Edit.changed = false;
}

function editable(node: { owner: Owner | null }, owner: Owner | null) {
  return owner !== null && node.owner === owner;
}

export function bitCount(i: number): number {
  i = i - ((i >>> 1) & 0x55555555);
  i = (i & 0x33333333) + ((i >>> 2) & 0x33333333);
  i = (i + (i >>> 4)) & 0x0f0f0f0f;
  return Math.imul(i, 0x01010101) >>> 24;
}

function indexOf(bitmap: number, bit: number) {
  return bitCount(bitmap & (bit - 1));
}

function keyEquals(cfg: Resolved, a: any, b: any) {
  return cfg.defaultEquals ? a === b || (Number.isNaN(a) && Number.isNaN(b)) : cfg.equals(a, b);
}

// --- reads ---

function flatGet(entries: any[], key: any, cfg: Resolved): any {
  const n = entries.length;
  if (cfg.defaultEquals) {
    if (!Number.isNaN(key)) {
      for (let i = 0; i < n; i += 2) if (entries[i] === key) return entries[i + 1];
    } else {
      for (let i = 0; i < n; i += 2) if (Number.isNaN(entries[i])) return entries[i + 1];
    }
    return NotFound;
  }
  const eq = cfg.equals;
  for (let i = 0; i < n; i += 2) if (eq(entries[i], key)) return entries[i + 1];
  return NotFound;
}

function trieGet(node: TrieNode, hash: number, key: any, cfg: Resolved): any {
  let shift = 0;
  while (true) {
    if (node.kind === KindChamp) {
      const bit = 1 << ((hash >>> shift) & 31);
      if (node.dataMap & bit) {
        const i = 2 * indexOf(node.dataMap, bit);
        return keyEquals(cfg, node.content[i], key) ? node.content[i + 1] : NotFound;
      }
      if (node.nodeMap & bit) {
        const c = node.content;
        node = c[c.length - 1 - indexOf(node.nodeMap, bit)];
        shift += 5;
        continue;
      }
      return NotFound;
    }
    if (node.hash !== hash) return NotFound;
    return flatGet(node.entries, key, cfg);
  }
}

/** The value of `key`, or `NotFound`. */
export function rootGet(root: Root, key: any, cfg: Resolved): any {
  if (root.kind === KindFlat) return flatGet(root.entries, key, cfg);
  const h = cfg.lookupHash(key);
  if (Number.isNaN(h)) return NotFound;
  return trieGet(root, h, key, cfg);
}

// --- set ---

function mergeTwo(
  k1: any,
  v1: any,
  h1: number,
  k2: any,
  v2: any,
  h2: number,
  shift: number,
  owner: Owner | null,
): TrieNode {
  if (h1 === h2) return new CollisionNode(owner, h1, [k1, v1, k2, v2]);
  // the hashes differ in some bit, so the fragments differ by shift 30 (its two top bits) at the latest
  const f1 = (h1 >>> shift) & 31;
  const f2 = (h2 >>> shift) & 31;
  if (f1 !== f2) {
    return new ChampNode(owner, (1 << f1) | (1 << f2), 0, f1 < f2 ? [k1, v1, k2, v2] : [k2, v2, k1, v1]);
  }
  return new ChampNode(owner, 0, 1 << f1, [mergeTwo(k1, v1, h1, k2, v2, h2, shift + 5, owner)]);
}

function mergeNodeAndEntry(
  sub: CollisionNode,
  subHash: number,
  key: any,
  value: any,
  hash: number,
  shift: number,
  owner: Owner | null,
): ChampNode {
  const f1 = (subHash >>> shift) & 31;
  const f2 = (hash >>> shift) & 31;
  if (f1 !== f2) return new ChampNode(owner, 1 << f2, 1 << f1, [key, value, sub]);
  return new ChampNode(owner, 0, 1 << f1, [mergeNodeAndEntry(sub, subHash, key, value, hash, shift + 5, owner)]);
}

function collisionSet(
  node: CollisionNode,
  shift: number,
  hash: number,
  key: any,
  value: any,
  cfg: Resolved,
  owner: Owner | null,
): TrieNode {
  if (hash !== node.hash) {
    Edit.added = true;
    return mergeNodeAndEntry(node, node.hash, key, value, hash, shift, owner);
  }
  const e = node.entries;
  for (let i = 0; i < e.length; i += 2) {
    if (!keyEquals(cfg, e[i], key)) continue;
    if (cfg.valueEquals(e[i + 1], value)) return node;
    Edit.changed = true;
    if (editable(node, owner)) {
      e[i + 1] = value;
      return node;
    }
    const c = e.slice();
    c[i + 1] = value;
    return new CollisionNode(owner, hash, c);
  }
  Edit.added = true;
  if (editable(node, owner)) {
    e.push(key, value);
    return node;
  }
  const c = e.slice();
  c.push(key, value);
  return new CollisionNode(owner, hash, c);
}

function champSet(
  node: ChampNode,
  shift: number,
  hash: number,
  key: any,
  value: any,
  cfg: Resolved,
  owner: Owner | null,
): ChampNode {
  const bit = 1 << ((hash >>> shift) & 31);
  const content = node.content;
  if (node.dataMap & bit) {
    const i = 2 * indexOf(node.dataMap, bit);
    const k = content[i];
    if (keyEquals(cfg, k, key)) {
      if (cfg.valueEquals(content[i + 1], value)) return node;
      Edit.changed = true;
      if (editable(node, owner)) {
        content[i + 1] = value;
        return node;
      }
      const c = content.slice();
      c[i + 1] = value;
      return new ChampNode(owner, node.dataMap, node.nodeMap, c);
    }
    // push the inline entry down with the new one: the data entry goes, a sub-node comes
    Edit.added = true;
    const sub = mergeTwo(k, content[i + 1], cfg.hash(k), key, value, hash, shift + 5, owner);
    const nodeMap = node.nodeMap | bit;
    const j = indexOf(nodeMap, bit);
    const c = content.slice();
    c.splice(i, 2);
    // sub-nodes from the back: the j-th is at `length - 1 - j`
    c.splice(c.length - j, 0, sub);
    if (editable(node, owner)) {
      node.dataMap ^= bit;
      node.nodeMap = nodeMap;
      node.content = c;
      return node;
    }
    return new ChampNode(owner, node.dataMap ^ bit, nodeMap, c);
  }
  if (node.nodeMap & bit) {
    const pos = content.length - 1 - indexOf(node.nodeMap, bit);
    const sub: TrieNode = content[pos];
    const next =
      sub.kind === KindChamp
        ? champSet(sub, shift + 5, hash, key, value, cfg, owner)
        : collisionSet(sub, shift + 5, hash, key, value, cfg, owner);
    if (next === sub) return node;
    if (editable(node, owner)) {
      content[pos] = next;
      return node;
    }
    const c = content.slice();
    c[pos] = next;
    return new ChampNode(owner, node.dataMap, node.nodeMap, c);
  }
  Edit.added = true;
  const i = 2 * indexOf(node.dataMap, bit);
  if (editable(node, owner)) {
    content.splice(i, 0, key, value);
    node.dataMap |= bit;
    return node;
  }
  const c = new Array(content.length + 2);
  for (let a = 0; a < i; a++) c[a] = content[a];
  c[i] = key;
  c[i + 1] = value;
  for (let a = i; a < content.length; a++) c[a + 2] = content[a];
  return new ChampNode(owner, node.dataMap | bit, node.nodeMap, c);
}

function promote(entries: any[], key: any, value: any, cfg: Resolved, owner: Owner | null): ChampNode {
  // built in place under a token of its own when the edit is persistent: no one else holds it
  const o = owner ?? { dead: false };
  let root = new ChampNode(o, 0, 0, []);
  for (let i = 0; i < entries.length; i += 2) {
    root = champSet(root, 0, cfg.hash(entries[i]), entries[i], entries[i + 1], cfg, o);
  }
  return champSet(root, 0, cfg.hash(key), key, value, cfg, o);
}

/** `root` with `key` set; `Edit` says whether it was added or changed. */
export function rootSet(root: Root, size: number, key: any, value: any, cfg: Resolved, owner: Owner | null): Root {
  resetEdit();
  if (root.kind === KindFlat) {
    const e = root.entries;
    const n = e.length;
    for (let i = 0; i < n; i += 2) {
      if (!keyEquals(cfg, e[i], key)) continue;
      if (cfg.valueEquals(e[i + 1], value)) return root;
      Edit.changed = true;
      if (editable(root, owner)) {
        e[i + 1] = value;
        return root;
      }
      const c = e.slice();
      c[i + 1] = value;
      return new FlatNode(owner, c);
    }
    if (size >= FlatThreshold) {
      const r = promote(e, key, value, cfg, owner);
      Edit.added = true;
      Edit.changed = false;
      return r;
    }
    Edit.added = true;
    if (editable(root, owner)) {
      e.push(key, value);
      return root;
    }
    const c = new Array(n + 2);
    for (let i = 0; i < n; i++) c[i] = e[i];
    c[n] = key;
    c[n + 1] = value;
    return new FlatNode(owner, c);
  }
  return champSet(root, 0, cfg.hash(key), key, value, cfg, owner);
}

// --- delete ---

/** A node holding one entry and nothing else, which its parent inlines (CHAMP's compaction). */
function isSingleton(node: TrieNode) {
  return node.kind === KindChamp ? node.nodeMap === 0 && bitCount(node.dataMap) === 1 : node.entries.length === 2;
}

function singletonContent(node: TrieNode) {
  return node.kind === KindChamp ? node.content : node.entries;
}

function collisionDelete(node: CollisionNode, hash: number, key: any, cfg: Resolved, owner: Owner | null): TrieNode {
  if (hash !== node.hash) return node;
  const e = node.entries;
  for (let i = 0; i < e.length; i += 2) {
    if (!keyEquals(cfg, e[i], key)) continue;
    Edit.removed = true;
    if (editable(node, owner)) {
      e.splice(i, 2);
      return node;
    }
    const c = e.slice();
    c.splice(i, 2);
    return new CollisionNode(owner, hash, c);
  }
  return node;
}

function champDelete(
  node: ChampNode,
  shift: number,
  hash: number,
  key: any,
  cfg: Resolved,
  owner: Owner | null,
): TrieNode {
  const bit = 1 << ((hash >>> shift) & 31);
  const content = node.content;
  if (node.dataMap & bit) {
    const i = 2 * indexOf(node.dataMap, bit);
    if (!keyEquals(cfg, content[i], key)) return node;
    Edit.removed = true;
    // the last entry beside a lone collision node: the collision node goes up (README §2)
    if (node.dataMap === bit && content.length === 3 && content[2].kind === KindCollision) return content[2];
    if (editable(node, owner)) {
      content.splice(i, 2);
      node.dataMap ^= bit;
      return node;
    }
    const c = content.slice();
    c.splice(i, 2);
    return new ChampNode(owner, node.dataMap ^ bit, node.nodeMap, c);
  }
  if (!(node.nodeMap & bit)) return node;
  const j = indexOf(node.nodeMap, bit);
  const pos = content.length - 1 - j;
  const sub: TrieNode = content[pos];
  const next =
    sub.kind === KindChamp
      ? champDelete(sub, shift + 5, hash, key, cfg, owner)
      : collisionDelete(sub, hash, key, cfg, owner);
  if (!Edit.removed) return node;
  const onlyChild = node.dataMap === 0 && content.length === 1;
  if (isSingleton(next)) {
    // a lone sub-node goes up whole; otherwise its entry moves inline at this slot
    if (onlyChild) return next;
    const s = singletonContent(next);
    const i = 2 * indexOf(node.dataMap, bit);
    const c = content.slice();
    c.splice(pos, 1);
    c.splice(i, 0, s[0], s[1]);
    if (editable(node, owner)) {
      node.dataMap |= bit;
      node.nodeMap ^= bit;
      node.content = c;
      return node;
    }
    return new ChampNode(owner, node.dataMap | bit, node.nodeMap ^ bit, c);
  }
  // a collision node alone in a node goes up too: the shape insertion makes (README §2)
  if (onlyChild && next.kind === KindCollision) return next;
  if (next === sub) return node;
  if (editable(node, owner)) {
    content[pos] = next;
    return node;
  }
  const c = content.slice();
  c[pos] = next;
  return new ChampNode(owner, node.dataMap, node.nodeMap, c);
}

function flatten(root: ChampNode, owner: Owner | null): FlatNode {
  const out: any[] = [];
  pushEntries(root, out);
  return new FlatNode(owner, out);
}

function pushEntries(node: TrieNode, out: any[]) {
  if (node.kind === KindCollision) {
    for (let i = 0; i < node.entries.length; i++) out.push(node.entries[i]);
    return;
  }
  const c = node.content;
  const d = 2 * bitCount(node.dataMap);
  for (let i = 0; i < d; i++) out.push(c[i]);
  for (let i = c.length - 1; i >= d; i--) pushEntries(c[i], out);
}

/** `root` without `key`; `Edit.removed` says whether it was there. */
export function rootDelete(root: Root, size: number, key: any, cfg: Resolved, owner: Owner | null): Root {
  resetEdit();
  if (root.kind === KindFlat) {
    const e = root.entries;
    for (let i = 0; i < e.length; i += 2) {
      if (!keyEquals(cfg, e[i], key)) continue;
      Edit.removed = true;
      if (editable(root, owner)) {
        e.splice(i, 2);
        return root;
      }
      const c = e.slice();
      c.splice(i, 2);
      return new FlatNode(owner, c);
    }
    return root;
  }
  const h = cfg.lookupHash(key);
  if (Number.isNaN(h)) return root;
  const next = champDelete(root, 0, h, key, cfg, owner);
  if (!Edit.removed) return root;
  // demoted to flat at the threshold (README §2)
  if (size - 1 <= FlatThreshold) return flatten(next as ChampNode, owner);
  if (next.kind === KindCollision) return new ChampNode(owner, 0, 1 << (h & 31), [next]);
  return next;
}

// --- walks ---

export function forEachEntry(node: Root | TrieNode, fn: (this: any, v: any, k: any) => void, ctx: any) {
  if (node.kind === KindChamp) {
    const c = node.content;
    const d = 2 * bitCount(node.dataMap);
    for (let i = 0; i < d; i += 2) fn.call(ctx, c[i + 1], c[i]);
    for (let i = c.length - 1; i >= d; i--) forEachEntry(c[i], fn, ctx);
    return;
  }
  const e = node.entries;
  for (let i = 0; i < e.length; i += 2) fn.call(ctx, e[i + 1], e[i]);
}

/** Calls `fn(k, v)` for every entry until it returns `false`; `false` when it did. */
export function everyEntry(node: Root | TrieNode, fn: (k: any, v: any) => boolean): boolean {
  if (node.kind === KindChamp) {
    const c = node.content;
    const d = 2 * bitCount(node.dataMap);
    for (let i = 0; i < d; i += 2) if (!fn(c[i], c[i + 1])) return false;
    for (let i = c.length - 1; i >= d; i--) if (!everyEntry(c[i], fn)) return false;
    return true;
  }
  const e = node.entries;
  for (let i = 0; i < e.length; i += 2) if (!fn(e[i], e[i + 1])) return false;
  return true;
}

/** Order-independent: a sum of the entries' hashes (README §4). */
export function hashEntries(root: Root, cfg: Resolved, withValues: boolean): number {
  let acc = 0;
  const vh = withValues ? cfg.valueHash : undefined;
  everyEntry(root, (k, v) => {
    acc = (acc + (vh ? hashEntry(cfg.hash(k), vh(v)) : cfg.hash(k))) | 0;
    return true;
  });
  return acc;
}

function lookupEquals(a: Root | TrieNode, b: Root, cfg: Resolved, values: boolean): boolean {
  return everyEntry(a, (k, v) => {
    const w = rootGet(b, k, cfg);
    return w !== NotFound && (!values || cfg.valueEquals(v, w));
  });
}

function nodeEquals(a: TrieNode, b: TrieNode, cfg: Resolved, values: boolean): boolean {
  if (a === b) return true;
  if (a.kind !== b.kind) return false;
  if (a.kind === KindCollision) {
    const e = (b as CollisionNode).entries;
    if (a.hash !== (b as CollisionNode).hash || a.entries.length !== e.length) return false;
    for (let i = 0; i < a.entries.length; i += 2) {
      const w = flatGet(e, a.entries[i], cfg);
      if (w === NotFound || (values && !cfg.valueEquals(a.entries[i + 1], w))) return false;
    }
    return true;
  }
  const bn = b as ChampNode;
  if (a.dataMap !== bn.dataMap || a.nodeMap !== bn.nodeMap) return false;
  const c = a.content;
  const d = bn.content;
  const n = 2 * bitCount(a.dataMap);
  for (let i = 0; i < n; i += 2) {
    if (!keyEquals(cfg, c[i], d[i])) return false;
    if (values && !cfg.valueEquals(c[i + 1], d[i + 1])) return false;
  }
  for (let i = n; i < c.length; i++) if (!nodeEquals(c[i], d[i], cfg, values)) return false;
  return true;
}

/** Equal contents, the sizes equal and the configs the same (README §4): a walk of the tries, lookups otherwise. */
export function rootEquals(a: Root, b: Root, cfg: Resolved, values: boolean): boolean {
  if (a === b) return true;
  if (a.kind === KindChamp && b.kind === KindChamp) return nodeEquals(a, b, cfg, values);
  return lookupEquals(a, b, cfg, values);
}
