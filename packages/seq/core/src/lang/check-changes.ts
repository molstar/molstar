/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The `seq.lang.*` checks of a change (_spec/README.md, Language, Checks of a change): `checkSeqLang` run on what a
// patch's `Changes` reach, as `repr/patch`'s validator, its cost in proportion to the changes. What a change reaches:
// the nodes it made, set or moved, the nodes that name them or named a node of a subtree it removed or moved (through a
// reverse index of the spec's references, kept per revision and derived from the last one along its changes), on
// through the user refs that name those and through refs with segments after them (`/refs/p/x`), a template's
// layers under a track whose params changed, and the sections and refs of a chain a changed one closes. Only what the
// change introduced is reported: the same checks on the same nodes before it are taken away, so an error already
// present never refuses an edit.

import { HashMap, HashSet, type MutableHashMap } from '@molstar/core/data/immutable';
import type { Changes, PatchValidator } from '@molstar/core/repr/patch';
import { type Address, type Issue, type Node, printAddress, type TreePath, type Value } from '@molstar/core/repr/tree';
import { checkSubtreeDepth, DefaultSpecBounds, type SpecBounds } from './bounds.js';
import {
  checkFollow,
  checkPlaced,
  fieldsOf,
  idRunOf,
  type Placed,
  type Report,
  reporter,
  type SeqLookup,
} from './check.js';
import { DefaultRefDepth, followRef, RefsId } from './refs.js';

/**
 * The spec's references by what they name: `exact`, the referring nodes' addresses by the id run they name (a
 * persistent set each, so a derive costs what changed), and `out`, each referring node's own targets, to take them back
 * out when it changes.
 */
interface RefIndex {
  readonly exact: HashMap<string, HashSet<string>>;
  readonly out: HashMap<string, readonly string[]>;
}

const Indexes = new WeakMap<Node, RefIndex>();

type Walkable = { readonly type: string; readonly of?: unknown; readonly fields?: unknown };

/** The address strings a value holds where its param type is a ref: in lists, unions, records and maps too. */
function collectRefs(type: Walkable | undefined, value: Value | undefined, out: string[]): void {
  if (!type || value === undefined || value === null) return;
  switch (type.type) {
    case 'ref':
      if (typeof value === 'string') out.push(value);
      return;
    case 'list':
      if (Array.isArray(value)) for (const v of value) collectRefs(type.of as Walkable, v, out);
      return;
    case 'union':
      for (const m of type.of as Walkable[]) {
        if (m.type === 'ref' && typeof value === 'string') return void out.push(value);
        if (m.type === 'list' && Array.isArray(value)) return collectRefs(m, value, out);
      }
      return;
    case 'record':
      if (typeof value === 'object' && !Array.isArray(value))
        for (const [k, f] of Object.entries(type.fields as Record<string, Walkable>))
          collectRefs(f, (value as Record<string, Value>)[k], out);
      return;
    case 'map':
      if (typeof value === 'object' && !Array.isArray(value))
        for (const v of Object.values(value as Record<string, Value>)) collectRefs(type.of as Walkable, v, out);
      return;
  }
}

/**
 * The id runs a node's references name, printed (`/data/p53`), once each; one through a ref with segments after it
 * (`/refs/p/x`) names that ref too, so a change of the ref reaches it.
 */
function targetsOf(node: Node): string[] {
  const fields = fieldsOf(node.kind) as Readonly<Record<string, Walkable>> | undefined;
  if (!fields || !node.params) return [];
  const raw: string[] = [];
  for (const key of Object.keys(node.params)) collectRefs(fields[key], node.params[key], raw);
  const out: string[] = [];
  const push = (ids: readonly string[]) => {
    const t = printAddress(ids);
    if (!out.includes(t)) out.push(t);
  };
  for (const a of raw) {
    const ids = idRunOf(a);
    if (!ids.length) continue;
    push(ids);
    if (ids.length > 2 && ids[0] === RefsId) push(ids.slice(0, 2));
  }
  return out;
}

const NoReferrers: HashSet<string> = HashSet.empty();

/** Edits an index in place (transients), each node's targets added or taken out. */
class IndexEdit {
  private readonly exact: MutableHashMap<string, HashSet<string>>;
  private readonly out: MutableHashMap<string, readonly string[]>;
  constructor(index: RefIndex) {
    this.exact = index.exact.asMutable();
    this.out = index.out.asMutable();
  }
  add(address: string, node: Node): void {
    const targets = targetsOf(node);
    if (!targets.length) return;
    for (const t of targets) this.exact.update(t, (s) => (s ?? NoReferrers).add(address));
    this.out.set(address, targets);
  }
  remove(address: string): void {
    const targets = this.out.get(address);
    if (!targets) return;
    for (const t of targets) {
      const s = this.exact.get(t)?.delete(address);
      if (!s?.size) this.exact.delete(t);
      else this.exact.set(t, s);
    }
    this.out.delete(address);
  }
  done(): RefIndex {
    return { exact: this.exact.asImmutable(), out: this.out.asImmutable() };
  }
}

/** Every addressed node of a subtree with its address (the subtree's root at `address`). */
function walk(node: Node, address: readonly string[], f: (address: string, node: Node) => void): void {
  f(printAddress(address), node);
  for (const c of node.children ?? []) if (c.id !== undefined) walk(c, [...address, c.id], f);
}

function indexOf(root: Node): RefIndex {
  let index = Indexes.get(root);
  if (!index) {
    // built once a spec, as plain maps of sets, then made persistent: what later revisions derive from
    const exact = new Map<string, string[]>();
    const out = HashMap.empty<string, readonly string[]>().asMutable();
    walk(root, [], (a, n) => {
      const targets = targetsOf(n);
      if (!targets.length) return;
      for (const t of targets) {
        const s = exact.get(t);
        if (s) s.push(a);
        else exact.set(t, [a]);
      }
      out.set(a, targets);
    });
    const ex = HashMap.empty<string, HashSet<string>>().asMutable();
    for (const [t, s] of exact) ex.set(t, HashSet.of(s));
    index = { exact: ex.asImmutable(), out: out.asImmutable() };
    Indexes.set(root, index);
  }
  return index;
}

const ChildIndices = new WeakMap<readonly Node[], Map<string, number>>();

/** A child's index by its id: a scan of a few children, else a map per children array, kept. */
function childIndex(children: readonly Node[] | undefined, id: string): number {
  if (!children) return -1;
  if (children.length <= 8) return children.findIndex((c) => c.id === id);
  let m = ChildIndices.get(children);
  if (!m) {
    m = new Map();
    children.forEach((c, i) => {
      if (c.id !== undefined && !m!.has(c.id)) m!.set(c.id, i);
    });
    ChildIndices.set(children, m);
  }
  return m.get(id) ?? -1;
}

/** A place in a tree by address (ids only): the node, its step path, its parent. */
function placeIn(root: Node, ids: readonly string[]): Placed | undefined {
  let node = root;
  let parent: Node | undefined;
  const path: (string | number)[] = [];
  for (const id of ids) {
    const i = childIndex(node.children, id);
    if (i < 0) return undefined;
    parent = node;
    node = node.children![i];
    path.push('children', i);
  }
  return { node, path: path as TreePath, address: printAddress(ids), ...(parent ? { parent } : {}) };
}

/** A lookup over a tree by exact address, each place found once. */
function lookupOf(spec: Node): SeqLookup {
  const found = new Map<string, Placed | undefined>();
  return {
    spec,
    get(address) {
      if (found.has(address)) return found.get(address);
      const p = address === '/' ? { node: spec, path: [] as TreePath, address: '/' } : placeIn(spec, split(address));
      found.set(address, p);
      return p;
    },
  };
}

const split = (address: string) => address.split('/').filter(Boolean);
/** A user ref's address (`/refs/<id>`). */
const isRef = (address: string) => {
  const ids = split(address);
  return ids.length === 2 && ids[0] === RefsId;
};
const under = (a: Address, prefix: Address) => a.length >= prefix.length && prefix.every((s, i) => a[i] === s);

/** An address through the moves, from the old root's to the new's (`forward`) or back. */
function throughMoves(address: string, changes: Changes, forward: boolean): string {
  let ids = split(address);
  const moves = forward ? changes.moved : [...changes.moved].reverse();
  for (const m of moves) {
    const [from, to] = forward ? [m.from, m.to] : [m.to, m.from];
    if (under(ids, from)) ids = [...to, ...ids.slice(from.length)];
  }
  return printAddress(ids);
}

/** The reverse index of `after`, derived from `before`'s along the changes when it is kept. */
function deriveIndex(before: Node, after: Node, changes: Changes): RefIndex {
  const kept = Indexes.get(after);
  if (kept) return kept;
  const prior = Indexes.get(before);
  if (!prior) return indexOf(after);
  const edit = new IndexEdit(prior);
  const old = (address: Address) => {
    const at = placeIn(before, address);
    if (at) walk(at.node, address, (a) => edit.remove(a));
  };
  const now = (address: Address) => {
    const at = placeIn(after, address);
    if (at) walk(at.node, address, (a, n) => edit.add(a, n));
  };
  for (const r of changes.removed) old(r);
  for (const m of changes.moved) old(m.from);
  for (const m of changes.made) old(split(throughMoves(printAddress(m), changes, false)));
  for (const c of changes.params) edit.remove(throughMoves(printAddress(c.at), changes, false));
  for (const m of changes.moved) now(m.to);
  for (const m of changes.made) now(m);
  for (const c of changes.params) {
    const at = placeIn(after, c.at);
    if (at) {
      edit.remove(printAddress(c.at));
      edit.add(printAddress(c.at), at.node);
    }
  }
  const index = edit.done();
  Indexes.set(after, index);
  return index;
}

/** An issue's key: its code, its node's address and where in the node (the path after the node's). */
function keyed(issues: Issue[], keys: string[], rename: (address: string) => string = (a) => a): Report {
  const push = reporter(issues);
  return (code, message, at, ...tail) => {
    push(code, message, at, ...tail);
    keys.push(`${code}|${rename(at.address ?? '/')}|${tail.join('/')}`);
  };
}

/**
 * The `seq.lang.*` issues a change introduced (_spec/README.md, Language, Checks of a change): `checkSeqLang`'s, on the
 * nodes `changes` reaches in `after`, less those the same nodes had in `before`. For a valid, normalized `after`.
 */
export function checkSeqChanges(
  before: Node,
  after: Node,
  changes: Changes,
  bounds: SpecBounds = DefaultSpecBounds,
): Issue[] {
  // the nesting bounds, as the whole spec's check holds them, before anything walks the changes (a subtree too deep to
  // walk is refused alone): the subtrees a change made or moved, the params it set
  const look0 = lookupOf(after);
  const deep: Issue[] = [];
  const nest = (address: Address, children: boolean) => {
    const at = look0.get(printAddress(address));
    if (at) deep.push(...checkSubtreeDepth(at.node, at.path, address.length, bounds, children));
  };
  for (const m of changes.made) nest(m, true);
  for (const m of changes.moved) nest(m.to, true);
  for (const c of changes.params) nest(c.at, false);
  if (deep.length) return deep;
  const index = deriveIndex(before, after, changes);
  const look = look0;
  const affected = new Set<string>();
  const queue: string[] = [];
  const touch = (address: string) => {
    if (affected.has(address) || !look.get(address)) return;
    affected.add(address);
    queue.push(address);
  };
  /**
   * What names `target`: the nodes whose refs name it, and through a ref with segments after it (`/refs/p/x`, where
   * `/refs/p` names a prefix of `target`) what names that, at most the ref depth deep.
   */
  const referrers = (target: string, depth = 0) => {
    index.exact.get(target)?.forEach((r) => touch(r));
    if (depth >= DefaultRefDepth) return;
    const ids = split(target);
    for (let n = 1; n < ids.length; n++) {
      const through = index.exact.get(printAddress(ids.slice(0, n)));
      if (!through) continue;
      const rest = ids.slice(n);
      through.forEach((r) => {
        if (isRef(r)) referrers(printAddress([...split(r), ...rest]), depth + 1);
      });
    }
  };
  /** What named any node of a subtree of `before` (gone, or moved away). */
  const referrersOfOld = (address: Address) => {
    const at = placeIn(before, address);
    if (at) walk(at.node, address, (a) => referrers(a));
  };

  for (const m of changes.made) {
    const at = look.get(printAddress(m));
    if (at) walk(at.node, m, (a) => touch(a));
  }
  for (const m of changes.moved) {
    const at = look.get(printAddress(m.to));
    if (at) walk(at.node, m.to, (a) => touch(a));
  }
  for (const c of changes.params) touch(printAddress(c.at));
  // what names a node that is gone (a ref left dangling under `allowDangling`, or one a move took)
  for (const r of changes.removed) referrersOfOld(r);
  for (const m of changes.moved) referrersOfOld(m.from);

  // what the changed nodes reach: what names them, on through refs; a template's layers; a chain a change closes
  const direct = queue.length;
  for (let i = 0; i < queue.length; i++) {
    const address = queue[i];
    const at = look.get(address)!;
    if (i < direct || at.node.kind === 'ref') referrers(address);
    if (at.node.kind === 'track') for (const c of at.node.children ?? []) touch(`${address}/${c.id}`);
    if (at.node.kind === 'section') {
      const chain: string[] = [address];
      for (let a = sectionFollow(look, address); a !== undefined; ) {
        const back = chain.indexOf(a);
        if (back >= 0) {
          for (const s of chain.slice(back)) touch(s);
          break;
        }
        chain.push(a);
        a = sectionFollow(look, a);
      }
    }
    if (at.node.kind === 'ref') {
      const f = followRef(after, address, DefaultRefDepth);
      if ('refused' in f) for (const via of f.via) touch(via);
    }
  }

  const issues: Issue[] = [];
  const keys: string[] = [];
  const report = keyed(issues, keys);
  for (const a of affected) {
    const at = look.get(a)!;
    checkPlaced(at, look, report);
    checkFollow(at, look, report);
  }
  if (!issues.length) return issues;

  // what the same nodes had before: not introduced by this change
  const was: Issue[] = [];
  const wasKeys: string[] = [];
  const old = lookupOf(before);
  const reportOld = keyed(was, wasKeys, (a) => throughMoves(a, changes, true));
  for (const a of affected) {
    const at = old.get(throughMoves(a, changes, false));
    if (!at) continue;
    checkPlaced(at, old, reportOld);
    checkFollow(at, old, reportOld);
  }
  const known = new Set(wasKeys);
  return issues.filter((_, i) => !known.has(keys[i]));
}

function sectionFollow(look: SeqLookup, address: string): string | undefined {
  const at = look.get(address);
  if (at?.node.kind !== 'section') return undefined;
  const f = at.node.params?.follow;
  if (typeof f !== 'string') return undefined;
  const ids = idRunOf(f);
  return ids.length ? printAddress(ids) : undefined;
}

/**
 * The language's checks as `repr/patch`'s validator (`ApplyOptions.validators`): `checkSeqChanges`, so a patch is
 * refused for a `seq.lang.*` error it introduced and never for one already there.
 */
export const SeqPatchValidator: PatchValidator = (before, after, changes) => checkSeqChanges(before, after, changes);

/** `SeqPatchValidator` with other bounds than the defaults (a composition's `bounds` configuration). */
export const seqPatchValidator =
  (bounds: SpecBounds): PatchValidator =>
  (before, after, changes) =>
    checkSeqChanges(before, after, changes, bounds);
