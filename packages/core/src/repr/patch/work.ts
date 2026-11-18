/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The ops at work (_spec/README.md "Patches", Applying): a draft over the root, made of the nodes along the paths the
// ops touch, each draft knowing where it came from (the old root's node and its address, or made by the patch). Every
// op is checked whole before it changes anything, so a refused op leaves the draft as it was. `finalize` gives the new
// root by structural sharing: an untouched subtree is the old object. The draft also gives the net `Changes` and
// records each op's inverse against the state just before it.

import type { Changes, ParamChange } from '../schema/index.js';
import {
  type Address,
  hashNode,
  IdHint,
  IdPattern,
  type Issue,
  isColumn,
  isNode,
  isRecord,
  isReservedKey,
  issueOf,
  KindPattern,
  MaxDepth,
  type Node,
  own,
  type Params,
  type TreePath,
  type Value,
  valuesEqual,
} from '../tree/index.js';
import { childIndex } from '../tree/address.js';
import { type Op, type PatchCode, PatchIssues, type Place, type Where } from './patch.js';

/** A node of the draft. Its unchanged parts are read from `base`; `kids` is materialized on first use. */
export interface Draft {
  /** The node this draft reads its unchanged parts from. */
  base: Node;
  /** The old root's node this one descends from (kept through moves and a replace); none when made from nothing. */
  origin: Node | undefined;
  /** The origin's address in the old root. */
  oldAddress: Address | undefined;
  /** Inserted or replaced, or inside a node that was. */
  made: boolean;
  /** Moved to another parent or renamed by an op. */
  relocated: boolean;
  kind: string;
  id: string | undefined;
  params: Params | undefined;
  custom: Params | undefined;
  /** Keys changed by ops or ref rewrites, by where. */
  touched: { params?: Set<string>; custom?: Set<string> } | undefined;
  /** Params keys an op set (not a ref rewrite): what the patch wrote. */
  written: Set<string> | undefined;
  kids: Draft[] | undefined;
}

function draftOf(node: Node, origin: Node | undefined, oldAddress: Address | undefined, made: boolean): Draft {
  return {
    base: node,
    origin,
    oldAddress,
    made,
    relocated: false,
    kind: node.kind,
    id: node.id,
    params: node.params,
    custom: node.custom,
    touched: undefined,
    written: undefined,
    kids: undefined,
  };
}

/** Normalization of what a patch makes, when a language is given (Applying, step 3). */
export interface Normalizer {
  node(node: Node): Node;
  value(kind: string, key: string, value: Value): Value;
}

export const isUnder = (a: Address, prefix: Address) => a.length >= prefix.length && prefix.every((s, i) => a[i] === s);

export const sameAddress = (a: Address, b: Address) => a.length === b.length && isUnder(a, b);

/** How many levels a value nests (lists, records and nodes count one each), or more than `budget` past it. */
function valueHeight(v: Value, budget: number): number {
  if (typeof v !== 'object' || v === null || isColumn(v)) return 0;
  if (isNode(v)) return nodeHeight(v, budget);
  if (budget <= 0) return 1;
  let h = 0;
  const items = Array.isArray(v) ? v : Object.values(v as Params);
  for (const x of items) {
    h = Math.max(h, valueHeight(x, budget - 1));
    if (h > budget - 1) return budget + 1;
  }
  return h + 1;
}

/**
 * A node's height, 1 for a leaf, its values' nesting counted; exact up to `budget`, and past it some number above it
 * (so a hostile depth costs no more than the bound).
 */
export function nodeHeight(n: Node, budget: number): number {
  if (budget <= 0) return 1;
  let h = 0;
  const parts: Value[] = [];
  for (const r of [n.params, n.custom]) if (r) parts.push(...Object.values(r));
  if (n.args) parts.push(...n.args);
  for (const x of parts) {
    h = Math.max(h, valueHeight(x, budget - 1));
    if (h > budget - 1) return budget + 1;
  }
  for (const c of n.children ?? []) {
    h = Math.max(h, nodeHeight(c, budget - 1));
    if (h > budget - 1) return budget + 1;
  }
  return h + 1;
}

const isAddress = (a: unknown): a is Address => Array.isArray(a) && a.every((s) => typeof s === 'string');
const isPlace = (p: unknown): p is Place =>
  isRecord(p as Value) &&
  ((Object.hasOwn(p as object, 'after') &&
    ((p as { after: unknown }).after === null || typeof (p as { after: unknown }).after === 'string')) ||
    (p as { last?: unknown }).last === true);
const isWhere = (w: unknown) => w === undefined || w === 'params' || w === 'custom';

class Refused {
  constructor(readonly issue: Issue) {}
}

/** The draft of a root and the ops applied to it so far. */
export class Work {
  readonly root: { draft: Draft };
  /** Each applied op's inverse, in the order applied (the patch's inverse is them reversed). */
  readonly inverse: Op[] = [];
  /** Address changes, in order: what refs follow. */
  readonly moves: { readonly from: Address; readonly to: Address }[] = [];
  /** Subtrees taken out (removed, or a replaced node's kids): where `removed` is read from. */
  private readonly gone: { draft: Draft; cover: Address | undefined }[] = [];
  /** Whether an op inserted, removed, replaced, renamed or moved a node to another parent. */
  structural = false;
  private index = 0;

  constructor(
    readonly old: Node,
    private readonly normalizer?: Normalizer,
  ) {
    this.root = { draft: draftOf(old, old, [], false) };
  }

  private refuse(code: PatchCode, message: string, ...at: (string | number)[]): never {
    throw new Refused(
      issueOf(PatchIssues, code, message, {
        path: ['children', this.index, ...(at.length ? ['params', ...at] : [])],
        hint: PatchIssues.codes[code].hint,
      }),
    );
  }

  kids(d: Draft): Draft[] {
    if (d.kids) return d.kids;
    const children = d.base.children ?? [];
    d.kids = children.map((c) => {
      // the old root's children are the origin's own: a made node's (a replacement's too) are new
      const inner = d.base === d.origin ? c : undefined;
      const oldAddress = inner && d.oldAddress && c.id !== undefined ? [...d.oldAddress, c.id] : undefined;
      return draftOf(c, inner, oldAddress, d.made);
    });
    return d.kids;
  }

  private kidAt(d: Draft, id: string): number {
    if (!d.kids) return childIndex(d.base.children, id);
    for (let i = 0; i < d.kids.length; i++) if (d.kids[i].id === id) return i;
    return -1;
  }

  /** The draft at an address (materialized along the path), or undefined. */
  find(address: Address): Draft | undefined {
    let d = this.root.draft;
    for (const id of address) {
      const i = this.kidAt(d, id);
      if (i < 0) return undefined;
      d = this.kids(d)[i];
    }
    return d;
  }

  /** The draft at a step path (`['children', 2, 'children', 0]`), materialized; the rest of the path is returned. */
  findPath(path: TreePath): { draft: Draft; rest: TreePath } | undefined {
    let d = this.root.draft;
    let i = 0;
    while (i + 1 < path.length && path[i] === 'children' && typeof path[i + 1] === 'number') {
      const kids = this.kids(d);
      d = kids[path[i + 1] as number];
      if (!d) return undefined;
      i += 2;
    }
    return { draft: d, rest: path.slice(i) };
  }

  /** Applies op number `index`, or returns why it is refused (nothing changed then). */
  apply(op: Op, index: number): Issue | undefined {
    this.index = index;
    try {
      this.op(op);
      return undefined;
    } catch (e) {
      if (e instanceof Refused) return e.issue;
      throw e;
    }
  }

  private op(op: Op) {
    if (!op || typeof op !== 'object') this.refuse('repr.patch.op', 'An op is an object with `op`');
    switch (op.op) {
      case 'set':
        return this.set(op.at, op.key, op.value, op.where);
      case 'unset':
        return this.unset(op.at, op.key, op.where);
      case 'insert':
        return this.insert(op.parent, op.place, op.node);
      case 'remove':
        return this.remove(op.at);
      case 'move':
        return this.move(op.at, op.parent, op.place, op.rename);
      case 'replace':
        return this.replace(op.at, op.node);
      default:
        this.refuse('repr.patch.op', 'Not an op: set, unset, insert, remove, move or replace', 'op');
    }
  }

  private target(at: unknown, key = 'at'): Draft {
    if (!isAddress(at)) this.refuse('repr.patch.op', 'An address is a list of ids', key);
    const d = this.find(at);
    if (!d) this.refuse('repr.patch.target', 'No node at this address', key);
    return d;
  }

  private keyOf(key: unknown, where: unknown) {
    if (typeof key !== 'string') this.refuse('repr.patch.op', 'A key is a string', 'key');
    if (!isWhere(where)) this.refuse('repr.patch.op', '`where` is params or custom', 'where');
    if (isReservedKey(key)) this.refuse('repr.patch.key', 'A reserved key', 'key');
  }

  private touch(d: Draft, where: Where, key: string, written: boolean) {
    d.touched ??= {};
    (d.touched[where] ??= new Set()).add(key);
    if (written && where === 'params') (d.written ??= new Set()).add(key);
  }

  private set(at: Address, key: string, value: Value, where: Where = 'params') {
    this.keyOf(key, where);
    const d = this.target(at);
    if (value === undefined) this.refuse('repr.patch.op', 'A set has a value', 'value');
    if (at.length + 2 + valueHeight(value, MaxDepth) > MaxDepth)
      this.refuse('repr.patch.depth', `The value nests past ${MaxDepth} levels`, 'value');
    const v = where === 'params' && this.normalizer ? this.normalizer.value(d.kind, key, value) : value;
    const record = d[where];
    const old = record ? own(record, key) : undefined;
    if (old !== undefined && valuesEqual(old, v)) return; // a no-op keeps the node
    this.inverse.push(old === undefined ? { op: 'unset', at, key, where } : { op: 'set', at, key, value: old, where });
    d[where] = { ...record, [key]: v };
    this.touch(d, where, key, true);
  }

  private unset(at: Address, key: string, where: Where = 'params') {
    this.keyOf(key, where);
    const d = this.target(at);
    const record = d[where];
    const old = record ? own(record, key) : undefined;
    if (old === undefined) return;
    this.inverse.push({ op: 'set', at, key, value: old, where });
    const next = { ...record };
    delete next[key];
    d[where] = Object.keys(next).length ? next : undefined;
    this.touch(d, where, key, false);
  }

  /** Where `place` puts a node among `kids` (the moved node already out of them); refuses a place naming no sibling. */
  private position(kids: readonly Draft[], place: Place, skip?: Draft): number {
    if ('last' in place && place.last === true) return kids.filter((k) => k !== skip).length;
    const after = (place as { after: string | null }).after;
    if (after === null) return 0;
    let n = 0;
    for (const k of kids) {
      if (k === skip) continue;
      n++;
      if (k.id === after) return n;
    }
    this.refuse('repr.patch.place', 'The place names no sibling', 'place');
  }

  private addressable(d: Draft, key: string) {
    for (const k of this.kids(d))
      if (k.id === undefined) this.refuse('repr.patch.unaddressable', 'A child of this parent has no id', key);
  }

  private checkNode(node: unknown, key: string): asserts node is Node {
    if (!isNode(node) || !KindPattern.test(node.kind)) this.refuse('repr.patch.op', 'Not a node', key);
  }

  private insert(parent: Address, place: Place, node: Node) {
    const p = this.target(parent, 'parent');
    if (!isPlace(place)) this.refuse('repr.patch.op', 'A place is {after} or {last: true}', 'place');
    this.checkNode(node, 'node');
    this.addressable(p, 'parent');
    if (node.id === undefined || !IdPattern.test(node.id))
      this.refuse('repr.patch.id', `An inserted node has an id: ${IdHint}`, 'node');
    const kids = this.kids(p);
    if (kids.some((k) => k.id === node.id)) this.refuse('repr.patch.id', 'A sibling has this id', 'node');
    const i = this.position(kids, place);
    if (parent.length + 1 + nodeHeight(node, MaxDepth) > MaxDepth)
      this.refuse('repr.patch.depth', `The node would sit past ${MaxDepth} levels`, 'node');
    const n = this.normalizer ? this.normalizer.node(node) : node;
    kids.splice(i, 0, draftOf(n, undefined, undefined, true));
    this.inverse.push({ op: 'remove', at: [...parent, node.id] });
    this.structural = true;
  }

  private remove(at: Address) {
    if (!isAddress(at)) this.refuse('repr.patch.op', 'An address is a list of ids', 'at');
    if (!at.length) this.refuse('repr.patch.op', 'The root is never removed', 'at');
    const parent = at.slice(0, -1);
    const p = this.target(parent);
    const kids = this.kids(p);
    const i = this.kidAt(p, at[at.length - 1]);
    if (i < 0) this.refuse('repr.patch.target', 'No node at this address', 'at');
    this.addressable(p, 'at');
    const d = kids[i];
    this.inverse.push({
      op: 'insert',
      parent,
      place: { after: i > 0 ? kids[i - 1].id! : null },
      node: this.finalize(d),
    });
    kids.splice(i, 1);
    this.gone.push({ draft: d, cover: undefined });
    this.structural = true;
  }

  private move(at: Address, parent: Address, place: Place, rename: string | undefined) {
    if (!isAddress(at)) this.refuse('repr.patch.op', 'An address is a list of ids', 'at');
    if (!at.length) this.refuse('repr.patch.op', 'The root is never moved', 'at');
    if (!isPlace(place)) this.refuse('repr.patch.op', 'A place is {after} or {last: true}', 'place');
    if (rename !== undefined && (typeof rename !== 'string' || !IdPattern.test(rename)))
      this.refuse('repr.patch.id', `A new id is an id: ${IdHint}`, 'rename');
    const from = at.slice(0, -1);
    const src = this.target(from);
    const srcKids = this.kids(src);
    const i = this.kidAt(src, at[at.length - 1]);
    if (i < 0) this.refuse('repr.patch.target', 'No node at this address', 'at');
    const dst = this.target(parent, 'parent');
    if (isUnder(parent, at)) this.refuse('repr.patch.cycle', 'A node never moves under itself', 'parent');
    this.addressable(src, 'at');
    this.addressable(dst, 'parent');
    const d = srcKids[i];
    const id = rename ?? d.id!;
    const dstKids = this.kids(dst);
    if (dstKids.some((k) => k !== d && k.id === id)) this.refuse('repr.patch.id', 'A sibling has this id', 'rename');
    const j = this.position(dstKids, place, d);
    const sameParent = src === dst;
    if (
      !sameParent &&
      parent.length > from.length &&
      parent.length + 1 + nodeHeight(this.finalize(d), MaxDepth) > MaxDepth
    )
      this.refuse('repr.patch.depth', `The node would sit past ${MaxDepth} levels`, 'parent');
    if (sameParent && id === d.id && i === j) return; // where it is: a no-op
    this.inverse.push({
      op: 'move',
      at: [...parent, id],
      parent: from,
      place: { after: i > 0 ? srcKids[i - 1].id! : null },
      ...(id !== d.id ? { rename: d.id! } : {}),
    });
    srcKids.splice(i, 1);
    dstKids.splice(j, 0, d);
    if (!sameParent || id !== d.id) {
      d.relocated = true;
      this.structural = true;
      this.moves.push({ from: at, to: [...parent, id] });
    }
    d.id = id;
  }

  private replace(at: Address, node: Node) {
    const d = this.target(at);
    this.checkNode(node, 'node');
    if (node.id !== d.id) this.refuse('repr.patch.id', 'A replacement keeps the id of what it replaces', 'node');
    if (at.length + nodeHeight(node, MaxDepth) > MaxDepth)
      this.refuse('repr.patch.depth', `The node would sit past ${MaxDepth} levels`, 'node');
    const n = this.normalizer ? this.normalizer.node(node) : node;
    const old = this.finalize(d);
    if (hashNode(n) === hashNode(old)) return; // equal by hash: keeps the old object
    this.inverse.push({ op: 'replace', at, node: old });
    const next = draftOf(n, d.origin, d.oldAddress, true);
    next.relocated = d.relocated;
    // what the old node held is made over at its address; a node moved into it from elsewhere is gone
    if (d.kids) for (const k of d.kids) this.gone.push({ draft: k, cover: d.oldAddress });
    if (!at.length) this.root.draft = next;
    else {
      const p = this.find(at.slice(0, -1))!;
      const kids = this.kids(p);
      kids[kids.indexOf(d)] = next;
    }
    this.structural = true;
  }

  /** Sets a value at a path inside a draft's params (a ref rewritten): not written by an op. */
  rewrite(path: TreePath, value: string) {
    const found = this.findPath(path);
    if (!found || found.rest[0] !== 'params' || typeof found.rest[1] !== 'string') return;
    const { draft: d, rest } = found;
    const key = rest[1] as string;
    const old = d.params ? own(d.params, key) : undefined;
    if (old === undefined) return;
    const next = setIn(old, rest.slice(2), value);
    if (next === old) return;
    d.params = { ...d.params, [key]: next };
    this.touch(d, 'params', key, false);
  }

  /** The node a draft stands for now. */
  finalize(d: Draft = this.root.draft): Node {
    const base = d.base;
    let children = base.children;
    if (d.kids) {
      const out = d.kids.map((k) => this.finalize(k));
      const same =
        !!base.children && out.length === base.children.length && out.every((c, i) => c === base.children![i]);
      children = same ? base.children : out;
    }
    if (
      children === base.children &&
      d.kind === base.kind &&
      d.id === base.id &&
      d.params === base.params &&
      d.custom === base.custom
    )
      return base;
    const n: { -readonly [K in keyof Node]: Node[K] } = { kind: d.kind };
    if (d.id !== undefined) n.id = d.id;
    if (d.params !== undefined) n.params = d.params;
    if (d.custom !== undefined) n.custom = d.custom;
    if (children !== undefined && (children.length || base.children)) n.children = children;
    if (base.args !== undefined) n.args = base.args;
    return n;
  }

  /** The net difference from the old root to the draft (Applying, the result). */
  changes(): Changes {
    const made: Address[] = [];
    const params: ParamChange[] = [];
    const children: Address[] = [];
    const moved: { from: Address; to: Address }[] = [];
    const walk = (d: Draft, address: Address | undefined, parentMade: boolean) => {
      if (d.made && !parentMade && address) made.push(address);
      if (d.relocated && d.oldAddress && address && !sameAddress(d.oldAddress, address))
        moved.push({ from: d.oldAddress, to: address });
      if (!d.made && d.origin && address) {
        for (const where of ['params', 'custom'] as const) {
          const keys = d.touched?.[where];
          if (!keys) continue;
          const before = d.origin[where];
          const now = d[where];
          const changed = [...keys].filter(
            (k) => !valuesEqual(before ? own(before, k) : undefined, now ? own(now, k) : undefined),
          );
          if (changed.length) params.push({ at: address, keys: changed.sort(), where });
        }
        if (d.kids && kidsChanged(d)) children.push(address);
      }
      if (d.kids)
        for (const k of d.kids) walk(k, address && k.id !== undefined ? [...address, k.id] : undefined, d.made);
    };
    walk(this.root.draft, [], false);
    const removed: Address[] = [];
    const collect = (d: Draft, cover: Address | undefined) => {
      let c = cover;
      if (d.oldAddress && !(cover && isUnder(d.oldAddress, cover))) {
        removed.push(d.oldAddress);
        c = d.oldAddress;
      }
      if (d.kids) for (const k of d.kids) collect(k, c);
    };
    for (const g of this.gone) collect(g.draft, g.cover);
    return { made, params, children, removed, moved };
  }
}

function kidsChanged(d: Draft): boolean {
  const before = d.origin!.children ?? [];
  const now = d.kids!;
  if (before.length !== now.length) return true;
  for (let i = 0; i < now.length; i++) if (now[i].origin !== before[i] || now[i].id !== before[i].id) return true;
  return false;
}

/** `value` with the leaf at `path` (keys and list indices) set to `leaf`; the same object when nothing changes. */
export function setIn(value: Value, path: TreePath, leaf: Value): Value {
  if (!path.length) return Object.is(value, leaf) ? value : leaf;
  const [head, ...rest] = path;
  if (Array.isArray(value) && typeof head === 'number') {
    const next = setIn(value[head], rest, leaf);
    if (next === value[head]) return value;
    const out = [...value];
    out[head] = next;
    return out;
  }
  if (typeof value === 'object' && value !== null && !isColumn(value) && typeof head === 'string') {
    // a record, or a node in a slot (its params, args and children are read the same way)
    const rec = value as unknown as Record<string, Value>;
    const inner = Object.hasOwn(rec, head) ? rec[head] : undefined;
    if (inner === undefined) return value;
    const next = setIn(inner, rest, leaf);
    return next === inner ? value : ({ ...rec, [head]: next } as Value);
  }
  return value;
}
