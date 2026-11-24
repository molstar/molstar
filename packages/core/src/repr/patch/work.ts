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
  resolveAddress,
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
  /** Params keys an op set (not a ref rewrite, not an unset): what the patch wrote. */
  written: Set<string> | undefined;
  /** Params keys an op set or unset. */
  byOp: Set<string> | undefined;
  kids: Draft[] | undefined;
  /** The kids by id (the first of two with one id), with `kids`. */
  byId: Map<string, Draft> | undefined;
  /** How many kids have no id, with `kids`. */
  anonymous: number;
  parent: Draft | undefined;
  /** Taken out by a remove or a replace. */
  detached: boolean;
}

function draftOf(
  node: Node,
  origin: Node | undefined,
  oldAddress: Address | undefined,
  made: boolean,
  parent: Draft | undefined,
): Draft {
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
    byOp: undefined,
    kids: undefined,
    byId: undefined,
    anonymous: 0,
    parent,
    detached: false,
  };
}

/** A parent's index of its kids by id, built again (the first of two with one id wins). */
function reindex(p: Draft) {
  const byId = new Map<string, Draft>();
  let anonymous = 0;
  for (const k of p.kids!) {
    if (k.id === undefined) anonymous++;
    else if (!byId.has(k.id)) byId.set(k.id, k);
    else Twins.add(p);
  }
  p.byId = byId;
  p.anonymous = anonymous;
}

/** A kid out of its parent's index; a sibling with the same id (a malformed tree) takes its place. */
function unindex(p: Draft, d: Draft) {
  if (d.id === undefined || p.byId!.get(d.id) !== d) return;
  p.byId!.delete(d.id);
  if (!Twins.has(p)) return;
  const twin = p.kids!.find((k) => k.id === d.id);
  if (twin) p.byId!.set(d.id, twin);
}

/** Parents that had two kids with one id when materialized (a malformed tree): only they look for a twin. */
const Twins = new WeakSet<Draft>();

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

/** What the ref index hears of the ops, so it is derived across them (refs.ts, `LiveRefs`). */
export interface WorkListener {
  /** An op set or unset a key of a draft's params. */
  key(d: Draft, key: string): void;
  /** A draft was made (inserted, or a replacement of `old`). */
  made(d: Draft, old?: Draft): void;
}

/** What an inverse op restores, for the inverse's restoring sets (apply.ts). */
export type InverseMeta =
  | { readonly kind: 'key'; readonly draft: Draft; readonly key: string }
  | { readonly kind: 'node'; readonly draft: Draft }
  | { readonly kind: 'move'; readonly from: Address; readonly to: Address }
  | undefined;

/** The draft of a root and the ops applied to it so far. */
export class Work {
  readonly root: { draft: Draft };
  /** Each applied op's inverse, in the order applied (the patch's inverse is them reversed), and what it restores. */
  readonly inverse: Op[] = [];
  readonly inverseMeta: InverseMeta[] = [];
  /** Address changes, in order: what refs follow. */
  readonly moves: { readonly from: Address; readonly to: Address }[] = [];
  /** Subtrees taken out (removed, or a replaced node's kids): where `removed` is read from. */
  private readonly gone: { draft: Draft; cover: Address | undefined }[] = [];
  /** Whether an op inserted, removed, replaced, renamed or moved a node to another parent. */
  structural = false;
  listener?: WorkListener;
  private index = 0;

  constructor(
    readonly old: Node,
    private readonly normalizer?: Normalizer,
  ) {
    this.root = { draft: draftOf(old, old, [], false, undefined) };
  }

  private refuse(code: PatchCode, message: string, ...at: (string | number)[]): never {
    throw new Refused(
      issueOf(PatchIssues, code, message, {
        // the patch tree's path, as patchToNode lays the op out: its kind, or its params by name
        path: ['children', this.index, ...(at[0] === 'kind' ? ['kind'] : at.length ? ['params', ...at] : [])],
        hint: PatchIssues.codes[code].hint,
      }),
    );
  }

  /** A draft's kids, materialized on first use with their index by id (the first of two with one id). */
  kids(d: Draft): Draft[] {
    if (d.kids) return d.kids;
    const children = d.base.children ?? [];
    d.kids = children.map((c) => {
      // the old root's children are the origin's own: a made node's (a replacement's too) are new
      const inner = d.base === d.origin ? c : undefined;
      const oldAddress = inner && d.oldAddress && c.id !== undefined ? [...d.oldAddress, c.id] : undefined;
      return draftOf(c, inner, oldAddress, d.made, d);
    });
    reindex(d);
    return d.kids;
  }

  /** The kid with `id`: O(1) once its parent is materialized (the child index of the base's array before). */
  kid(d: Draft, id: string): Draft | undefined {
    if (!d.kids) {
      const i = childIndex(d.base.children, id);
      return i < 0 ? undefined : this.kids(d)[i];
    }
    return d.byId!.get(id);
  }

  /** The draft at an address (materialized along the path), or undefined. */
  find(address: Address): Draft | undefined {
    let d: Draft | undefined = this.root.draft;
    for (const id of address) {
      d = this.kid(d, id);
      if (!d) return undefined;
    }
    return d;
  }

  /** The draft `tail` (child positions) below `d`, materialized; undefined past the end. */
  below(d: Draft, tail: readonly number[]): Draft | undefined {
    let at: Draft | undefined = d;
    for (const i of tail) {
      at = this.kids(at)[i];
      if (!at) return undefined;
    }
    return at;
  }

  /** Whether a draft is in the tree now (no ancestor taken out). */
  attached(d: Draft): boolean {
    let at: Draft | undefined = d;
    while (at.parent) {
      if (at.detached) return false;
      at = at.parent;
    }
    return !at.detached && at === this.root.draft;
  }

  /** A draft's address now, or undefined (out of the tree, or under a node without an id). */
  addressOf(d: Draft): Address | undefined {
    const ids: string[] = [];
    let at: Draft | undefined = d;
    while (at.parent) {
      if (at.detached || at.id === undefined) return undefined;
      ids.push(at.id);
      at = at.parent;
    }
    return at === this.root.draft && !at.detached ? ids.reverse() : undefined;
  }

  /**
   * The old root's address of the node now at `address`, or undefined when a node made by the patch is there (or
   * none): whether an address still names the node it named (link integrity, a ref's target's identity).
   */
  originAt(address: Address): Address | undefined {
    // a replacement is an edit in place: below it, a node at an id the old node had is that node; an inserted node
    // is none of the old root's
    let d: Draft | undefined = this.root.draft;
    let old: Address | undefined = d.oldAddress;
    let i = 0;
    for (; i < address.length && d.kids; i++) {
      const parent: Draft = d;
      d = d.byId!.get(address[i]);
      if (!d) return undefined;
      if (d.oldAddress) old = d.oldAddress;
      else if (d.made && !parent.made) old = undefined;
      else old = old && [...old, address[i]];
    }
    if (old && i < address.length) old = [...old, ...address.slice(i)];
    return old && resolveAddress(this.old, old) ? old : undefined;
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
        this.refuse('repr.patch.op', 'Not an op: set, unset, insert, remove, move or replace', 'kind');
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

  private touch(d: Draft, where: Where, key: string, op: 'set' | 'unset' | 'rewrite') {
    d.touched ??= {};
    (d.touched[where] ??= new Set()).add(key);
    if (where !== 'params' || op === 'rewrite') return;
    (d.byOp ??= new Set()).add(key);
    if (op === 'set') (d.written ??= new Set()).add(key);
    else d.written?.delete(key);
    this.listener?.key(d, key);
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
    this.push(old === undefined ? { op: 'unset', at, key, where } : { op: 'set', at, key, value: old, where }, {
      kind: 'key',
      draft: d,
      key,
    });
    d[where] = { ...record, [key]: v };
    this.touch(d, where, key, 'set');
  }

  private unset(at: Address, key: string, where: Where = 'params') {
    this.keyOf(key, where);
    const d = this.target(at);
    const record = d[where];
    const old = record ? own(record, key) : undefined;
    if (old === undefined) return;
    this.push({ op: 'set', at, key, value: old, where }, { kind: 'key', draft: d, key });
    const next = { ...record };
    delete next[key];
    d[where] = Object.keys(next).length ? next : undefined;
    this.touch(d, where, key, 'unset');
  }

  private push(op: Op, meta: InverseMeta) {
    this.inverse.push(op);
    this.inverseMeta.push(meta);
  }

  /** Where `place` puts a node among a parent's kids, `skip` (a moved node) out of them; refuses a place naming none. */
  private position(p: Draft, place: Place, skip?: Draft): number {
    const kids = p.kids!;
    const out = skip && skip.parent === p ? 1 : 0;
    if ('last' in place && place.last === true) return kids.length - out;
    const after = (place as { after: string | null }).after;
    if (after === null) return 0;
    const k = p.byId!.get(after);
    if (!k || k === skip) this.refuse('repr.patch.place', 'The place names no sibling', 'after');
    const i = kids.indexOf(k);
    return out && kids.indexOf(skip!) < i ? i : i + 1;
  }

  private addressable(d: Draft, key: string) {
    this.kids(d);
    if (d.anonymous) this.refuse('repr.patch.unaddressable', 'A child of this parent has no id', key);
  }

  private checkNode(node: unknown, key: string): asserts node is Node {
    if (!isNode(node) || !KindPattern.test(node.kind)) this.refuse('repr.patch.op', 'Not a node', key);
  }

  private insert(parent: Address, place: Place, node: Node) {
    const p = this.target(parent, 'parent');
    if (!isPlace(place)) this.refuse('repr.patch.op', 'A place is {after} or {last: true}');
    this.checkNode(node, 'node');
    this.addressable(p, 'parent');
    if (node.id === undefined || !IdPattern.test(node.id))
      this.refuse('repr.patch.id', `An inserted node has an id: ${IdHint}`, 'node');
    if (p.byId!.has(node.id)) this.refuse('repr.patch.id', 'A sibling has this id', 'node');
    const i = this.position(p, place);
    if (parent.length + 1 + nodeHeight(node, MaxDepth) > MaxDepth)
      this.refuse('repr.patch.depth', `The node would sit past ${MaxDepth} levels`, 'node');
    const n = this.normalizer ? this.normalizer.node(node) : node;
    const d = draftOf(n, undefined, undefined, true, p);
    p.kids!.splice(i, 0, d);
    p.byId!.set(node.id, d);
    this.push({ op: 'remove', at: [...parent, node.id] }, undefined);
    this.structural = true;
    this.listener?.made(d);
  }

  private remove(at: Address) {
    if (!isAddress(at)) this.refuse('repr.patch.op', 'An address is a list of ids', 'at');
    if (!at.length) this.refuse('repr.patch.op', 'The root is never removed', 'at');
    const parent = at.slice(0, -1);
    const p = this.target(parent);
    const d = this.kid(p, at[at.length - 1]);
    if (!d) this.refuse('repr.patch.target', 'No node at this address', 'at');
    this.addressable(p, 'at');
    const kids = p.kids!;
    const i = kids.indexOf(d);
    this.push(
      { op: 'insert', parent, place: { after: i > 0 ? kids[i - 1].id! : null }, node: this.finalize(d) },
      { kind: 'node', draft: d },
    );
    kids.splice(i, 1);
    unindex(p, d);
    d.detached = true;
    this.gone.push({ draft: d, cover: undefined });
    this.structural = true;
  }

  private move(at: Address, parent: Address, place: Place, rename: string | undefined) {
    if (!isAddress(at)) this.refuse('repr.patch.op', 'An address is a list of ids', 'at');
    if (!at.length) this.refuse('repr.patch.op', 'The root is never moved', 'at');
    if (!isPlace(place)) this.refuse('repr.patch.op', 'A place is {after} or {last: true}');
    if (rename !== undefined && (typeof rename !== 'string' || !IdPattern.test(rename)))
      this.refuse('repr.patch.id', `A new id is an id: ${IdHint}`, 'rename');
    const from = at.slice(0, -1);
    const src = this.target(from);
    const d = this.kid(src, at[at.length - 1]);
    if (!d) this.refuse('repr.patch.target', 'No node at this address', 'at');
    const dst = this.target(parent, 'parent');
    if (isUnder(parent, at)) this.refuse('repr.patch.cycle', 'A node never moves under itself', 'parent');
    this.addressable(src, 'at');
    this.addressable(dst, 'parent');
    const id = rename ?? d.id!;
    const taken = dst.byId!.get(id);
    if (taken && taken !== d) this.refuse('repr.patch.id', 'A sibling has this id', 'rename');
    const j = this.position(dst, place, d);
    const sameParent = src === dst;
    if (
      !sameParent &&
      parent.length > from.length &&
      parent.length + 1 + nodeHeight(this.finalize(d), MaxDepth) > MaxDepth
    )
      this.refuse('repr.patch.depth', `The node would sit past ${MaxDepth} levels`, 'parent');
    const srcKids = src.kids!;
    const i = srcKids.indexOf(d);
    if (sameParent && id === d.id && i === j) return; // where it is: a no-op
    const to = [...parent, id];
    const changes = !sameParent || id !== d.id;
    this.push(
      {
        op: 'move',
        at: to,
        parent: from,
        place: { after: i > 0 ? srcKids[i - 1].id! : null },
        ...(id !== d.id ? { rename: d.id! } : {}),
      },
      changes ? { kind: 'move', from: to, to: at } : undefined,
    );
    srcKids.splice(i, 1);
    unindex(src, d);
    dst.kids!.splice(j, 0, d);
    d.id = id;
    d.parent = dst;
    dst.byId!.set(id, d);
    if (changes) {
      d.relocated = true;
      this.structural = true;
      this.moves.push({ from: at, to });
    }
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
    this.push({ op: 'replace', at, node: old }, { kind: 'node', draft: d });
    const next = draftOf(n, d.origin, d.oldAddress, true, d.parent);
    next.relocated = d.relocated;
    // what the old node held is made over at its address; a node moved into it from elsewhere is gone
    if (d.kids) for (const k of d.kids) this.gone.push({ draft: k, cover: d.oldAddress });
    d.detached = true;
    if (!d.parent) this.root.draft = next;
    else {
      const kids = d.parent.kids!;
      kids[kids.indexOf(d)] = next;
      if (d.id !== undefined && d.parent.byId!.get(d.id) === d) d.parent.byId!.set(d.id, next);
    }
    this.structural = true;
    this.listener?.made(next, d);
  }

  /** Sets a value at `rel` (`['params', key, …]`) in a draft (a ref rewritten): not written by an op. */
  rewrite(d: Draft, rel: TreePath, value: string) {
    if (rel[0] !== 'params' || typeof rel[1] !== 'string') return;
    const key = rel[1];
    const old = d.params ? own(d.params, key) : undefined;
    if (old === undefined) return;
    const next = setIn(old, rel.slice(2), value);
    if (next === old) return;
    d.params = { ...d.params, [key]: next };
    this.touch(d, 'params', key, 'rewrite');
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
