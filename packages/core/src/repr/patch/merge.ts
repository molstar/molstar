/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The three-way merge (_spec/README.md "Merge"). Both sides are diffed from the base and each change is keyed
// by the node it touches, an address in the base (diff infers no move across parents and no rename, so every address it
// writes is one). Changes to different nodes, to different keys of one node, or equal on both sides apply; the rest
// conflict by the table of six ops. Theirs is then applied on ours as one patch under the options, so link integrity
// and the validators judge the merged root against ours. Orders are compared as orders over the children on both
// sides (kept, or inserted alike on both), and a child only one side inserted goes after its nearest predecessor on
// that side that still stands.

import { NoChanges } from '../schema/changes.js';
import type { Changes, Language, ValidateOptions } from '../schema/index.js';
import {
  type Address,
  hashNode,
  type Issue,
  issueOf,
  type Node,
  parseAddress,
  resolveAddress,
  valuesEqual,
} from '../tree/index.js';
import { applyPatch } from './apply.js';
import { diff, longestIncreasing } from './diff.js';
import { type InsertOp, type MoveOp, type Op, PatchIssues, type PatchValidator } from './patch.js';
import { parsedRuns, RefPass } from './refs.js';

export interface MergeOptions {
  /** With a language, the merge is normalized, its refs checked for links and validated; without one, it says so. */
  readonly language?: Language;
  /** The slots' languages and the units hook, as validation and normalization take them. */
  readonly validate?: ValidateOptions;
  readonly validators?: readonly PatchValidator[];
  readonly allowDangling?: boolean;
  readonly allowRebind?: boolean;
}

export type MergeResult =
  | {
      readonly ok: true;
      readonly root: Node;
      /** The net change from ours to the root. */
      readonly changes: Changes;
      readonly issues: readonly Issue[];
    }
  | {
      readonly ok: false;
      readonly conflicts: readonly Conflict[];
      /** Ours with every change of theirs independent of every conflict, valid under the options; else ours. */
      readonly partial: Node;
      readonly issues: readonly Issue[];
    };

/**
 * What two changes conflict on. `move` is reserved for a merge of written patches: `merge3` works from `diff`, which
 * infers no move across parents and no rename, so it never gives one.
 */
export type ConflictKind = 'set' | 'remove' | 'insert' | 'move' | 'replace' | 'order' | 'invalid';

export interface Conflict {
  readonly kind: ConflictKind;
  /**
   * In the base: the node removed or replaced, the node set, the inserted node, the reordered parent. For `invalid`, in
   * the merged root: where its first error is (`[]` when it says nowhere).
   */
  readonly at: Address;
  readonly ours: readonly Op[];
  /** Theirs in conflict, then the changes of theirs that depend on them. */
  readonly theirs: readonly Op[];
  /** The issues refusing the merged root, for `invalid`. */
  readonly issues?: readonly Issue[];
}

const key = (a: Address) => a.join('/');

/** The address an op's change is keyed by: its target; an insert's new node; a reorder's parent. */
function touched(op: Op): Address {
  if (op.op === 'insert') return [...op.parent, op.node.id!];
  if (op.op === 'move') return op.parent;
  return op.at;
}

/** One side's diff from the base, keyed by node. */
class Side {
  readonly ops: readonly Op[];
  readonly order = new Map<Op, number>();
  /** Removes and replaces by their target. */
  readonly structural = new Map<string, Op>();
  /** Sets and unsets by target, `where` and key. */
  readonly sets = new Map<string, Op>();
  /** Inserts by the new node's address. */
  readonly inserts = new Map<string, InsertOp>();
  /** Reorders by parent. */
  readonly moves = new Map<string, MoveOp[]>();
  /** Parents of the inserts, by parent. */
  readonly insertsUnder = new Map<string, InsertOp[]>();

  constructor(
    readonly root: Node,
    base: Node,
  ) {
    this.ops = diff(base, root).ops;
    this.ops.forEach((op, i) => {
      this.order.set(op, i);
      switch (op.op) {
        case 'remove':
        case 'replace':
          this.structural.set(key(op.at), op);
          break;
        case 'set':
        case 'unset':
          this.sets.set(setKey(op), op);
          break;
        case 'insert': {
          this.inserts.set(key(touched(op)), op);
          push(this.insertsUnder, key(op.parent), op);
          break;
        }
        case 'move':
          push(this.moves, key(op.parent), op);
          break;
      }
    });
  }
}

const setKey = (op: Op & { op: 'set' | 'unset' }) => `${key(op.at)}\u0000${op.where ?? 'params'}\u0000${op.key}`;

function push<K, V>(map: Map<K, V[]>, k: K, v: V) {
  const list = map.get(k);
  if (list) list.push(v);
  else map.set(k, [v]);
}

interface Open {
  readonly kind: ConflictKind;
  readonly at: Address;
  readonly ours: Set<Op>;
  readonly theirs: Set<Op>;
}

/**
 * Merges what `ours` and `theirs` changed from `base`: the merged root with its net change from ours, or the
 * conflicts and a partial merge.
 */
export function merge3(base: Node, ours: Node, theirs: Node, options: MergeOptions = {}): MergeResult {
  const info: Issue[] = options.language
    ? []
    : [
        issueOf(PatchIssues, 'repr.patch.unvalidated', 'The merge was not validated: no language was given', {
          hint: PatchIssues.codes['repr.patch.unvalidated'].hint,
        }),
      ];
  const o = new Side(ours, base);
  const t = new Side(theirs, base);
  const conflicts = new Map<string, Open>();
  /** The conflict each of theirs' changes left out is listed in. */
  const outOf = new Map<Op, Open>();
  /** Theirs' changes ours already has. */
  const dropped = new Set<Op>();

  const conflict = (kind: ConflictKind, at: Address, ourOp: Op | undefined, theirOp: Op | undefined, k = '') => {
    const id = `${kind}\u0000${key(at)}\u0000${k}`;
    let c = conflicts.get(id);
    if (!c) conflicts.set(id, (c = { kind, at, ours: new Set(), theirs: new Set() }));
    if (ourOp) c.ours.add(ourOp);
    if (theirOp) {
      c.theirs.add(theirOp);
      if (!outOf.has(theirOp)) outOf.set(theirOp, c);
    }
  };

  // remove and replace: a change at or under what the other side removed or replaced (the table's rows, both ways)
  const structural = (op: Op, mine: 'ours' | 'theirs', other: Side) => {
    const at = touched(op);
    for (let n = 0; n <= at.length; n++) {
      const q = n === at.length ? at : at.slice(0, n);
      const r = other.structural.get(key(q));
      if (!r) continue;
      const pair = (kind: ConflictKind) => (mine === 'ours' ? conflict(kind, q, op, r) : conflict(kind, q, r, op));
      const self = n === at.length;
      if (op.op === 'remove') {
        // removed on both sides, or under what the other removed: already in it
        if (r.op === 'remove') {
          if (mine === 'theirs') dropped.add(op);
        } else pair(self ? 'remove' : 'replace');
      } else if (op.op === 'replace' && self && r.op === 'replace') {
        if (hashNode(op.node) === hashNode((r as Op & { op: 'replace' }).node)) {
          if (mine === 'theirs') dropped.add(op);
        } else pair('replace');
      } else pair(r.op === 'remove' ? 'remove' : 'replace');
      return;
    }
  };
  for (const op of o.ops) structural(op, 'ours', t);
  for (const op of t.ops) structural(op, 'theirs', o);
  const open = (op: Op) => !outOf.has(op) && !dropped.has(op);

  // sets: one key of one node; inserts: one id under one parent
  for (const op of t.ops) {
    if (!open(op)) continue;
    if (op.op === 'set' || op.op === 'unset') {
      const mine = o.sets.get(setKey(op));
      if (!mine) continue;
      const same = mine.op === op.op && (op.op === 'unset' || valuesEqual((mine as typeof op).value, op.value));
      if (same) dropped.add(op);
      else conflict('set', op.at, mine, op, `${op.where ?? 'params'}\u0000${op.key}`);
    } else if (op.op === 'insert') {
      const mine = o.inserts.get(key(touched(op)));
      if (!mine) continue;
      if (hashNode(mine.node) === hashNode(op.node)) dropped.add(op);
      else conflict('insert', touched(op), mine, op);
    }
  }

  // orders, compared as orders over the children on both sides, per parent theirs changed the children of
  const plans: Plan[] = [];
  const parents = new Set([...t.moves.keys(), ...t.insertsUnder.keys()]);
  for (const pk of parents) {
    const moves = t.moves.get(pk) ?? [];
    const inserts = t.insertsUnder.get(pk) ?? [];
    if ([...moves, ...inserts].some((op) => outOf.has(op) && outOf.get(op)!.kind !== 'insert')) {
      // the parent is removed or replaced on one side: every change of theirs under it is in that conflict
      continue;
    }
    const at = (moves[0] ?? inserts[0]).parent;
    const b = resolveAddress(base, at)!;
    const x = resolveAddress(ours, at)!;
    const y = resolveAddress(theirs, at)!;
    const ids = (n: Node) => (n.children ?? []).map((c) => c.id!);
    const inX = new Set(ids(x));
    const inY = new Set(ids(y));
    // present on both sides: kept from the base, or inserted on both with equal content (its theirs' insert dropped)
    const both = inserts.filter((op) => dropped.has(op));
    const kept = new Set(ids(b).filter((id) => inX.has(id) && inY.has(id)));
    for (const op of both) kept.add(op.node.id!);
    const kb = ids(b).filter((id) => kept.has(id));
    const kx = ids(x).filter((id) => kept.has(id));
    const ky = ids(y).filter((id) => kept.has(id));
    const inBase = (id: string) => resolveAddress(b, [id]) !== undefined;
    const oursChanged = !sameList(kx.filter(inBase), kb);
    const theirsChanged = !sameList(ky.filter(inBase), kb);
    let takeTheirs = false;
    if (!sameList(kx, ky) && (both.length || (oursChanged && theirsChanged))) {
      // two different orders, or one node both inserted placed two ways: one conflict, ours kept
      for (const m of moves) conflict('order', at, undefined, m);
      for (const m of o.moves.get(pk) ?? []) conflict('order', at, m, undefined);
      for (const op of both) conflict('order', at, o.inserts.get(key(touched(op))), op);
    } else if (theirsChanged && !oursChanged) takeTheirs = true;
    else for (const m of moves) dropped.add(m);
    plans.push({ at, x, y, order: takeTheirs ? ky : kx, takeTheirs, moves, inserts });
  }

  // independence by node: a change of theirs depends on theirs' changes that made what it names
  leaveOutDependents(t, base, theirs, options, outOf);

  const applied = t.ops.filter(open);
  const patch = buildPatch(t, applied, plans, outOf, dropped);
  const listed = [...conflicts.values()];
  const result = patch.length
    ? applyPatch(
        ours,
        { ops: patch },
        {
          ...(options.language ? { language: options.language } : {}),
          ...(options.validate ? { validate: options.validate } : {}),
          ...(options.validators ? { validators: options.validators } : {}),
          ...(options.allowDangling ? { allowDangling: true } : {}),
          ...(options.allowRebind ? { allowRebind: true } : {}),
          maxOps: Number.POSITIVE_INFINITY,
        },
      )
    : undefined;
  const sorted = (ops: Set<Op>, side: Side) => [...ops].sort((a, b) => side.order.get(a)! - side.order.get(b)!);
  const out: Conflict[] = listed.map((c) => ({
    kind: c.kind,
    at: c.at,
    ours: sorted(c.ours, o),
    theirs: sorted(c.theirs, t),
  }));

  if (result && !result.ok) {
    const first = result.issues.find((i) => i.severity === 'error');
    out.push({ kind: 'invalid', at: addressOf(first?.address), ours: [], theirs: applied, issues: result.issues });
    return { ok: false, conflicts: out, partial: ours, issues: info };
  }
  const root = result?.root ?? ours;
  const issues = result ? [...result.issues, ...info] : info;
  if (out.length) return { ok: false, conflicts: out, partial: root, issues };
  return { ok: true, root, changes: result?.changes ?? NoChanges, issues };
}

const sameList = (a: readonly string[], b: readonly string[]) => a.length === b.length && a.every((x, i) => x === b[i]);

/** The address of an issue as ids, or the root when it has none (or has forms). */
function addressOf(text: string | undefined): Address {
  if (!text) return [];
  const p = parseAddress(text).pattern;
  if (!p || p.subtree || p.segments.some((s) => typeof s !== 'string')) return [];
  return p.segments as string[];
}

interface Plan {
  readonly at: Address;
  /** The parent in ours and in theirs. */
  readonly x: Node;
  readonly y: Node;
  /** The children on both sides, in the order taken. */
  readonly order: readonly string[];
  readonly takeTheirs: boolean;
  readonly moves: readonly MoveOp[];
  readonly inserts: readonly InsertOp[];
}

/**
 * Leaves out with a conflict every change of theirs that depends on one in it, transitively: an insert placed after a
 * sibling theirs inserted, and a change writing a ref to what theirs inserted or replaced (with a language), by every id
 * run of the ref, those of the addresses nested in its forms included. A stored patch's refs name nodes of its own dry
 * run, not of theirs: they make no dependency (link integrity still judges them).
 */
function leaveOutDependents(t: Side, base: Node, theirs: Node, options: MergeOptions, outOf: Map<Op, Open>) {
  if (!outOf.size) return;
  // what made each address: theirs' inserts and replaces
  const makers = new Map<string, Op>();
  for (const op of t.ops) if (op.op === 'insert' || op.op === 'replace') makers.set(key(touched(op)), op);
  const dependents = new Map<Op, Op[]>();
  const depend = (op: Op, on: Op | undefined) => {
    if (on && on !== op) push(dependents, on, op);
  };
  const named = (run: Address, op: Op) => {
    for (let n = 1; n <= run.length; n++) depend(op, makers.get(key(run.slice(0, n))));
    depend(op, makers.get(''));
  };
  const pass = options.language ? new RefPass(theirs, options.language, options.validate) : undefined;
  for (const op of t.ops) {
    if (op.op === 'insert') {
      const after = 'after' in op.place ? op.place.after : null;
      if (after !== null) {
        const sib = [...op.parent, after];
        if (!resolveAddress(base, sib)) depend(op, makers.get(key(sib)));
      }
    }
    if (!pass) continue;
    if (op.op === 'insert' || op.op === 'replace')
      pass.walkSubtree(op.node, (_tail, _rel, value, step) => {
        if (!step) for (const run of parsedRuns(value) ?? []) named(run, op);
      });
    else if (op.op === 'set' && (op.where ?? 'params') === 'params') {
      const node = resolveAddress(theirs, op.at);
      if (node)
        pass.walkKey(node.kind, op.key, op.value, (_rel, value, step) => {
          if (!step) for (const run of parsedRuns(value) ?? []) named(run, op);
        });
    }
  }
  const queue = [...outOf.keys()];
  while (queue.length) {
    const op = queue.pop()!;
    const c = outOf.get(op)!;
    for (const d of dependents.get(op) ?? []) {
      if (outOf.has(d)) continue;
      c.theirs.add(d);
      outOf.set(d, c);
      queue.push(d);
    }
  }
}

/**
 * Theirs' changes as one patch on ours: removes, replaces and sets as they are (their addresses stand in ours), then
 * per parent the children in the merged order, by the fewest moves that keep the longest run (as `diff` emits them).
 */
function buildPatch(t: Side, applied: readonly Op[], plans: readonly Plan[], outOf: Map<Op, Open>, dropped: Set<Op>) {
  const ops: Op[] = applied.filter((op) => op.op !== 'insert' && op.op !== 'move');
  const removed = new Set<string>();
  for (const op of applied) if (op.op === 'remove') removed.add(key(op.at));
  for (const plan of plans) {
    const accepted = plan.inserts.filter((op) => !outOf.has(op) && !dropped.has(op));
    if (!accepted.length && !plan.takeTheirs) continue;
    const ids = (n: Node) => (n.children ?? []).map((c) => c.id!);
    const current = ids(plan.x).filter((id) => !removed.has(key([...plan.at, id])));
    const kept = new Set(plan.order);
    const stands = new Set(current);
    // after which sibling each child goes that is not in the kept order; null for first
    const groups = new Map<string | null, string[]>();
    for (let i = 0; i < current.length; i++)
      if (!kept.has(current[i])) push(groups, i > 0 ? current[i - 1] : null, current[i]);
    const byId = new Map<string, InsertOp>();
    for (const op of accepted) byId.set(op.node.id!, op);
    const theirIds = ids(plan.y);
    for (let j = 0; j < theirIds.length; j++) {
      const id = theirIds[j];
      if (!byId.has(id)) continue;
      let k = j - 1;
      while (k >= 0 && !stands.has(theirIds[k])) k--;
      push(groups, k >= 0 ? theirIds[k] : null, id);
      stands.add(id);
    }
    const target: string[] = [];
    const emit = (anchor: string | null) => {
      const stack: string[] = [];
      const g = groups.get(anchor);
      if (g) for (let i = g.length - 1; i >= 0; i--) stack.push(g[i]);
      while (stack.length) {
        const id = stack.pop()!;
        target.push(id);
        const h = groups.get(id);
        if (h) for (let i = h.length - 1; i >= 0; i--) stack.push(h[i]);
      }
    };
    emit(null);
    for (const id of plan.order) {
      target.push(id);
      emit(id);
    }
    if (sameList(target, current)) continue;
    const pos = new Map<string, number>();
    current.forEach((id, i) => pos.set(id, i));
    const inC: number[] = [];
    for (let j = 0; j < target.length; j++) if (pos.has(target[j])) inC.push(j);
    const stay = new Set(longestIncreasing(inC.map((j) => pos.get(target[j])!)).map((k) => inC[k]));
    for (let j = 0; j < target.length; j++) {
      if (stay.has(j)) continue;
      const place = { after: j > 0 ? target[j - 1] : null };
      const id = target[j];
      if (pos.has(id)) ops.push({ op: 'move', at: [...plan.at, id], parent: plan.at, place });
      else ops.push({ ...byId.get(id)!, place });
    }
  }
  return ops;
}
