/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A patch and its operations (_spec/README.md "Patches"): revisions (immutable roots, named by their hash) are the
// truth, and a patch is the form a change travels in, derived by `diff` or written by an agent, the wire or a stored
// story step. Six operations by addresses of ids; applied all or nothing, with the inverse of what was applied.

import type { Changes, Language, ValidateOptions, Where } from '../schema/index.js';
import { type Address, defineIssues, type Issue, type Node, type Value } from '../tree/index.js';

export type { Changes, Where };

/** Where a node goes among its new siblings: after one (by id), first (`after: null`), or last. */
export type Place = { readonly after: string | null } | { readonly last: true };

export interface SetOp {
  readonly op: 'set';
  readonly at: Address;
  readonly key: string;
  readonly value: Value;
  readonly where?: Where;
}
export interface UnsetOp {
  readonly op: 'unset';
  readonly at: Address;
  readonly key: string;
  readonly where?: Where;
}
export interface InsertOp {
  readonly op: 'insert';
  readonly parent: Address;
  readonly place: Place;
  readonly node: Node;
}
export interface RemoveOp {
  readonly op: 'remove';
  readonly at: Address;
}
export interface MoveOp {
  readonly op: 'move';
  readonly at: Address;
  readonly parent: Address;
  readonly place: Place;
  readonly rename?: string;
}
export interface ReplaceOp {
  readonly op: 'replace';
  readonly at: Address;
  readonly node: Node;
}
export type Op = SetOp | UnsetOp | InsertOp | RemoveOp | MoveOp | ReplaceOp;

export interface Patch {
  readonly ops: readonly Op[];
  /** The revision it was made against: refused on any other, unless the root is already `result`. */
  readonly base?: string;
  /** The revision it gives, normalized, when the writer knows it: a mismatch after applying refuses. */
  readonly result?: string;
}

/** A language's own check of a change: what the change introduced, in proportion to `changes`. */
export type PatchValidator = (before: Node, after: Node, changes: Changes) => readonly Issue[];

export interface ApplyOptions {
  /** With a language, a patch is normalized, its refs follow moves, and it is validated and checked for links. */
  readonly language?: Language;
  /** The slots' languages and the units hook, as validation and normalization take them. */
  readonly validate?: ValidateOptions;
  readonly validators?: readonly PatchValidator[];
  /** At most this many ops (default 10,000); a patch from `diff` is applied with `Infinity`. */
  readonly maxOps?: number;
  /** Refuse a patch without a `base` (default false). */
  readonly requireBase?: boolean;
  /** A ref left dangling is a warning, not a refusal (default false). */
  readonly allowDangling?: boolean;
  /** A dangling ref that comes to name a node is allowed (default false). */
  readonly allowRebind?: boolean;
}

export type PatchResult =
  | {
      readonly ok: true;
      readonly changed: boolean;
      readonly root: Node;
      readonly revision: string;
      /** Absent when nothing changed. */
      readonly inverse?: Patch;
      /** What applying the inverse needs: the allowances this patch used, turned around. */
      readonly inverseNeeds?: { readonly allowDangling?: true; readonly allowRebind?: true };
      /** The net difference, ref rewrites included as `params`. */
      readonly changes: Changes;
      /** Said once when there was no language to rewrite refs with. */
      readonly refs?: 'not-rewritten';
      /** Warnings and infos. */
      readonly issues: readonly Issue[];
    }
  | { readonly ok: false; readonly issues: readonly Issue[] };

/** The default bound on a written patch's ops. */
export const MaxOps = 10_000;

export const PatchIssues = defineIssues('repr.patch', {
  'repr.patch.stale': {
    severity: 'error',
    description: 'The patch was made against another revision than the root (its `base`)',
    hint: 'Make the patch again against the current revision, or rebase it on purpose',
  },
  'repr.patch.base': {
    severity: 'error',
    description: 'The patch has no `base`, and the caller requires one',
    hint: 'Give the revision the patch was made against',
  },
  'repr.patch.result': {
    severity: 'error',
    description: 'The patch gave another revision than its `result`',
    hint: 'Leave `result` out when the writer cannot normalize',
  },
  'repr.patch.target': {
    severity: 'error',
    description: 'An op names a node the tree does not have (its `at` or `parent`)',
  },
  'repr.patch.id': {
    severity: 'error',
    description:
      'An inserted or renamed node has no id, a malformed one or one a sibling has, or a replacement has another id',
    hint: 'An id is unique among its siblings; a replacement keeps the id of what it replaces',
  },
  'repr.patch.place': {
    severity: 'error',
    description: 'A place names no sibling (its `after`)',
  },
  'repr.patch.cycle': {
    severity: 'error',
    description: 'A move puts a node under itself',
  },
  'repr.patch.key': {
    severity: 'error',
    description: 'A `set` or `unset` names a reserved key (`kind`, `__proto__`, `$…`)',
  },
  'repr.patch.unaddressable': {
    severity: 'error',
    description: 'An insert, remove or move under a parent with a child that has no id',
    hint: 'Change that parent by `replace`',
  },
  'repr.patch.op': {
    severity: 'error',
    description: 'Not an op: an unknown op, a field missing or malformed, or a remove or move of the root',
  },
  'repr.patch.limit': {
    severity: 'error',
    description: 'More ops than the caller allows (`maxOps`)',
  },
  'repr.patch.depth': {
    severity: 'error',
    description: 'An op would make the tree, or a value in it, deeper than MaxDepth levels',
  },
  'repr.patch.dangling': {
    severity: 'error',
    description:
      'A ref names no node after the patch, and named one before it or was written by it (a warning with `allowDangling`)',
    hint: 'Give the ref another target, or allow it to dangle',
  },
  'repr.patch.rebound': {
    severity: 'error',
    description: 'A ref the patch did not change named no node before it and names one after it',
    hint: 'Give the ref another target, or allow the rebinding',
  },
});

export type PatchCode = keyof typeof PatchIssues.codes;
