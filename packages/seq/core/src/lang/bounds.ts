/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { CategoryColumn, type Issue, isColumn, issueOf, type Node, type TreePath } from '@molstar/core/repr/tree';
import { BoundIssues } from './issues.js';

// The bounds on what a spec sizes (_spec/README.md, Language, Bounds): a spec comes from a URL, a link or an agent, so
// the work its size asks for is refused with a coded issue before it starts. `checkSpecBounds` is the first gate of
// `checkSeqSpec`: one iterative walk that counts and allocates no tree (a 20,000-deep spec costs no stack), stops at the
// first bound a whole spec cannot pass (its size, its depth) and reports the first place of the others. A patch's
// incremental check (`checkSeqChanges`) refuses the same nesting depth through `checkSubtreeDepth`.

/** The bounds, each a configuration item with its default (`DefaultSpecBounds`). */
export interface SpecBounds {
  /** A spec's decoded size: strings by their length, numbers 8 bytes, columns by their bytes: 64 MB. */
  readonly maxBytes: number;
  /** Levels of nodes and of the values nested in their params: 256. */
  readonly maxDepth: number;
  /** Levels of nested query forms (`(kind :key (kind …))`) in an address: 64. */
  readonly maxQueryDepth: number;
  /** A column's length: 2^31 − 1 (and the spec's size above). */
  readonly maxColumnLength: number;
  /** Issues kept a node: 100. */
  readonly maxNodeIssues: number;
  /** Places listed a code, the first ones: 10. */
  readonly maxPlaces: number;
}

export const DefaultSpecBounds: SpecBounds = {
  maxBytes: 64 * 1024 * 1024,
  maxDepth: 256,
  maxQueryDepth: 64,
  maxColumnLength: 0x7fffffff,
  maxNodeIssues: 100,
  maxPlaces: 10,
};

/** The defaults with `over` laid on them. */
export function specBounds(over: Partial<SpecBounds> = {}): SpecBounds {
  return { ...DefaultSpecBounds, ...over };
}

/** The deepest nesting of `(` in a string outside double quotes (a query form's depth). */
export function queryDepth(s: string): number {
  let depth = 0;
  let max = 0;
  let quoted = false;
  for (let i = 0; i < s.length; i++) {
    const c = s.charCodeAt(i);
    if (quoted) {
      if (c === 0x5c) i++;
      else if (c === 0x22) quoted = false;
    } else if (c === 0x22) quoted = true;
    else if (c === 0x28) {
      if (++depth > max) max = depth;
    } else if (c === 0x29 && depth > 0) depth--;
  }
  return max;
}

/** What the size estimate charges for each node, array, record, key and element: about JSON's punctuation. */
export const ValueOverhead = 8;

/** A place in the walk: its parent frame and the steps from it (built only when an issue names the place). */
interface Frame {
  readonly value: unknown;
  readonly node: boolean;
  readonly depth: number;
  readonly parent: Frame | undefined;
  readonly a: string | number | undefined;
  readonly b: string | number | undefined;
  /** Whether strings here are addresses (queries are counted); not in `custom`. */
  readonly query: boolean;
}

interface WalkOptions {
  readonly bounds: SpecBounds;
  /** The root's step path and its level (an address's length). */
  readonly path: TreePath;
  readonly at: number;
  /** Walk the root's children (false: the node's own params, custom and args). */
  readonly children: boolean;
  /** Only these keys of the root's `params` or `custom`. */
  readonly keys?: { readonly where: 'params' | 'custom'; readonly keys: readonly string[] };
  /** The bytes counted so far, shared by walks that add up (a patch's); a fresh count when absent. */
  readonly budget?: { bytes: number };
}

function pathOf(f: Frame, root: TreePath): TreePath {
  const out: (string | number)[] = [];
  for (let x: Frame | undefined = f; x; x = x.parent) {
    if (x.parent === undefined) out.push(...[...root].reverse());
    else {
      if (x.b !== undefined) out.push(x.b);
      if (x.a !== undefined) out.push(x.a);
    }
  }
  return out.reverse() as unknown as TreePath;
}

const isPlain = (v: unknown): v is Record<string, unknown> =>
  typeof v === 'object' && v !== null && !Array.isArray(v) && !ArrayBuffer.isView(v) && !(v instanceof CategoryColumn);

function columnBytes(c: unknown): number {
  if (c instanceof CategoryColumn) {
    let n = c.codes.byteLength;
    for (const s of c.dict) n += s.length;
    return n;
  }
  return (c as ArrayBufferView).byteLength;
}

/**
 * One iterative walk over a node, its children, params, custom data and args, as values alike: the depth, the bytes
 * (a fixed overhead for each node, record, array, key and element, charged before the children are pushed or
 * visited, so a million empty records cost what they say), a column's length, a query's depth in params and args;
 * the first fatal one (size, depth, shape) ends it.
 */
function walk(root: unknown, o: WalkOptions): { issues: Issue[]; bytes: number } {
  const { bounds } = o;
  const issues: Issue[] = [];
  const budget = o.budget ?? { bytes: 0 };
  let column = false;
  let query = false;
  const stack: Frame[] = [
    { value: root, node: true, depth: o.at, parent: undefined, a: undefined, b: undefined, query: true },
  ];
  const say = (code: keyof typeof BoundIssues.codes, message: string, f: Frame) =>
    issues.push(issueOf(BoundIssues, code, message, { hint: BoundIssues.codes[code].hint, path: pathOf(f, o.path) }));
  const sizeIssue = (f: Frame) => {
    say('seq.lang.bound.size', `The spec is larger than ${bounds.maxBytes} bytes decoded`, f);
    return { issues, bytes: budget.bytes };
  };
  const frame = (
    v: unknown,
    f: Frame,
    a: string | number | undefined,
    b?: string | number,
    node = false,
    q = f.query,
  ): Frame => ({ value: v, node, depth: f.depth + 1, parent: f, a, b, query: q });
  /** A scalar is counted in place, no frame made. */
  const leaf = (v: unknown, f: Frame, a: string | number): boolean => {
    if (typeof v === 'object' && v !== null) return false;
    if (typeof v === 'string') {
      budget.bytes += v.length;
      if (f.query && !query && v.charCodeAt(0) === 0x2f && v.includes('(') && queryDepth(v) > bounds.maxQueryDepth) {
        query = true;
        say('seq.lang.bound.query_depth', `A query nests deeper than ${bounds.maxQueryDepth} forms`, frame(v, f, a));
      }
    } else if (typeof v === 'bigint') budget.bytes += v.toString().length;
    else budget.bytes += 8;
    return true;
  };
  const fatal = (code: keyof typeof BoundIssues.codes, message: string, f: Frame) => {
    say(code, message, f);
    return { issues, bytes: budget.bytes };
  };
  while (stack.length) {
    const f = stack.pop()!;
    if (f.depth > bounds.maxDepth)
      return fatal('seq.lang.bound.depth', `The spec nests deeper than ${bounds.maxDepth} levels`, f);
    const v = f.value;
    if (f.node) {
      if (!isPlain(v) || typeof v.kind !== 'string')
        return fatal('seq.lang.bound.shape', 'A node is not a record with a kind', f);
      const n = v as unknown as Node;
      budget.bytes += ValueOverhead + n.kind.length + (typeof n.id === 'string' ? n.id.length : 0);
      if (o.keys && f.parent === undefined) {
        const rec = n[o.keys.where];
        if (rec !== undefined) {
          if (!isPlain(rec))
            return fatal('seq.lang.bound.shape', `\`${o.keys.where}\` is not a record`, frame(rec, f, o.keys.where));
          const holder = frame(rec, f, o.keys.where, undefined, false, o.keys.where === 'params');
          for (const k of o.keys.keys) {
            const x = rec[k];
            if (x === undefined) continue;
            budget.bytes += ValueOverhead + k.length;
            if (!leaf(x, holder, k)) stack.push(frame(x, holder, k));
          }
        }
      } else {
        if (o.children || f.parent) {
          const kids = n.children;
          if (kids !== undefined) {
            if (!Array.isArray(kids))
              return fatal('seq.lang.bound.shape', '`children` is not a list', frame(kids, f, 'children'));
            budget.bytes += kids.length * ValueOverhead;
            if (budget.bytes > bounds.maxBytes) return sizeIssue(f);
            for (let i = kids.length - 1; i >= 0; i--) stack.push(frame(kids[i], f, 'children', i, true));
          }
        }
        const { params, custom, args } = n;
        if (params !== undefined) stack.push(frame(params, f, 'params'));
        if (custom !== undefined) stack.push(frame(custom, f, 'custom', undefined, false, false));
        if (args !== undefined) stack.push(frame(args, f, 'args'));
      }
    } else if (typeof v === 'object' && v !== null) {
      if (isColumn(v)) {
        if (v.length > bounds.maxColumnLength && !column) {
          column = true;
          say('seq.lang.bound.column', `A column is longer than ${bounds.maxColumnLength}`, f);
        }
        budget.bytes += columnBytes(v);
      } else if (ArrayBuffer.isView(v)) budget.bytes += v.byteLength;
      else if (Array.isArray(v)) {
        budget.bytes += (v.length + 1) * ValueOverhead;
        if (budget.bytes > bounds.maxBytes) return sizeIssue(f);
        for (let i = v.length - 1; i >= 0; i--) {
          if (!leaf(v[i], f, i)) stack.push(frame(v[i], f, i));
          if (budget.bytes > bounds.maxBytes) return sizeIssue(f);
        }
      } else {
        budget.bytes += ValueOverhead;
        for (const k in v) {
          if (!Object.hasOwn(v, k)) continue;
          budget.bytes += ValueOverhead + k.length;
          const x = (v as Record<string, unknown>)[k];
          if (!leaf(x, f, k)) stack.push(frame(x, f, k));
          if (budget.bytes > bounds.maxBytes) return sizeIssue(f);
        }
      }
    }
    if (budget.bytes > bounds.maxBytes) return sizeIssue(f);
  }
  return { issues, bytes: budget.bytes };
}

/**
 * The bounds a whole spec must pass, before the description reads it (nothing is allocated past a bound: the walk
 * stops at the first it cannot pass). The issues: `seq.lang.bound.size`, `.depth` or `.shape` alone when the spec
 * breaks one, else the first place of `.column` and `.query_depth`.
 */
export function checkSpecBounds(spec: Node, bounds: SpecBounds = DefaultSpecBounds): Issue[] {
  return walk(spec, { bounds, path: [], at: 0, children: true }).issues;
}

/**
 * What a node weighs by the size estimate (and the issues of its walk), for what is held by weight (a dynamic node of
 * an answer): stops once past `bounds.maxBytes`.
 */
export function measureNode(node: unknown, bounds: SpecBounds = DefaultSpecBounds): { bytes: number; issues: Issue[] } {
  return walk(node, { bounds, path: [], at: 0, children: true });
}

/**
 * The bounds of a subtree a patch made or moved, `at` levels below the root (an address's length): the first place
 * deeper than `maxDepth`, or the first query deeper than `maxQueryDepth`, as `checkSpecBounds` reports them, for the
 * cost of the subtree alone. `children` false walks the node's own params, custom data and args; `keys` only those of
 * its `params` or `custom` a change set. `budget` adds the bytes of this walk to a total kept by the caller (a patch's:
 * `maxBytes` bounds what one patch makes and sets).
 */
export function checkSubtreeDepth(
  node: Node,
  path: TreePath,
  at: number,
  bounds: SpecBounds = DefaultSpecBounds,
  children = true,
  extra: Pick<WalkOptions, 'keys' | 'budget'> = {},
): Issue[] {
  return walk(node, { bounds, path, at, children, ...extra }).issues;
}

/**
 * What a reporter keeps: the first `maxPlaces` issues of a code, at most `maxNodeIssues` on a node; the rest counted.
 * `summary()` is one `seq.lang.bound.issues` a code that lost any.
 */
export class IssueCaps {
  private readonly perCode = new Map<string, number>();
  private readonly perNode = new Map<object | string, number>();
  private readonly lost = new Map<string, number>();
  constructor(private readonly bounds: SpecBounds = DefaultSpecBounds) {}

  /** Whether an issue of `code` on `node` is kept; counts it when it is not. */
  keep(code: string, node: object | string): boolean {
    const byCode = this.perCode.get(code) ?? 0;
    const byNode = this.perNode.get(node) ?? 0;
    if (byCode >= this.bounds.maxPlaces || byNode >= this.bounds.maxNodeIssues) {
      this.lost.set(code, (this.lost.get(code) ?? 0) + 1);
      return false;
    }
    this.perCode.set(code, byCode + 1);
    this.perNode.set(node, byNode + 1);
    return true;
  }

  summary(): Issue[] {
    const out: Issue[] = [];
    for (const [code, n] of this.lost)
      out.push(
        issueOf(BoundIssues, 'seq.lang.bound.issues', `${n} more issues of ${code} are not listed`, {
          hint: BoundIssues.codes['seq.lang.bound.issues'].hint,
        }),
      );
    return out;
  }
}
