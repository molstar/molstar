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

/** A place in the walk: its parent frame, its steps from it, how deep it is. */
interface Frame {
  readonly value: unknown;
  readonly node: boolean;
  readonly depth: number;
  readonly parent: Frame | undefined;
  readonly steps: TreePath;
}

function pathOf(f: Frame): TreePath {
  const parts: (readonly (string | number)[])[] = [];
  for (let x: Frame | undefined = f; x; x = x.parent) parts.push(x.steps);
  const out: (string | number)[] = [];
  for (let i = parts.length - 1; i >= 0; i--) out.push(...parts[i]);
  return out as TreePath;
}

const bound = (code: keyof typeof BoundIssues.codes, message: string, f: Frame): Issue =>
  issueOf(BoundIssues, code, message, { hint: BoundIssues.codes[code].hint, path: pathOf(f) });

function columnBytes(c: unknown): number {
  if (c instanceof CategoryColumn) {
    let n = c.codes.byteLength;
    for (const s of c.dict) n += s.length;
    return n;
  }
  return (c as ArrayBufferView).byteLength;
}

/**
 * The bounds a whole spec must pass, before the description reads it (nothing is allocated past a bound: the walk
 * stops at the first it cannot pass). The issues: `seq.lang.bound.size` or `.depth` alone when the spec breaks one,
 * else the first place of `.column` and `.query_depth`.
 */
export function checkSpecBounds(spec: Node, bounds: SpecBounds = DefaultSpecBounds): Issue[] {
  const issues: Issue[] = [];
  let bytes = 0;
  let column = false;
  let query = false;
  const stack: Frame[] = [{ value: spec, node: true, depth: 0, parent: undefined, steps: [] }];
  /** A scalar is counted in place, no frame made (a list of a million numbers pushes nothing). */
  const leaf = (v: unknown, parent: Frame, ...steps: (string | number)[]): boolean => {
    if (typeof v === 'object' && v !== null) return false;
    if (typeof v !== 'string') bytes += 8;
    else {
      bytes += v.length;
      if (!query && v.charCodeAt(0) === 0x2f && v.includes('(') && queryDepth(v) > bounds.maxQueryDepth) {
        query = true;
        const at: Frame = { value: v, node: false, depth: parent.depth + 1, parent, steps: steps as TreePath };
        issues.push(bound('seq.lang.bound.query_depth', `A query nests deeper than ${bounds.maxQueryDepth} forms`, at));
      }
    }
    return true;
  };
  while (stack.length) {
    const f = stack.pop()!;
    if (f.depth > bounds.maxDepth) {
      issues.push(bound('seq.lang.bound.depth', `The spec nests deeper than ${bounds.maxDepth} levels`, f));
      return issues;
    }
    const v = f.value;
    if (f.node) {
      const n = v as Node;
      bytes += (n.kind?.length ?? 0) + (n.id?.length ?? 0);
      const kids = n.children;
      if (kids)
        for (let i = kids.length - 1; i >= 0; i--)
          stack.push({ value: kids[i], node: true, depth: f.depth + 1, parent: f, steps: ['children', i] });
      const params = n.params;
      if (params)
        for (const k of Object.keys(params)) {
          bytes += k.length;
          if (!leaf(params[k], f, 'params', k))
            stack.push({ value: params[k], node: false, depth: f.depth + 1, parent: f, steps: ['params', k] });
        }
    } else if (typeof v === 'object' && v !== null) {
      if (isColumn(v)) {
        if (v.length > bounds.maxColumnLength && !column) {
          column = true;
          issues.push(bound('seq.lang.bound.column', `A column is longer than ${bounds.maxColumnLength}`, f));
        }
        bytes += columnBytes(v);
      } else if (ArrayBuffer.isView(v)) bytes += v.byteLength;
      else if (Array.isArray(v)) {
        bytes += 8;
        for (let i = v.length - 1; i >= 0; i--)
          if (!leaf(v[i], f, i)) stack.push({ value: v[i], node: false, depth: f.depth + 1, parent: f, steps: [i] });
      } else {
        for (const k of Object.keys(v)) {
          bytes += k.length;
          const x = (v as Record<string, unknown>)[k];
          if (!leaf(x, f, k)) stack.push({ value: x, node: false, depth: f.depth + 1, parent: f, steps: [k] });
        }
      }
    } else bytes += 8;
    if (bytes > bounds.maxBytes) {
      issues.push(bound('seq.lang.bound.size', `The spec is larger than ${bounds.maxBytes} bytes decoded`, f));
      return issues;
    }
  }
  return issues;
}

/**
 * The nesting bounds of a subtree a patch made or moved, `at` levels below the root (an address's length): the first
 * place deeper than `maxDepth`, or the first query deeper than `maxQueryDepth`, as `checkSpecBounds` reports them, for
 * the cost of the subtree alone (`children` false: the node's own params only).
 */
export function checkSubtreeDepth(
  node: Node,
  path: TreePath,
  at: number,
  bounds: SpecBounds = DefaultSpecBounds,
  children = true,
): Issue[] {
  const issues: Issue[] = [];
  const stack: Frame[] = [{ value: node, node: true, depth: at, parent: undefined, steps: path }];
  let query = false;
  const leaf = (v: unknown, parent: Frame, ...steps: (string | number)[]): boolean => {
    if (typeof v === 'object' && v !== null) return false;
    if (
      typeof v === 'string' &&
      !query &&
      v.charCodeAt(0) === 0x2f &&
      v.includes('(') &&
      queryDepth(v) > bounds.maxQueryDepth
    ) {
      query = true;
      const at: Frame = { value: v, node: false, depth: parent.depth + 1, parent, steps: steps as TreePath };
      issues.push(bound('seq.lang.bound.query_depth', `A query nests deeper than ${bounds.maxQueryDepth} forms`, at));
    }
    return true;
  };
  while (stack.length) {
    const f = stack.pop()!;
    if (f.depth > bounds.maxDepth) {
      issues.push(bound('seq.lang.bound.depth', `The spec nests deeper than ${bounds.maxDepth} levels`, f));
      return issues;
    }
    const v = f.value;
    if (f.node) {
      const n = v as Node;
      for (let i = children || f.parent ? (n.children?.length ?? 0) - 1 : -1; i >= 0; i--)
        stack.push({ value: n.children![i], node: true, depth: f.depth + 1, parent: f, steps: ['children', i] });
      const params = n.params ?? {};
      for (const k of Object.keys(params)) {
        if (!leaf(params[k], f, 'params', k))
          stack.push({ value: params[k], node: false, depth: f.depth + 1, parent: f, steps: ['params', k] });
      }
    } else if (typeof v === 'object' && v !== null && !ArrayBuffer.isView(v) && !(v instanceof CategoryColumn)) {
      if (Array.isArray(v))
        for (let i = v.length - 1; i >= 0; i--) {
          if (!leaf(v[i], f, i)) stack.push({ value: v[i], node: false, depth: f.depth + 1, parent: f, steps: [i] });
        }
      else
        for (const k of Object.keys(v)) {
          const x = (v as Record<string, unknown>)[k];
          if (!leaf(x, f, k)) stack.push({ value: x, node: false, depth: f.depth + 1, parent: f, steps: [k] });
        }
    }
  }
  return issues;
}

/**
 * What a reporter keeps: the first `maxPlaces` issues of a code, at most `maxNodeIssues` on a node; the rest counted.
 * `summary()` is one `seq.lang.bound.issues` a code that lost any.
 */
export class IssueCaps {
  private readonly perCode = new Map<string, number>();
  private readonly perNode = new WeakMap<object, number>();
  private readonly lost = new Map<string, number>();
  constructor(private readonly bounds: SpecBounds = DefaultSpecBounds) {}

  /** Whether an issue of `code` on `node` is kept; counts it when it is not. */
  keep(code: string, node: object): boolean {
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
