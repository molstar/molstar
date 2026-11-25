/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The three result types (_spec/README.md "Queries", Result types): Nodes (exact addresses of spec nodes and held
// dynamic nodes), Keys (keys of what data makes) and Location (positions on systems). A query gives one of them;
// conversions between them are forms of their own (`extent`, `items`), never implicit. In a selection they travel as
// the tree's `Named`: nodes as nodes, the others as `{ kind: 'other', type, value }`.

import type { Named, Selected } from '@molstar/core/repr/tree';
import type { Keys } from './keys.js';
import type { Location } from './location.js';

export type QueryResult =
  | { readonly type: 'nodes'; readonly nodes: readonly Selected[] }
  | { readonly type: 'keys'; readonly keys: Keys }
  | { readonly type: 'location'; readonly location: Location };

export type QueryResultType = QueryResult['type'];

/** The tree's `Named` of a result. */
export function toNamed(r: QueryResult): Named {
  switch (r.type) {
    case 'nodes':
      return { kind: 'nodes', matches: r.nodes };
    case 'keys':
      return { kind: 'other', type: 'keys', value: r.keys };
    case 'location':
      return { kind: 'other', type: 'location', value: r.location };
  }
}

const isObject = (v: unknown): v is Record<string, unknown> => typeof v === 'object' && v !== null;

/** Whether a value has a Location's shape: parts and rects as arrays, each part a system and an Int32Array. */
function isLocation(v: unknown): v is Location {
  if (!isObject(v) || !Array.isArray(v.parts) || !Array.isArray(v.rects)) return false;
  const part = (p: unknown) => isObject(p) && typeof p.system === 'string' && p.ranges instanceof Int32Array;
  return v.parts.every(part) && v.rects.every((r: unknown) => isObject(r) && part(r.a) && part(r.b));
}

/** Whether a value has Keys' shape: keys, each an owner and a string or number key. */
function isKeys(v: unknown): v is Keys {
  if (!isObject(v) || !Array.isArray(v.keys)) return false;
  return v.keys.every(
    (k: unknown) =>
      isObject(k) && typeof k.owner === 'string' && (typeof k.key === 'string' || typeof k.key === 'number'),
  );
}

/**
 * A result of a `Named`, or undefined when another resolver made it: a type the queries do not know, or a value
 * whose shape is not the type's (a resolver of another language using the same type name).
 */
export function fromNamed(n: Named): QueryResult | undefined {
  if (n.kind === 'nodes') return { type: 'nodes', nodes: n.matches };
  if (n.type === 'keys' && isKeys(n.value)) return { type: 'keys', keys: n.value };
  if (n.type === 'location' && isLocation(n.value)) return { type: 'location', location: n.value };
  return undefined;
}

/** The result type a `Named` holds, as an issue names it. */
export function namedType(n: Named): string {
  return n.kind === 'nodes' ? 'nodes' : n.type;
}
