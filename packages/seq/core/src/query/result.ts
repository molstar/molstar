/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The three result types (_spec/README.md "Queries", Result types): Nodes (exact addresses of spec nodes and held
// dynamic nodes), Elements (keys of what data makes) and Location (positions on systems). A query gives one of them;
// conversions between them are forms of their own (`extent`, `items`), never implicit. In a selection they travel as
// the tree's `Named`: nodes as nodes, the others as `{ kind: 'other', type, value }`.

import type { Named, Selected } from '@molstar/core/repr/tree';
import type { Elements } from './elements.js';
import type { Location } from './location.js';

export type QueryResult =
  | { readonly type: 'nodes'; readonly nodes: readonly Selected[] }
  | { readonly type: 'elements'; readonly elements: Elements }
  | { readonly type: 'location'; readonly location: Location };

export type QueryResultType = QueryResult['type'];

/** The tree's `Named` of a result. */
export function toNamed(r: QueryResult): Named {
  switch (r.type) {
    case 'nodes':
      return { kind: 'nodes', matches: r.nodes };
    case 'elements':
      return { kind: 'other', type: 'elements', value: r.elements };
    case 'location':
      return { kind: 'other', type: 'location', value: r.location };
  }
}

/** A result of a `Named`, or undefined when another resolver made it (a type the queries do not know). */
export function fromNamed(n: Named): QueryResult | undefined {
  if (n.kind === 'nodes') return { type: 'nodes', nodes: n.matches };
  if (n.type === 'elements') return { type: 'elements', elements: n.value as Elements };
  if (n.type === 'location') return { type: 'location', location: n.value as Location };
  return undefined;
}

/** The result type a `Named` holds, as an issue names it. */
export function namedType(n: Named): string {
  return n.kind === 'nodes' ? 'nodes' : n.type;
}
