/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Node } from '@molstar/core/repr/tree';

/** A spec's section or track (by id) with params added: the story's own composition over an example's. */
export function withParams(spec: Node, kind: 'section' | 'track', id: string, params: Record<string, unknown>): Node {
  const patch = (n: Node): Node =>
    n.kind === kind && n.id === id
      ? { ...n, params: { ...n.params, ...params } as never }
      : n.children
        ? { ...n, children: n.children.map(patch) }
        : n;
  return patch(spec);
}
