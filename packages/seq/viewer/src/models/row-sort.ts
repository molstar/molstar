/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { derive, type ReadonlyState } from '@molstar/core/reactive';
import type { Node } from '@molstar/core/repr/tree';
import type { StaticModel } from '@molstar/seqx-core/resolve';

/** A rows layer the spec orders by a column (`SeqViewer.rowSort`). */
export interface RowSortEntry {
  /** The layer's address. */
  readonly layer: string;
  /** The track's address. */
  readonly track: string;
  readonly field: string;
  readonly direction: 'ascending' | 'descending';
}

const None: readonly RowSortEntry[] = [];

/** An `order` param as a sort by a column (not `data` or `tree`), or undefined. */
export function sortOf(order: unknown): { field: string; direction: 'ascending' | 'descending' } | undefined {
  if (typeof order !== 'object' || order === null) return undefined;
  const o = order as { field?: unknown; direction?: unknown };
  return typeof o.field === 'string'
    ? { field: o.field, direction: o.direction === 'descending' ? 'descending' : 'ascending' }
    : undefined;
}

/** Every rows layer of a spec ordered by a column, in the spec's order. */
export function rowSortOf(root: Node): readonly RowSortEntry[] {
  let out: RowSortEntry[] | undefined;
  const walk = (node: Node, path: readonly string[]) => {
    for (const c of node.children ?? []) {
      if (c.id === undefined) continue;
      if (c.kind === 'layer') {
        const s = c.params?.type === 'rows' ? sortOf(c.params?.order) : undefined;
        if (s) (out ??= []).push({ layer: `/${[...path, c.id].join('/')}`, track: `/${path.join('/')}`, ...s });
      } else walk(c, [...path, c.id]);
    }
  };
  walk(root, []);
  return out ?? None;
}

function same(a: readonly RowSortEntry[], b: readonly RowSortEntry[]): boolean {
  return (
    a === b ||
    (a.length === b.length &&
      a.every((x, i) => x.layer === b[i].layer && x.field === b[i].field && x.direction === b[i].direction))
  );
}

/** The sorts of the spec shown (_spec/README.md "The rows' commands"): the application's to list and clear. */
export function rowSortState(model: ReadonlyState<StaticModel>): ReadonlyState<readonly RowSortEntry[]> {
  return derive([model], (m) => rowSortOf(m.spec), same);
}
