/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A categorical grid axis (_spec/README.md "Resolved view", Grids: categorical axes; EC-35): a grid section's `of` or
// `of_y` a rows system (an alignment's rows, a matrix's states), each position a named member, in the order the
// section's `member_order` (x) or `member_order_y` (y) gives with the rows layer's values.

import type { Value } from '@molstar/core/repr/tree';
import { memberSystem } from '@molstar/seqx-core/resolve';

/** The members' order: as the data lists them, a tree's leaves, or a field of the rows table. */
export type MemberOrder = 'data' | 'tree' | { readonly field: string; readonly direction: 'ascending' | 'descending' };

/** A categorical axis: the node its members belong to, which of its systems, and their order. */
export interface MemberAxis {
  /** The alignment (its rows) or the matrix (its states). */
  readonly owner: string;
  readonly kind: 'rows' | 'states';
  readonly order: MemberOrder;
}

/** A member order as written (`member_order`, `member_order_y`), `data` when absent or unread. */
export function memberOrderOf(v: Value | undefined): MemberOrder {
  if (v === 'tree') return 'tree';
  if (v && typeof v === 'object' && !Array.isArray(v)) {
    const r = v as Record<string, Value>;
    if (typeof r.field === 'string')
      return { field: r.field, direction: r.direction === 'descending' ? 'descending' : 'ascending' };
  }
  return 'data';
}

/** The categorical axis of a grid's system and its order param, else undefined (a positions axis). */
export function memberAxis(system: string, order: Value | undefined): MemberAxis | undefined {
  const m = memberSystem(system);
  return m && { ...m, order: memberOrderOf(order) };
}
