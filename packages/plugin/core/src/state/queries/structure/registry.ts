/**
 * Copyright (c) 2019-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { StructureSelectionQuery } from './query.js';

export class StructureSelectionQueryRegistry {
  list: StructureSelectionQuery[] = [];
  options: [StructureSelectionQuery, string, string][] = [];
  version = 1;
  /** Registration counts by query identity */
  private counts = new Map<StructureSelectionQuery, number>();

  /** Queries have no key, so registration never conflicts. */
  findConflict(_q: StructureSelectionQuery): string | undefined {
    return undefined;
  }

  /** Registering the same query again increments its count. */
  add(q: StructureSelectionQuery) {
    const count = this.counts.get(q);
    if (count !== undefined) {
      this.counts.set(q, count + 1);
      return;
    }

    this.counts.set(q, 1);
    this.list.push(q);
    this.options.push([q, q.label, q.category]);
    this.version += 1;
  }

  /** Decrements the count and removes the query at zero; no-op for an unknown query. */
  remove(q: StructureSelectionQuery) {
    const count = this.counts.get(q);
    if (count === undefined) return;
    if (count > 1) {
      this.counts.set(q, count - 1);
      return;
    }

    this.counts.delete(q);
    const idx = this.list.indexOf(q);
    if (idx !== -1) {
      this.list.splice(idx, 1);
      this.options.splice(idx, 1);
      this.version += 1;
    }
  }
}
