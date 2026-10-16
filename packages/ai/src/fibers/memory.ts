/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { RowConflict, type RowStore, type RowWrite } from './storage.js';

/** Rows kept in this process only: tests, and fibers that should not outlive the page. */
export function createMemoryRows<T extends string = string>(): RowStore<T> {
  const scopes = new Map<string, { generation: number; rows: Map<string, RowWrite<T>> }>();
  const id = (w: RowWrite<T>) => `${w.table}\u0000${w.key}`;
  return {
    async load(scope) {
      const s = scopes.get(scope);
      return s ? { generation: s.generation, rows: [...s.rows.values()] } : undefined;
    },
    async apply(scope, expected, writes) {
      const s = scopes.get(scope) ?? { generation: 0, rows: new Map() };
      if (s.generation !== expected) throw new RowConflict(scope);
      for (const w of writes) {
        if (w.row) s.rows.set(id(w), w);
        else s.rows.delete(id(w));
      }
      s.generation = expected + 1;
      scopes.set(scope, s);
    },
    async drop(scope) {
      scopes.delete(scope);
    },
    async generation(scope) {
      return scopes.get(scope)?.generation ?? 0;
    },
  };
}
