/**
 * Copyright (c) 2019-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { StructureSelectionQuery } from './query.js';
import { StructureSelectionQueries } from './catalog.js';
import { AminoAcidSelectionQueries, NucleicBaseSelectionQueries } from './residue.js';

export class StructureSelectionQueryRegistry {
  list: StructureSelectionQuery[] = [];
  options: [StructureSelectionQuery, string, string][] = [];
  version = 1;

  add(q: StructureSelectionQuery) {
    this.list.push(q);
    this.options.push([q, q.label, q.category]);
    this.version += 1;
  }

  remove(q: StructureSelectionQuery) {
    const idx = this.list.indexOf(q);
    if (idx !== -1) {
      this.list.splice(idx, 1);
      this.options.splice(idx, 1);
      this.version += 1;
    }
  }

  constructor() {
    // add built-in
    this.list.push(
      ...Object.values(StructureSelectionQueries),
      ...AminoAcidSelectionQueries,
      ...NucleicBaseSelectionQueries,
    );
    this.options.push(...this.list.map((q) => [q, q.label, q.category] as [StructureSelectionQuery, string, string]));
  }
}
