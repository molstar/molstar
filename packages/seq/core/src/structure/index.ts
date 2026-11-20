/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// @molstar/seqx-core/structure (_spec/README.md "Structures"): the translation between Seq*'s Locations and a
// structure's units, for a host wiring a structure viewer beside a sequence view. A library, not a port: plain data in
// and out, no structure viewer's types, nothing registered. Named re-exports only.

export {
  type AuthNumbering,
  type ResidueRows,
  type StructureFacts,
  structureFacts,
  type UnitFacts,
  type UnitKind,
} from './facts.js';
export { StructureSyncIssues } from './issues.js';
export type { StructureBinding, StructureModel } from './model.js';
export {
  DefaultSyncLimits,
  type FromUnitsOptions,
  type FromUnitsResult,
  fromUnits,
  type StructureTarget,
  type SyncLimits,
  type ToUnitsOptions,
  toUnits,
  type UnitsResult,
} from './units.js';
