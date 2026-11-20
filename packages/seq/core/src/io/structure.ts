/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// @molstar/seqx-core/io/structure: Seq io (_spec/README.md "Seq io"): a structure's host-neutral description and the
// one converter from it to the nodes molecules are read as, shared by every adapter of a structure (a file's runtime,
// a host's model). `@molstar/core` only. Named re-exports only.

export { StructureIssues } from './issues.js';
export type {
  AtomTable,
  ChainTable,
  ComponentTable,
  CompoundTable,
  EntityTable,
  EntityType,
  IdentifierTable,
  LinkTable,
  ResidueTable,
  StructureDescription,
} from './structure-description.js';
export { type StructureNodes, type StructureNodesOptions, structureNodes } from './structure-nodes.js';
export { DefaultStructureLimits, type StructureLimits } from './structure-read.js';
