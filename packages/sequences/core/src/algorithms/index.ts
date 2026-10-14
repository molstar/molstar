/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Aligner } from "./aligner.ts";
import { biomsa } from "./biomsa-aligner.ts";
import { myers } from "./myers.ts";
import { needlemanWunsch } from "./needleman-wunsch.ts";

export {
  type Aligner,
  type AlignerOptions,
  AlignerOptionsError,
  type AlignerParam,
  type AlignerParams,
  alignable,
  cigarFromRows,
  cigarOf,
  detectType,
  pairType,
  type ResolvedOptions,
  readOptions,
  type SequenceType,
} from "./aligner.ts";
export { BiomsaMatrices, type BiomsaMatrix, biomsa } from "./biomsa-aligner.ts";
export { AlignerSizeError, checkCells, MAX_CELLS } from "./limits.ts";
export { type MyersMode, type MyersResult, myers, myersAlign } from "./myers.ts";
export { blosum62, needlemanWunsch } from "./needleman-wunsch.ts";
export {
  type Clusters,
  identityOf,
  Redundancy,
  RedundancySizeError,
  type RowCodes,
  type RowStats,
  redundancy,
  rowStats,
} from "./row-stats.ts";
export {
  alignmentDistances,
  type GuideTree,
  MAX_TREE_ROWS,
  midpointRoot,
  neighborJoining,
  TreeSizeError,
  upgma,
} from "./tree.ts";

/** Every aligner of this package, `biomsa` (the kernel's default, seqspec:format#mappings) first. */
export const Aligners: readonly Aligner[] = [biomsa, needlemanWunsch, myers];
