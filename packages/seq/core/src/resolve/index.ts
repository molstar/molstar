/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// @molstar/seqx-core/resolve: the static model (_spec/README.md "Resolve"): coordinate systems, maps in their
// encodings, paths, projection and closure, sources bound to runtimes, annotations' access plans; the query forms'
// model over a real spec. Named re-exports only.

export {
  codingLength,
  codonLinear,
  codonPairs,
  complementCode,
  codonRuns,
  codonStatus,
  DefaultGeneticCode,
  type GeneticCode,
  GeneticCodes,
  hasUnsplitExon,
  residuesOf,
  translateCodon,
} from './codon.js';
export { type GridData, type PairsData } from './held.js';
export { ResolveIssues } from './issues.js';
export {
  canonicalRuns,
  clipRun,
  composeLinear,
  type Encoding,
  ExpandedLinear,
  expand,
  ImplicitLinear,
  invertRun,
  isEmptyLinear,
  type Linear,
  linearDir,
  linearsEqual,
  type Run,
  RunsLinear,
  runsEqual,
  runTargets,
  scaleLinear,
  toRuns,
  wrapRun,
} from './linear.js';
export {
  type AccessPlan,
  type AnnotationEntry,
  type BoundSourceEntry,
  createResolver,
  DefaultResolveConfig,
  type Lengths,
  type MapEntry,
  type NumberingEntry,
  OpenLength,
  type PathResult,
  type Projection,
  type ResolveConfig,
  type ResolveOptions,
  type Resolver,
  resolveSpec,
  rowsSystemOf,
  type SourceBinding,
  type SourceParams,
  type StaticModel,
  type SystemEntry,
  sourceLengths,
} from './model.js';
export {
  type StructureMapData,
  type StructureStatus,
  structureProvenance,
  structureStatus,
} from './maps-segments.js';
export { viaList } from './maps-via.js';
export type { CompoundEntry, Molecular } from './molecule-query.js';
export type { AtomsEntry, ComplexEntry, ModifiedIndex, MoleculeEntry, UnitEntry } from './molecules.js';
export {
  buildOverlay,
  EmptyOverlay,
  itemAddress,
  type NestedSource,
  type NodeReader,
  nodeAddress,
  type Overlay,
  OverlayBuilder,
  type OverlayNode,
  type OverlayOptions,
} from './overlay.js';
export { overlayQueryModel, type RowsSource } from './query-model.js';
export {
  composeEdges,
  DefaultPathBounds,
  type Edge,
  type ElementStatus,
  findPath,
  identityPath,
  inversePath,
  type Path,
  type PathBounds,
  type PathGraph,
  type PathSearch,
  codingRanges,
  pathDir,
  pathDirection,
  pathRank,
  pathStatus,
  projectBoundary,
  projectEnd,
  projectPoint,
  projectRanges,
  type Step,
} from './paths.js';
