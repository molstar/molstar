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
export { type AssemblyEntry, type AssemblySequence, genomeLoaded, genomeLoadFailure, loadGenome } from './genome.js';
export type { AssemblyNames, IncludeValue, NamedSequence, NameMatch } from './assembly-names.js';
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
  type LinkEnds,
  type BoundSourceEntry,
  createResolver,
  DefaultResolveConfig,
  isSystemInfo,
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
export { type DescribedMap, describedAddress, pendingMaps } from './maps-described.js';
export { structureStatus } from './maps-segments.js';
export { viaList } from './maps-via.js';
export type { CompoundEntry, Molecular } from './molecule-query.js';
export {
  type AtomsEntry,
  /** @internal A molecule's atoms system's address: read by the structure translation (`structure/`). */
  atomsSystemOf,
  type ComplexEntry,
  /** @internal The atoms group map at a path's ends: read by the structure translation (`structure/`). */
  collapseAtoms,
  /** @internal The atoms group map at a path's ends: read by the structure translation (`structure/`). */
  expandAtoms,
  type ModifiedIndex,
  type MoleculeEntry,
  type UnitEntry,
} from './molecules.js';
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
export { openedModel } from './opened.js';
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
  pathStatus,
  projectBoundary,
  projectEnd,
  projectPoint,
  projectRanges,
  type Step,
} from './paths.js';
