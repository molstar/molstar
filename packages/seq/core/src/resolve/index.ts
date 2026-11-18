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
  codonRuns,
  codonStatus,
  DefaultGeneticCode,
  type GeneticCode,
  GeneticCodes,
  hasUnsplitExon,
  residuesOf,
  translateCodon,
} from './codon.js';
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
  type MapEntry,
  OpenLength,
  type PathResult,
  type Projection,
  type ResolveConfig,
  type ResolveOptions,
  type Resolver,
  resolveSpec,
  type SourceBinding,
  type SourceParams,
  type StaticModel,
  type SystemEntry,
} from './model.js';
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
export { overlayQueryModel } from './query-model.js';
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
