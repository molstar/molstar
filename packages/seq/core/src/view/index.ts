/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// @molstar/seqx-core/view: the resolved view (_spec/README.md "Resolved view"), the public, renderer-independent form
// other renderers read. Named re-exports only.

export {
  type AxisShape,
  axesEqual,
  axisBreaks,
  axisKey,
  axisTotal,
  breakHolding,
  cellOffset,
  elementAt,
  hiddenLength,
  monotonic,
  nearestOffsetOf,
  offsetOf,
  pieceAtOffset,
  pieceOfPosition,
  pieceOffsets,
  piecesTotal,
  positionAt,
  wholePieces,
} from './axis.js';
export type { ResolvedBridge, StackEntry } from './bridges.js';
export { type BlockGroup, blockOverridesOf, blockUnits, opensComplexes } from './blocks.js';
export {
  type BlockPlace,
  blockElementAt,
  blockPositionAt,
  blockSpan,
  blockUnitOf,
  placeOfUnit,
  placesInto,
  placesOf,
  sameBlocks,
  type SystemPosition,
  unitOfPlace,
} from './places.js';
export { ViewIssues } from './issues.js';
export { type LiftData, regionMaps } from './lift.js';
export { type MemberAxis, type MemberOrder, memberAxis, memberOrderOf } from './members.js';
export {
  type Axis,
  type Follow,
  type FollowWindow,
  followRange,
  followSpan,
  followWindow,
  type LayerData,
  type DataForms,
  splitDataForms,
  type LayerType,
  LayerTypes,
  type LayerKinds,
  type LayerShape,
  BuiltInLayerKinds,
  type ResolvedLayer,
  type GridAxis,
  type ResolvedBlock,
  type ResolvedSection,
  type ResolvedTrack,
  type ResolvedNavigator,
  type ResolvedView,
  BridgePass,
  resolveView,
  resolveViewNode,
  resolveNavigatorNode,
  resolveViews,
  type TrackTemplate,
  type ViewLink,
} from './view.js';
export { expandTemplate, templateNeeds, templateNodes } from './templates.js';
export { type ProjectOptions, projectOntoView } from './project.js';
