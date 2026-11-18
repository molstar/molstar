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
  positionAt,
  wholePieces,
} from './axis.js';
export { ViewIssues } from './issues.js';
export {
  type Axis,
  type Follow,
  followRange,
  type LayerData,
  type LayerType,
  LayerTypes,
  type ResolvedLayer,
  type ResolvedSection,
  type ResolvedTrack,
  type ResolvedView,
  resolveView,
  resolveViews,
  type TrackTemplate,
  type ViewLink,
} from './view.js';
export { expandTemplate, templateNeeds, templateNodes } from './templates.js';
export { type ProjectOptions, projectOntoView } from './project.js';
