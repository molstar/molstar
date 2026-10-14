/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The renderer API (seqspec:viewer#extension-points, seqspec:viewer#the-render-backend): everything a renderer,
// built-in or a host's, may use, versioned with the package and evolving with the built-ins, never frozen apart.
// Renderer modules import the viewer through this module (or `/draw`) and `/testing` only; the composer, surface,
// layout manager, MarkSet internals, region cache and graph are not part of it (Biome: lib/seqstar/layers.grit).
//   - drawing: the batches and the backend's draw methods (not its pass, clip and copy methods, the composer's),
//     `DrawView` and the placement helpers;
//   - headless: the axis view and its converters, hits, themes and legends, the prepares and their helpers (segment
//     tables, pyramids, density), text, marks as the overlay reads them, the stream view, the renderer entry types.
// Every name is listed (VL-16): what the built-in renderers, the interactions renderer, the conformance kit and a
// backend use, and every type their signatures name (a navigator's `NavView`, a sink, a pyramid: LC-02). Navigator,
// ring-label, matrix and dot-plot internals are not part of it (tests reach them by path).
import type { MarkSet } from "../marks.ts";

export { type Axis, AxisGeometry, type AxisView, cellGap, cellSpan, type ShownIndex } from "../axis.ts";
export type { Skin } from "../graph.ts";
export { copyAt, Hit, HoverKind, sourceAt } from "../hit.ts";
export type { HoverState } from "../hover.ts";
export type { LabelContext, LabelSink } from "../labels.ts";
export { letterCap, type RenderProps } from "../layout.ts";
export { MarkBit, type MarkSet, type OverlayMarks } from "../marks.ts";
export { CellIO, cellLayout, endsIn, litPairs, pairsAt, pairsInView } from "../pairs-index.ts";
export {
  type EmphasisState,
  type Features,
  featureAt,
  keyOf,
  LetterForm,
  type LettersState,
  type OfLetters,
  ofPartAt,
  type PairsState,
  type PerPosition,
  type PreparedPart,
  PreparedParts,
  positionKey,
  type RulerSource,
  type RulerState,
  rowValues,
  SecondaryKind,
  type SplitRows,
  type StreamFields,
} from "../prepare.ts";
export { binValues, PYRAMID_BASE, segWindow, type ValuePyramid } from "../pyramids.ts";
export type { CellColoring } from "../reductions.ts";
export type {
  ColumnRenderer,
  LayerInput,
  LayerRenderer,
  NamedEntry,
  PrepareContext,
  SeqViewerRegistryEntry,
} from "../registry.ts";
export {
  Arcs,
  hitItem,
  kernelEntry,
  LETTER_LINE,
  laneHeight,
  laneInset,
  lanePad,
  legendOf,
  letterLines,
  wedgeLabel,
} from "../renderer-helpers.ts";
export { MatrixBands, type MatrixState, prepareMatrix, RowsHeader, type RowsState } from "../rows.ts";
export {
  GROUP_LABEL_PAD,
  groupLabelPx,
  type LaneSegments,
  type LineTable,
  laneLineHeights,
  lineLanes,
  type PreparedRun,
  runXInto,
} from "../runs.ts";
export {
  blockLevels,
  LETTERS_PER_PX,
  type LocalCells,
  overviewLevels,
  RegionState,
  type StreamItems,
  type StreamLayer,
  type StreamView,
  streamedBlock,
  streamStats,
} from "../streams.ts";
export { type FontId, fitCount, fitIn, fitted, type PreparedLabel, type TextManagerApi } from "../text.ts";
export {
  CellInk,
  cellInk,
  dataIndex,
  type FieldTable,
  type GlyphTable,
  Granularity,
  glyphTableOf,
  type LegendSink,
  type LegendSizeOf,
  NO_GLYPH,
  NO_SIZE,
  OPAQUE,
  SIZE_ONE,
  SizeUnit,
  sizeOf,
  type Theme,
  type ThemeContext,
  ThemeLocation,
  type ThemeProvider,
} from "../themes.ts";
export { rowAt, rowLabel, type TreeGeo, type TreeRows } from "../tree.ts";
export { Lod } from "../viewport.ts";
export { inkOf, pushArc } from "./arcs.ts";
export {
  BackendCaps,
  type BackendImages,
  CHARS,
  type DrawResources,
  GlyphBatch,
  type GlyphRecord,
  GlyphShape,
  NO_COLOR,
  Pass,
  PathBatch,
  PathMode,
  PathVerb,
  PathVerbCoords,
  RectBatch,
  type RenderBackend,
  SectorBatch,
  SHAPE_COUNT,
  SkinToken,
  TextAlign,
  TextBaseline,
  TextBatch,
} from "./backend.ts";
export {
  type CellSource,
  cellBand,
  columnBins,
  columnEnd,
  columnLevel,
  drawCells,
  shownHull,
  visibleRange,
} from "./cells.ts";
export type { ComposerFrame } from "./composer-frame.ts";
export {
  DensityKind,
  densityAdd,
  densityColumns,
  densityDraw,
  densityMode,
  densityShown,
  densitySlot,
  densityStart,
  drawDensity,
  featureGlyph,
  featureKey,
  itemsInView,
  keyAlpha,
  packAlpha,
} from "./density.ts";
export {
  type DotPlotGrid,
  type DotPlotLevel,
  type DotPlotState,
  type DotPlotTileSet,
  DotPlotTiles,
  dotPlotAt,
  dotPlotValue,
  drawDotPlot,
  gridAt,
  overlayDotPlot,
  prepareDotPlot,
  type TileLevel,
} from "./dot-plot.ts";
export {
  densityBand,
  drawLineLanes,
  featureHoverKey,
  featureInk,
  hitLineLanes,
  hoveredAt,
  laneRows,
  outlineLineLanes,
} from "./features.ts";
export {
  HOVER_WIDTH,
  HoverOutline,
  hoverEmphasis,
  hoveredItem,
  hoverFill,
  hoverInk,
  hoverWidth,
  type ItemOrder,
  outlineItem,
  ringItem,
} from "./hover.ts";
export { elementSymbol, letterFontPx, letterRoom, textShown } from "./letter-text.ts";
export {
  drawMatrixTiles,
  type MatrixTileSet,
  MatrixTiles,
  matrixBandRows,
  matrixFrame,
  matrixHitBand,
  outlineMatrixCell,
} from "./matrix-tiles.ts";
export type { Miniature } from "./miniature.ts";
export { KeyedImage, type NavView, WINDOW_WASH, windowBox, windowOverMap } from "./nav.ts";
export { inlineGroupLabels, numberText } from "./numbers.ts";
export { band, boundaryEdge, mark, placeRun } from "./place.ts";
export { drawPoints, hitPoints, outlinePoints } from "./points.ts";
export type { RingLabelSink } from "./ring-labels.ts";
export {
  bar,
  barBand,
  DOT_R,
  featureSize,
  hitLayer,
  hitWindow,
  type SizeView,
  sizePx,
  sizeTheme,
} from "./sizes.ts";
export {
  bandOver,
  drawStreamedItems,
  hitStreamed,
  lettersGate,
  outlineStreamed,
  owned,
  type RegionDraw,
  streamColumn,
} from "./streamed.ts";
export { drawTree, drawTreeHover } from "./tree-paths.ts";
export { DrawView, type RingView, type TracksHost } from "./view.ts";

/** A channel's lit runs as a renderer reads them (`DrawView.marks.runs(bit)`): per block, sorted disjoint runs. */
export type MarkRuns = Pick<MarkSet, "start" | "stop" | "count" | "from" | "to" | "has" | "any" | "first">;
