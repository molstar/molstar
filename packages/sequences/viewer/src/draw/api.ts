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
// layout manager, MarkSet internals, region cache and graph are not part of it (Biome: packages/sequences/layers.grit).
//   - drawing: the batches and the backend's draw methods (not its pass, clip and copy methods, the composer's),
//     `DrawView` and the placement helpers;
//   - headless: the axis view and its converters, hits, themes and legends, the prepares and their helpers (segment
//     tables, pyramids, density), text, marks as the overlay reads them, the stream view, the renderer entry types.
// Every name is listed (VL-16): what the built-in renderers, the interactions renderer, the conformance kit and a
// backend use, and every type their signatures name (a navigator's `NavView`, a sink, a pyramid: LC-02). Navigator,
// ring-label, matrix and dot-plot internals are not part of it (tests reach them by path).
import type { MarkSet } from '../marks.js';

export { type Axis, AxisGeometry, type AxisView, cellGap, cellSpan, type ShownIndex } from '../axis.js';
export type { Skin } from '../graph.js';
export { copyAt, Hit, HoverKind, sourceAt } from '../hit.js';
export type { HoverState } from '../hover.js';
export type { LabelContext, LabelSink } from '../labels.js';
export { letterCap, type RenderProps } from '../layout.js';
export { MarkBit, type MarkSet, type OverlayMarks } from '../marks.js';
export { CellIO, cellLayout, endsIn, litPairs, pairsAt, pairsInView } from '../pairs-index.js';
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
} from '../prepare.js';
export { binValues, PYRAMID_BASE, segWindow, type ValuePyramid } from '../pyramids.js';
export type { CellColoring } from '../reductions.js';
export type {
  ColumnRenderer,
  LayerInput,
  LayerRenderer,
  NamedEntry,
  PrepareContext,
  SeqViewerRegistryEntry,
} from '../registry.js';
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
} from '../renderer-helpers.js';
export { MatrixBands, type MatrixState, prepareMatrix, RowsHeader, type RowsState } from '../rows.js';
export {
  CHIP_PAD,
  chipFontPx,
  GROUP_LABEL_PAD,
  groupLabelPx,
  type LaneSegments,
  type LineTable,
  laneLineHeights,
  lineLanes,
  type PreparedRun,
  runXInto,
  type StripSource,
  stripLineHeights,
} from '../runs.js';
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
} from '../streams.js';
export { type FontId, fitCount, fitIn, fitted, type PreparedLabel, type TextManagerApi } from '../text.js';
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
  LineStyle,
  lineStyleId,
  NO_GLYPH,
  NO_LINE_STYLE,
  NO_SIZE,
  OPAQUE,
  SIZE_ONE,
  SizeUnit,
  sizeOf,
  type Theme,
  type ThemeContext,
  ThemeLocation,
  type ThemeProvider,
} from '../themes.js';
export { rowAt, rowLabel, type TreeGeo, type TreeRows } from '../tree.js';
export { Lod } from '../viewport.js';
export { inkOf, pushArc } from './arcs.js';
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
} from './backend.js';
export {
  type CellSource,
  cellBand,
  columnBins,
  columnEnd,
  columnLevel,
  drawCells,
  shownHull,
  visibleRange,
} from './cells.js';
export type { ComposerFrame } from './composer-frame.js';
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
} from './density.js';
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
} from './dot-plot.js';
export {
  densityBand,
  drawLineLanes,
  drawLineText,
  featureHoverKey,
  featureInk,
  hitLineLanes,
  hoveredAt,
  laneRows,
  outlineLineLanes,
} from './features.js';
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
} from './hover.js';
export { elementSymbol, letterFontPx, letterRoom, textShown } from './letter-text.js';
export { featureLineStyle, StyledLines, styledSegment } from './line-styles.js';
export {
  drawMatrixTiles,
  type MatrixTileSet,
  MatrixTiles,
  matrixBandRows,
  matrixFrame,
  matrixHitBand,
  outlineMatrixCell,
} from './matrix-tiles.js';
export type { Miniature } from './miniature.js';
export { KeyedImage, type NavView, WINDOW_WASH, windowBox, windowOverMap } from './nav.js';
export { inlineGroupLabels, numberText } from './numbers.js';
export { band, boundaryEdge, mark, placeRun } from './place.js';
export { drawPoints, hitPoints, outlinePoints } from './points.js';
export type { RingLabelSink } from './ring-labels.js';
export { drawStrip, outlineStrip, type StripColors, stripBand } from './strips.js';
export {
  bar,
  barBand,
  DOT_R,
  featureSize,
  hitLayer,
  hitReach,
  hitWindow,
  hitWindowAt,
  type SizeView,
  sizePx,
  sizeTheme,
} from './sizes.js';
export {
  bandOver,
  drawStreamedItems,
  hitStreamed,
  lettersGate,
  outlineStreamed,
  owned,
  type RegionDraw,
  streamColumn,
} from './streamed.js';
export { drawTree, drawTreeHover } from './tree-paths.js';
export { DrawView, type RingView, type TracksHost } from './view.js';

/** A channel's lit runs as a renderer reads them (`DrawView.marks.runs(bit)`): per block, sorted disjoint runs. */
export type MarkRuns = Pick<MarkSet, 'start' | 'stop' | 'count' | 'from' | 'to' | 'has' | 'any' | 'first'>;
