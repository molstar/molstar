/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// @molstar/seqx-core/render: tier 3 (_spec/README.md "Render"): the layer-type contract and the built-in types, the
// section layouts (strips, warps, placements, navigation) and the tracks layout, the window mark, draw lists in passes, hit-testing, the token table, the text metrics port, the recording and SVG
// backends. Named re-exports only; what a minimal composition does not ask for is on its own entry:
// `render/render-svg` (`renderSvg`, the built-ins), `render/export-png`, `render/built-ins`, `render/navigator`,
// `render/navigator-picture` (the built-in pictures), `render/navigator-export`.

export { AxisTransform, type BandViewport, Cut, type Viewport } from './axis.js';
export {
  DefaultRenderConfig,
  type DensityStep,
  DensitySteps,
  type DensityThresholds,
  type DensityUnits,
  densityStep,
  fontSize,
  perPosition,
  type RenderConfig,
  type RenderConfigOverrides,
  renderConfig,
} from './config.js';
export {
  Anchor,
  ClipList,
  closePath,
  cubicTo,
  DrawLists,
  EmphasisPasses,
  GlyphList,
  type ImageData,
  ImageList,
  lineTo,
  moveTo,
  NoColor,
  Op,
  type OpKind,
  PassList,
  type PassName,
  Passes,
  type PassVisitor,
  PassWriter,
  PathCommand,
  PathCommandPoints,
  PathList,
  RectList,
  Round,
  TextList,
  walkPass,
} from './draw-lists.js';
export {
  type AttachOptions,
  attachSources,
  NoAssets,
  type OpenedModel,
  type OpenModelOptions,
  openModel,
  type PositionsViewport,
  type RenderSvgOptions,
  type RenderSvgResult,
  exportSvg,
  viewSources,
} from './export.js';
export { RenderIssues } from './issues.js';
export {
  DefaultThemeContext,
  Extents,
  LayerChannels,
  PlainDatum,
  RowDatum,
  type ChannelContext,
  type ThemeContext,
} from './channels.js';
export { type LayerColors, layerColors, type ThemeInput, viewTheme } from './layer-color.js';
export { LegendBand, layoutLegends, legendsHeight, viewLegends } from './legends.js';
export { OverlaySources } from './overlay-sources.js';
export { LaneIndex, StickyLanes } from './lanes.js';
export { type LayerSource, layerSource, readTiles, type TileContext, type TileRead } from './layer-data.js';
export {
  type Box,
  type Emphasis,
  type EmphasisHighlight,
  type HighlightStyle,
  HighlightStyles,
  EmptyOverlay,
  type Hit,
  type HitMark,
  type HitTarget,
  type LabelSink,
  type ItemDescription,
  type LayerType,
  type ContentContext,
  type LayoutContext,
  type Overlay,
  type StoreReader,
  sectionId,
} from './layer-type.js';
export { BlockMetrics, BlockPlacement } from './layout/blocks.js';
export {
  openGroupOf,
  RingMetrics,
  type RingParams,
  RingPlacement,
  RingsLayout,
  RingsNav,
  ringParams,
} from './layout/rings.js';
export { PlacePoint, StripPlacement, Strips, WarpKind } from './layout/strips.js';
export {
  lineEnd,
  lineOf,
  lineStart,
  lineWidthOf,
  type TextGeometry,
  TextLayout,
  TextNav,
  type TextParams,
  TextPlacement,
  textGeometry,
  textMinScale,
  textParams,
} from './layout/text.js';
export { clampWindow, sameWindow, TracksLayout, TracksNav } from './layout/tracks.js';
export { ContentSizing, type LayerContent, type SizedTrack, type TrackSizing, type TrackSizingHook } from './sizing.js';
export { GridLayout, GridNav, GridPlacement } from './layout/grid.js';
export type {
  NavContext,
  PlaceContext,
  Placement,
  PointKind,
  SectionLayoutType,
  SectionNav,
  StepDir,
  WindowPart,
} from './layout/types.js';
export { PolarWarp, TranslateWarp } from './layout/warp.js';
export {
  DefaultZoomStrip,
  drawNavigatorWindow,
  drawWindowArgs,
  floorExtent,
  navigatorPart,
  NavWin,
  NavWindow,
  NavWindowFlags,
  WindowFill,
  windowPart,
} from './window.js';
export type { NavigatorPicture, NavigatorPictureContext, NavigatorPictureTracks } from './navigator-picture.js';
export { type NarrowColumns, narrowColumns, overviewWidth } from './narrow.js';
export { type SourceLabelsOptions, withSourceNumberings } from './numberings.js';
export {
  drawRail,
  hitRail,
  MaxRowScale,
  MinRowScale,
  placeRail,
  placeRowsRail,
  Rail,
  RowsRail,
  type PositionsRailMode,
  type RailPlace,
  type RailRow,
  railRowList,
  railRows,
  railSections,
  railUnitAt,
  rowsRailTracks,
  type RowsRailMode,
  SlimRailPx,
} from './rail.js';
export { type BandsLay, type BandsPrep, BandsType } from './bands.js';
export { type LettersLay, type LettersPrep, LettersType, shownLetters } from './letters.js';
export { type RangesLay, type RangesPrep, RangesType } from './ranges.js';
export { fmt, recordDrawLists } from './record.js';
export { labelStep, type RulerLay, type RulerPrep, RulerType } from './ruler.js';
export { type Frame, type FrameInput, Renderer, type RendererOptions } from './renderer.js';
export { type SecondaryLay, type SecondaryPrep, SecondaryType, secondaryKind } from './secondary.js';
export { SquircleHandle, squirclePath, squirclePathArgs, squircleRadius } from './squircle.js';
export { drawListsToSvg, escapeXml, type SvgAccessible, type SvgOptions } from './svg.js';
export {
  accessibleOf,
  itemText,
  locationText,
  MaxNamedItems,
  MaxNamedItemsPerTrack,
  rangesText,
  spanText,
} from './accessibility.js';
export { bundledTextMetrics, bundledWidth, createTextMetrics, parseFont, type TextMetrics } from './text.js';
export {
  builtInTokens,
  CategoryCount,
  colorHex,
  createTokens,
  DarkTokens,
  DefaultFonts,
  type FontStacks,
  inkOn,
  LightTokens,
  mixColors,
  parseColor,
  resolveColor,
  type SeqTokenName,
  SeqTokenNames,
  type SeqTokens,
  withAlpha,
} from './tokens.js';
export {
  type BridgeLayout,
  bridgeAt,
  drawView,
  type HeldDomains,
  hitView,
  type LayerLayout,
  type LayoutInput,
  LayoutState,
  layerKey,
  layoutView,
  multiSystem,
  type SectionLayout,
  sectionLayoutName,
  type TrackLayout,
  type ViewLayout,
} from './tracks.js';
export {
  actionColumns,
  builtInHeaderStates,
  collapsedSection,
  type ActionColumns,
  type HeaderInputFacts,
  type HeaderIcon,
  type HeaderLine,
  HeaderMenuName,
  HeaderMenuPart,
  type HeaderParts,
  type HeaderPartsOf,
  type HeaderSlotBox,
  type HeaderSlotPart,
  type HeaderStateBox,
  HeaderStateNames,
  type HeaderStatePart,
  type HiddenTrack,
  hiddenTrackName,
  type SectionHeaderLine,
} from './headers.js';
export { type ValuesLay, type ValuesPrep, ValuesType } from './values.js';
export { type VariantsLay, type VariantsPrep, VariantsType } from './variants.js';
export { type FrameLay, type FramePrep, FrameType } from './frames.js';
export { LogoLay, type LogoPrep, LogoType, MinLogoGlyphPx } from './logo.js';
export { GridLay, type GridPrep, type GridSummary, GridType, gridSummary } from './grid.js';
export { MatrixLay, type MatrixPrep, MatrixType } from './matrix.js';
export {
  binsText,
  genotypeText,
  gridCellHover,
  matrixCellHover,
  pairsBinHover,
  type RegionHover,
  type SystemName,
} from './region-hover.js';
export { crosses, PairsLay, type PairsPrep, PairsType, SplitExact, splitCrossing } from './pairs.js';
export { readHeld, RowsLay, type RowsOrder, type RowsPrep, RowsType, rowsSource } from './rows.js';
export {
  allClades,
  buildClades,
  type Clades,
  cladeNamed,
  cladesAtDepth,
  cladesHolding,
  cladesWithin,
  type ColumnMeasures,
  columnMeasures,
  consensus,
  DisplayRows,
  type HeldCells,
  logoHeights,
  measureEntries,
  type Profile,
  profileOf,
} from './rows-model.js';
