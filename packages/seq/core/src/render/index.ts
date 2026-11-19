/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// @molstar/seqx-core/render: tier 3 (_spec/README.md "Render"): the layer-type contract and the built-in types, the
// section layouts (strips, warps, placements, navigation) and the tracks layout, the window mark, draw lists in passes, hit-testing, the token table, the text metrics port, the recording and SVG
// backends, and `renderSvg`. Named re-exports only.

export { AxisTransform, Cut, type Viewport } from './axis.js';
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
  type PositionsViewport,
  type RenderSvgOptions,
  type RenderSvgResult,
  renderSvg,
  viewSources,
} from './export.js';
export { RenderIssues } from './issues.js';
export {
  DefaultThemeContext,
  Extents,
  LayerChannels,
  PlainDatum,
  RowDatum,
  type ThemeContext,
} from './channels.js';
export { LegendBand, layoutLegends, legendsHeight, viewLegends } from './legends.js';
export { OverlaySources } from './overlay-sources.js';
export { LaneIndex, StickyLanes } from './lanes.js';
export { type LayerSource, layerSource } from './layer-data.js';
export {
  type Box,
  type Emphasis,
  type EmphasisHighlight,
  EmptyOverlay,
  type Hit,
  type LabelSink,
  type LayerType,
  type LayoutContext,
  type Overlay,
  type StoreReader,
  sectionId,
} from './layer-type.js';
export {
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
export { drawWindow, WindowMetrics, windowExtent, windowPart } from './window.js';
export { type LettersLay, type LettersPrep, LettersType, shownLetters } from './letters.js';
export { type RangesLay, type RangesPrep, RangesType } from './ranges.js';
export { fmt, recordDrawLists } from './record.js';
export { labelStep, type RulerLay, type RulerPrep, RulerType } from './ruler.js';
export { type Frame, type FrameInput, Renderer, type RendererOptions } from './renderer.js';
export { type SecondaryLay, type SecondaryPrep, SecondaryType, secondaryKind } from './secondary.js';
export { SquircleHandle, squirclePath, squircleRadius } from './squircle.js';
export { drawListsToSvg, escapeXml, type SvgOptions } from './svg.js';
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
  BuiltInLayerTypes,
  BuiltInLayouts,
  drawView,
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
export { type ValuesLay, type ValuesPrep, ValuesType } from './values.js';
export { type VariantsLay, type VariantsPrep, VariantsType } from './variants.js';
export { type FrameLay, type FramePrep, FrameType } from './frames.js';
