/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// @molstar/seqx-core/theme (_spec/README.md "Themes"): the theme sublanguage's evaluation without the layout, so another
// reader evaluates a theme as the canvas and the SVG do: forms compiled to channels, palettes over Mol*'s value
// palettes, the scheme registry and the built-in schemes, measures, legends' data and their converters. Named
// re-exports only.

export {
  Black,
  colorHex,
  contrastRgba,
  fromRgb,
  inkFor,
  isColorSpelling,
  parseColor,
  resolveColor,
  type TokenTable,
  toRgb,
  White,
} from './colors.js';
export {
  type Channel,
  compileChannel,
  GreyChroma,
  glyphAgainst,
  inkAgainst,
  type LegendSource,
  type MeasureEntry,
  type ThemeDatum,
  type ThemeEnv,
  type ThemeIssueCode,
} from './compile.js';
export {
  BuiltInMeasures,
  type ChannelKind,
  type ChannelName,
  ChannelWords,
  ClauseKinds,
  type ClauseKind,
  isThemeNode,
  PaletteKinds,
  type PaletteKind,
  param,
  type ScaleParams,
  type ThemeBody,
  ThemeChannels,
  ThemeForms,
  type ThemeFormKind,
  type ThemeNode,
} from './forms.js';
export {
  formatValue,
  type Legend,
  type LegendEntries,
  type LegendNote,
  type LegendOptions,
  legendOf,
  legendsEqual,
  legendTitle,
  qualifyLegends,
  mergeLegends,
  niceTicks,
  toMolLegend,
} from './legend.js';
export {
  compilePalette,
  continuousLut,
  defaultPalette,
  fromMvsPalette,
  fromMvsScale,
  LutSize,
  namedList,
  normalize,
  type Palette,
  scaled,
  TokenLists,
  toMvsPalette,
} from './palettes.js';
export {
  BuiltInSchemes,
  builtInSchemes,
  DarkBand,
  darkVariant,
  lettersDefault,
  type SchemeEntry,
  SchemeRegistry,
  type SchemeSubject,
  schemeColors,
} from './schemes.js';
