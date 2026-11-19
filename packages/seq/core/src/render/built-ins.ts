/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { BandsType } from './bands.js';
import { FrameType } from './frames.js';
import { GridType } from './grid.js';
import type { LayerType } from './layer-type.js';
import { GridLayout } from './layout/grid.js';
import { RingsLayout } from './layout/rings.js';
import { TextLayout } from './layout/text.js';
import { TracksLayout } from './layout/tracks.js';
import type { SectionLayoutType } from './layout/types.js';
import { LettersType } from './letters.js';
import { LogoType } from './logo.js';
import { PairsType } from './pairs.js';
import { RangesType } from './ranges.js';
import { RowsType } from './rows.js';
import { RulerType } from './ruler.js';
import { SecondaryType } from './secondary.js';
import { ValuesType } from './values.js';
import { VariantsType } from './variants.js';

// The core's layer types and section layouts (_spec/README.md "Context", Registry entries): the default lists a
// composition registers, never a table a renderer falls back to. A composition that lists fewer ships fewer: nothing
// imports this module but the default compositions (`createSeqContext`, `renderSvg`, the viewer's default spec).

/** The layer types the core ships, by name: what the default compositions register. */
export const BuiltInLayerTypes: Readonly<Record<string, LayerType<any, any>>> = {
  letters: LettersType,
  ranges: RangesType,
  values: ValuesType,
  ruler: RulerType,
  secondary: SecondaryType,
  variants: VariantsType,
  frame: FrameType,
  rows: RowsType,
  logo: LogoType,
  bands: BandsType,
  pairs: PairsType,
  grid: GridType,
};

/** The section layouts the core ships, by name (_spec/README.md "Render", Section layouts). */
export const BuiltInLayouts: Readonly<Record<string, SectionLayoutType<any>>> = {
  tracks: TracksLayout,
  text: TextLayout,
  rings: RingsLayout,
  grid: GridLayout,
};

/** What a renderer draws with, the built-ins: spread into a `Renderer`'s or a layout's options. */
export const BuiltInRender = { types: BuiltInLayerTypes, layouts: BuiltInLayouts } as const;
