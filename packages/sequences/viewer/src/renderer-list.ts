/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Every kernel type's renderer and `seqv:map`, one module each under `renderers/` (seqspec:viewer#extension-points), in
// the format's type order; `DefaultRegistry` registers them. A host that replaces one registers a copy of
// `DefaultRegistry`.
import type { ColumnRenderer, LayerRenderer } from "./registry.ts";
import { ColorRenderer } from "./renderers/color.ts";
import { CategoryColumn, FractionColumn, NumberColumn, TextColumn } from "./renderers/columns.ts";
import { DotPlotRenderer } from "./renderers/dot-plot.ts";
import { EmphasisRenderer } from "./renderers/emphasis.ts";
import { LettersRenderer } from "./renderers/letters.ts";
import { LogoRenderer } from "./renderers/logo.ts";
import { MapRenderer } from "./renderers/map.ts";
import { MatrixRenderer } from "./renderers/matrix.ts";
import { PairsRenderer } from "./renderers/pairs.ts";
import { PointsRenderer } from "./renderers/points.ts";
import { RangesRenderer } from "./renderers/ranges.ts";
import { RowsRenderer } from "./renderers/rows.ts";
import { RulerRenderer } from "./renderers/ruler.ts";
import { SecondaryRenderer } from "./renderers/secondary.ts";
import { ValuesRenderer } from "./renderers/values.ts";
import { VariantsRenderer } from "./renderers/variants.ts";

export const KernelRenderers: readonly LayerRenderer[] = [
  LettersRenderer,
  RangesRenderer,
  PointsRenderer,
  EmphasisRenderer,
  VariantsRenderer,
  SecondaryRenderer,
  ValuesRenderer,
  ColorRenderer,
  RulerRenderer,
  PairsRenderer,
  RowsRenderer,
  LogoRenderer,
  MatrixRenderer,
  DotPlotRenderer,
  MapRenderer,
] as readonly LayerRenderer[];

/** The built-in data columns; `DefaultRegistry` registers them. */
export const KernelColumns: readonly ColumnRenderer[] = [NumberColumn, FractionColumn, TextColumn, CategoryColumn];
