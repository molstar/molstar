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
import type { ColumnRenderer, LayerRenderer } from "./registry.js";
import { ColorRenderer } from "./renderers/color.js";
import { CategoryColumn, FractionColumn, NumberColumn, TextColumn } from "./renderers/columns.js";
import { DotPlotRenderer } from "./renderers/dot-plot.js";
import { EmphasisRenderer } from "./renderers/emphasis.js";
import { LettersRenderer } from "./renderers/letters.js";
import { LogoRenderer } from "./renderers/logo.js";
import { MapRenderer } from "./renderers/map.js";
import { MatrixRenderer } from "./renderers/matrix.js";
import { PairsRenderer } from "./renderers/pairs.js";
import { PointsRenderer } from "./renderers/points.js";
import { RangesRenderer } from "./renderers/ranges.js";
import { RowsRenderer } from "./renderers/rows.js";
import { RulerRenderer } from "./renderers/ruler.js";
import { SecondaryRenderer } from "./renderers/secondary.js";
import { ValuesRenderer } from "./renderers/values.js";
import { VariantsRenderer } from "./renderers/variants.js";

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
