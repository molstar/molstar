/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { InteractionsRepresentationProvider } from '@molstar/graphics/props/computed/representations/interactions';
import { InteractionTypeColorThemeProvider } from '@molstar/graphics/props/computed/themes/interaction-type';
import { ElementSymbolColorThemeProvider } from '@molstar/graphics/theme/color/element-symbol';
import { UniformColorThemeProvider } from '@molstar/graphics/theme/color/uniform';
import { PhysicalSizeThemeProvider } from '@molstar/graphics/theme/size/physical';
import { UncertaintySizeThemeProvider } from '@molstar/graphics/theme/size/uncertainty';
import { UniformSizeThemeProvider } from '@molstar/graphics/theme/size/uniform';
import { mergeRegistryEntries } from '@molstar/plugin/registry/merge';
import { Backbone } from '@molstar/plugin/registry/structure/backbone';
import { BallAndStick } from '@molstar/plugin/registry/structure/ball-and-stick';
import { Carbohydrate } from '@molstar/plugin/registry/structure/carbohydrate';
import { Cartoon } from '@molstar/plugin/registry/structure/cartoon';
import { GaussianSurface } from '@molstar/plugin/registry/structure/gaussian-surface';
import { Line } from '@molstar/plugin/registry/structure/line';
import { MolecularSurface } from '@molstar/plugin/registry/structure/molecular-surface';
import { Putty } from '@molstar/plugin/registry/structure/putty';
import { Spacefill } from '@molstar/plugin/registry/structure/spacefill';
import { Isosurface } from '@molstar/plugin/registry/volume/isosurface';
import { Slice } from '@molstar/plugin/registry/volume/slice';
import type { PluginRegistryEntry } from '@molstar/plugin/spec';

/**
 * A catalog module: it lists every representation and color/size theme that MolViewSpec documents and the MVS loading
 * extensions can name, so that an application can include exactly what MVS needs. MVS documents select representations
 * and themes by name, which makes MVS a deliberate consumer of catalogs.
 *
 * - Structure representations: the nine of the MVS `representation` node (cartoon, backbone, ball-and-stick, line,
 *   spacefill, carbohydrate, gaussian-surface, molecular-surface, putty), each through its plugin entry, which carries
 *   the representation and the providers of its default themes; plus `interactions`, which the
 *   non-covalent-interactions loading extension names.
 * - Volume representations: the two of the MVS `volume_representation` node (isosurface and slice), through their
 *   plugin entries.
 * - Structure color themes: `uniform` (MVS colors), `element-symbol` and `interaction-type` (the
 *   non-covalent-interactions extension), plus the default color themes brought in by the representation entries.
 * - Structure size themes: `uniform`, `physical`, and `uncertainty`, plus the default size themes brought in by the
 *   representation entries.
 * - Volume color themes: `uniform`, plus the default themes brought in by the representation entries.
 *
 * What it does not list, because the `MolViewSpec` behavior registers it with `plugin.register`: the MVS label
 * representations (`custom-label`, `mvs-annotation-label`) and the MVS color themes (`mvs-split-uniform`,
 * `mvs-annotation`, and `mvs-multilayer`, which is built from the plugin's color theme registry). A document that names
 * a theme with the `molstar_color_theme_name` custom property, or sets representation parameters with
 * `molstar_representation_params`, can name any provider; the app registers those itself.
 *
 * An MVS-only app lists this entry, the format entries and markdown extensions MVS uses, and the `MolViewSpec`
 * behavior.
 */
export const MVSRuntimeRegistry: PluginRegistryEntry = mergeRegistryEntries(
  Cartoon,
  Backbone,
  BallAndStick,
  Line,
  Spacefill,
  Carbohydrate,
  GaussianSurface,
  MolecularSurface,
  Putty,
  {
    structure: {
      representations: [InteractionsRepresentationProvider],
      themes: {
        color: [UniformColorThemeProvider, ElementSymbolColorThemeProvider, InteractionTypeColorThemeProvider],
        size: [UniformSizeThemeProvider, PhysicalSizeThemeProvider, UncertaintySizeThemeProvider],
      },
    },
  },
  Isosurface,
  Slice,
  { volume: { themes: { color: [UniformColorThemeProvider] } } },
);
