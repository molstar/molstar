/**
 * Copyright (c) 2019-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Adam Midlik <midlik@gmail.com>
 * @author Gianluca Tomasello <giagitom@gmail.com>
 */

import { StateObjectRef } from '@molstar/core/state';
import { mergeRegistryEntries } from '@molstar/plugin/registry/merge';
import type { PluginRegistryEntry } from '@molstar/plugin/spec';
import { SpacefillRepresentationProvider } from '@molstar/graphics/repr/structure/representation/spacefill';
import { IllustrativeColorThemeProvider } from '@molstar/graphics/theme/color/illustrative';
import { Spacefill } from '@molstar/plugin/registry/structure/spacefill';
import { StructureRepresentationPresetProvider, presetStaticComponent } from './types.js';

const CommonParams = StructureRepresentationPresetProvider.CommonParams;
const reprBuilder = StructureRepresentationPresetProvider.reprBuilder;
const updateFocusRepr = StructureRepresentationPresetProvider.updateFocusRepr;

export const IllustrativePreset = StructureRepresentationPresetProvider({
  id: 'preset-structure-representation-illustrative',
  alias: 'illustrative',
  display: {
    name: 'Illustrative',
    group: 'Miscellaneous',
    description: 'Show everything in spacefill representation with illustrative colors and ignore light.',
  },
  params: () => CommonParams,
  async apply(ref, params, plugin) {
    const structureCell = StateObjectRef.resolveAndCheck(plugin.state.data, ref);
    if (!structureCell) return {};

    const components = {
      all: await presetStaticComponent(plugin, structureCell, 'all'),
      branched: undefined,
    };

    const structure = structureCell.obj!.data;

    const { update, builder, typeParams, color } = reprBuilder(plugin, params, structure);

    const representations = {
      all: builder.buildRepresentation(
        update,
        components.all,
        {
          type: SpacefillRepresentationProvider,
          typeParams: { ...typeParams, ignoreLight: true },
          color: IllustrativeColorThemeProvider,
          // the params are merged over the defaults of the theme, so only the changed ones are given
          colorParams: { style: { name: 'entity-id', params: { overrideWater: true } } as any },
        },
        { tag: 'all' },
      ),
    };
    await update.commit({ revertOnError: true });
    await updateFocusRepr(plugin, structure, params.theme?.focus?.name ?? color, params.theme?.focus?.params);

    return { components, representations };
  },
});

/** The IllustrativePreset preset with the representations and themes it builds. */
export const IllustrativePresetEntry: PluginRegistryEntry = mergeRegistryEntries(
  { structure: { presets: { representation: [IllustrativePreset] } } },
  Spacefill,
  { structure: { themes: { color: [IllustrativeColorThemeProvider] } } },
);
