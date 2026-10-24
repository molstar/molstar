/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { StateObjectRef } from '@molstar/core/state';
import { LineRepresentationProvider } from '@molstar/graphics/repr/structure/representation/line';
import { mergeRegistryEntries } from '@molstar/plugin/registry/merge';
import { Line } from '@molstar/plugin/registry/structure/line';
import type { PluginRegistryEntry } from '@molstar/plugin/spec';
import {
  presetStaticComponent,
  StructureRepresentationPresetProvider,
} from '@molstar/plugin/state/builder/structure/representation-presets/types';
import { colorTheme } from '../prism-viewer/reps.js';

// VMD's default look as a representation preset (_spec §5.3): one rep, Lines of everything, colored by Name (VMD's
// element colors, carbon cyan). VMD's dialect loads with it, so `mol delrep 0 top` then `mol addrep top` scripts end as
// in VMD; registered with the plugin, it is also in Prism's Preset menu for anyone.

export const VmdLinesPresetId = 'command-language-vmd-lines';

export const VmdLinesPreset = StructureRepresentationPresetProvider({
  id: VmdLinesPresetId,
  display: {
    name: 'VMD lines',
    group: 'Other',
    description: 'VMD’s default: lines, colored by element with VMD’s colors.',
  },
  params: () => StructureRepresentationPresetProvider.CommonParams,
  async apply(ref, _params, plugin) {
    const structureCell = StateObjectRef.resolveAndCheck(plugin.state.data, ref);
    if (!structureCell) return {};
    const components = { all: await presetStaticComponent(plugin, structureCell, 'all') };
    const update = plugin.state.data.build();
    const theme = colorTheme({ theme: 'element-symbol', params: { palette: 'vmd' } });
    const representations = {
      all: plugin.builders.structure.representation.buildRepresentation(
        update,
        components.all,
        { type: LineRepresentationProvider, color: theme.name as 'element-symbol', colorParams: theme.params },
        { tag: 'all' },
      ),
    };
    await update.commit({ revertOnError: true });
    return { components, representations };
  },
});

/** The preset with the representation and theme it builds. */
export const VmdLinesPresetEntry: PluginRegistryEntry = mergeRegistryEntries(
  { structure: { presets: { representation: [VmdLinesPreset] } } },
  Line,
);
