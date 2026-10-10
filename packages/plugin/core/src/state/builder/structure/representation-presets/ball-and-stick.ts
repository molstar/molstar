/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { StateObjectRef } from '@molstar/core/state';
import { BallAndStickRepresentationProvider } from '@molstar/graphics/repr/structure/representation/ball-and-stick';
import { BallAndStick } from '@molstar/plugin/registry/structure/ball-and-stick';
import { mergeRegistryEntries } from '@molstar/plugin/registry/merge';
import type { PluginRegistryEntry } from '@molstar/plugin/spec';
import { StructureRepresentationPresetProvider, presetStaticComponent, BuiltInPresetGroupName } from './types.js';

const CommonParams = StructureRepresentationPresetProvider.CommonParams;
const reprBuilder = StructureRepresentationPresetProvider.reprBuilder;
const updateFocusRepr = StructureRepresentationPresetProvider.updateFocusRepr;

/**
 * Shows the whole structure as ball-and-stick, regardless of its size. It is not part of the default presets; apps
 * list `BallAndStickPresetEntry`, which brings the representation and the themes it builds with.
 */
export const BallAndStickPreset = StructureRepresentationPresetProvider({
  id: 'preset-structure-representation-ball-and-stick',
  alias: 'ball-and-stick',
  display: {
    name: 'Ball & Stick',
    group: BuiltInPresetGroupName,
    description: 'Shows everything as Ball & Stick.',
  },
  params: () => CommonParams,
  async apply(ref, params, plugin) {
    const structureCell = StateObjectRef.resolveAndCheck(plugin.state.data, ref);
    if (!structureCell) return {};

    const components = {
      all: await presetStaticComponent(plugin, structureCell, 'all'),
    };

    const structure = structureCell.obj!.data;

    const { update, builder, typeParams, color, ballAndStickColor } = reprBuilder(plugin, params, structure);

    const representations = {
      all: builder.buildRepresentation(
        update,
        components.all,
        { type: BallAndStickRepresentationProvider, typeParams, color, colorParams: ballAndStickColor },
        { tag: 'all' },
      ),
    };

    await update.commit({ revertOnError: true });
    await updateFocusRepr(
      plugin,
      structure,
      params.theme?.focus?.name ?? color,
      params.theme?.focus?.params ?? ballAndStickColor,
    );

    return { components, representations };
  },
});

/** The ball-and-stick preset with the representation and themes it builds. */
export const BallAndStickPresetEntry: PluginRegistryEntry = mergeRegistryEntries(
  { structure: { presets: { representation: [BallAndStickPreset] } } },
  BallAndStick,
);
