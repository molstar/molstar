/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Adam Midlik <midlik@gmail.com>
 * @author Ludovic Autin <ludovic.autin@gmail.com>
 */

import { PluginStateTransform, PluginStateObject as SO } from '@molstar/plugin/state/objects';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { Task } from '@molstar/core/task';
import { trajectoryFromCifCore } from '@molstar/model/formats/structure/cif-core';
import { trajectoryProps } from './helpers.js';
import { TrajectoryFormatProvider, defaultVisuals } from './provider.js';
import { TrajectoryFormatCategory } from './category.js';
import { guessCifVariant, applyTransformerRaw, rawDataObject } from '@molstar/plugin/state/formats/provider';
import { ParseCif } from '@molstar/plugin/state/formats/cif';

export { TrajectoryFromCifCore };
type TrajectoryFromCifCore = typeof TrajectoryFromCifCore;
const TrajectoryFromCifCore = PluginStateTransform.BuiltIn({
  name: 'trajectory-from-cif-core',
  display: {
    name: 'Parse CIF Core',
    description: 'Identify and create all separate models in the specified CIF data block',
  },
  from: SO.Format.Cif,
  to: SO.Molecule.Trajectory,
  params(a) {
    if (!a) {
      return {
        blockHeader: PD.Optional(
          PD.Text(void 0, {
            description: 'Header of the block to parse. If none is specifed, the 1st data block in the file is used.',
          }),
        ),
      };
    }
    const { blocks } = a.data;
    return {
      blockHeader: PD.Optional(
        PD.Select(
          blocks[0] && blocks[0].header,
          blocks.map((b) => [b.header, b.header] as [string, string]),
          { description: 'Header of the block to parse' },
        ),
      ),
    };
  },
})({
  apply({ a, params }) {
    return Task.create('Parse CIF Core', async (ctx) => {
      const header = params.blockHeader || a.data.blocks[0].header;
      const block = a.data.blocks.find((b) => b.header === header);
      if (!block) throw new Error(`Data block '${[header]}' not found.`);
      const models = await trajectoryFromCifCore(block).runInContext(ctx);
      if (models.frameCount === 0) throw new Error('No models found.');
      const props = trajectoryProps(models);
      return new SO.Molecule.Trajectory(models, props);
    });
  },
});

export const CifCoreProvider = TrajectoryFormatProvider({
  name: 'cifCore',
  label: 'cifCore',
  description: 'CIF Core',
  category: TrajectoryFormatCategory,
  stringExtensions: ['cif'],
  isApplicable: (info, data) => {
    if (info.ext === 'cif') return guessCifVariant(info, data) === 'coreCif';
    return false;
  },
  parse: async (plugin, data, params?) => {
    const state = plugin.state.data;
    const cif = state
      .build()
      .to(data)
      .apply(ParseCif, void 0, { state: { isGhost: true } });
    const trajectory = await cif
      .apply(TrajectoryFromCifCore, void 0, { tags: params?.trajectoryTags })
      .commit({ revertOnError: true });

    if ((cif.selector.cell?.obj?.data.blocks.length || 0) > 1) {
      plugin.state.data.updateCellState(cif.ref, { isGhost: false });
    }
    return { trajectory };
  },
  parseRaw: async (plugin, ctx, data) => {
    const cif = await applyTransformerRaw(plugin, ctx, ParseCif, rawDataObject(data));
    const trajectory = await applyTransformerRaw(plugin, ctx, TrajectoryFromCifCore, cif);
    return { trajectory: trajectory.data };
  },
  visuals: defaultVisuals,
});
