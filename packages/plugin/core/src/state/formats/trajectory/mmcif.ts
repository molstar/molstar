/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Adam Midlik <midlik@gmail.com>
 * @author Ludovic Autin <ludovic.autin@gmail.com>
 */

import { PluginStateTransform, PluginStateObject as SO } from '@molstar/plugin/state/objects';
import { Task } from '@molstar/core/task';
import { Model, ArrayTrajectory, type Trajectory } from '@molstar/model/model/structure';
import { trajectoryFromMmCIF, trajectoryFromCCD } from '@molstar/model/formats/structure/mmcif';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { trajectoryProps } from './helpers.js';
import { type TrajectoryFormatProvider, defaultVisuals } from './provider.js';
import { TrajectoryFormatCategory } from './category.js';
import { guessCifVariant, applyTransformerRaw, rawDataObject } from '@molstar/plugin/state/formats/provider';
import { ParseCif } from '@molstar/plugin/state/formats/cif';

export { TrajectoryFromBlob };
type TrajectoryFromBlob = typeof TrajectoryFromBlob;
const TrajectoryFromBlob = PluginStateTransform.BuiltIn({
  name: 'trajectory-from-blob',
  display: { name: 'Parse Blob', description: 'Parse format blob into a single trajectory.' },
  from: SO.Format.Blob,
  to: SO.Molecule.Trajectory,
})({
  apply({ a }) {
    return Task.create('Parse Format Blob', async (ctx) => {
      const models: Model[] = [];
      for (const e of a.data) {
        if (e.kind !== 'cif') continue;
        const block = e.data.blocks[0];
        const xs = await trajectoryFromMmCIF(block).runInContext(ctx);
        if (xs.frameCount === 0) throw new Error('No models found.');

        for (let i = 0; i < xs.frameCount; i++) {
          const x = await Task.resolveInContext(xs.getFrameAtIndex(i), ctx);
          models.push(x);
        }
      }

      for (let i = 0; i < models.length; i++) {
        Model.TrajectoryInfo.set(models[i], { index: i, size: models.length });
      }

      const props = { label: 'Trajectory', description: `${models.length} model${models.length === 1 ? '' : 's'}` };
      return new SO.Molecule.Trajectory(new ArrayTrajectory(models), props);
    });
  },
});

export { TrajectoryFromMmCif };
type TrajectoryFromMmCif = typeof TrajectoryFromMmCif;
const TrajectoryFromMmCif = PluginStateTransform.BuiltIn({
  name: 'trajectory-from-mmcif',
  display: {
    name: 'Trajectory from mmCIF',
    description: 'Identify and create all separate models in the specified CIF data block',
  },
  from: SO.Format.Cif,
  to: SO.Molecule.Trajectory,
  params(a) {
    if (!a) {
      return {
        loadAllBlocks: PD.Optional(
          PD.Boolean(false, {
            description:
              'If True, ignore Block Header and Block Index parameters and parse all datablocks into a single trajectory.',
          }),
        ),
        blockHeader: PD.Optional(
          PD.Text(void 0, {
            description: 'Header of the block to parse. If not specifed, Block Index parameter applies.',
            hideIf: (p) => p.loadAllBlocks === true,
          }),
        ),
        blockIndex: PD.Optional(
          PD.Numeric(
            0,
            { min: 0, step: 1 },
            {
              description:
                'Zero-based index of the block to parse. Only applies when Block Header parameter is not specified.',
              hideIf: (p) => p.loadAllBlocks === true || p.blockHeader,
            },
          ),
        ),
      };
    }
    const { blocks } = a.data;
    const headers = blocks.map((b) => [b.header, b.header] as [string, string]);
    headers.push(['', '[Use Block Index]']);
    return {
      loadAllBlocks: PD.Optional(
        PD.Boolean(false, {
          description:
            'If True, ignore Block Header and Block Index parameters and parse all data blocks into a single trajectory.',
        }),
      ),
      blockHeader: PD.Optional(
        PD.Select(blocks[0] && blocks[0].header, headers, {
          description: 'Header of the block to parse. If not specifed, Block Index parameter applies.',
          hideIf: (p) => p.loadAllBlocks === true,
        }),
      ),
      blockIndex: PD.Optional(
        PD.Numeric(
          0,
          { min: 0, step: 1, max: blocks.length - 1 },
          {
            description:
              'Zero-based index of the block to parse. Only applies when Block Header parameter is not specified.',
            hideIf: (p) => p.loadAllBlocks === true || p.blockHeader,
          },
        ),
      ),
    };
  },
})({
  isApplicable: (a) => a.data.blocks.length > 0,
  apply({ a, params }) {
    return Task.create('Parse mmCIF', async (ctx) => {
      let trajectory: Trajectory;
      if (params.loadAllBlocks) {
        const models: Model[] = [];
        for (const block of a.data.blocks) {
          if (ctx.shouldUpdate) {
            await ctx.update(`Parsing ${block.header}...`);
          }
          const t = await trajectoryFromMmCIF(block).runInContext(ctx);
          for (let i = 0; i < t.frameCount; i++) {
            models.push(await Task.resolveInContext(t.getFrameAtIndex(i), ctx));
          }
        }
        trajectory = new ArrayTrajectory(models);
      } else {
        const header = params.blockHeader || a.data.blocks[params.blockIndex ?? 0].header;
        const block = a.data.blocks.find((b) => b.header === header);
        if (!block) throw new Error(`Data block '${[header]}' not found.`);
        const isCcd =
          block.categoryNames.includes('chem_comp_atom') &&
          !block.categoryNames.includes('atom_site') &&
          !block.categoryNames.includes('ihm_sphere_obj_site') &&
          !block.categoryNames.includes('ihm_gaussian_obj_site');
        trajectory = isCcd
          ? await trajectoryFromCCD(block).runInContext(ctx)
          : await trajectoryFromMmCIF(block, a.data).runInContext(ctx);
      }
      if (trajectory.frameCount === 0) throw new Error('No models found.');
      const props = trajectoryProps(trajectory);
      return new SO.Molecule.Trajectory(trajectory, props);
    });
  },
});

export const MmcifProvider: TrajectoryFormatProvider = {
  label: 'mmCIF',
  description: 'mmCIF',
  category: TrajectoryFormatCategory,
  stringExtensions: ['cif', 'mmcif', 'mcif'],
  binaryExtensions: ['bcif'],
  isApplicable: (info, data) => {
    if (info.ext === 'mmcif' || info.ext === 'mcif') return true;
    // assume undetermined cif/bcif files are mmCIF
    if (info.ext === 'cif' || info.ext === 'bcif') return guessCifVariant(info, data) === -1;
    return false;
  },
  parse: async (plugin, data, params) => {
    const state = plugin.state.data;
    const cif = state
      .build()
      .to(data)
      .apply(ParseCif, void 0, { state: { isGhost: true } });
    const trajectory = await cif
      .apply(TrajectoryFromMmCif, void 0, { tags: params?.trajectoryTags })
      .commit({ revertOnError: true });

    if ((cif.selector.cell?.obj?.data.blocks.length || 0) > 1) {
      plugin.state.data.updateCellState(cif.ref, { isGhost: false });
    }

    return { trajectory };
  },
  parseRaw: async (plugin, ctx, data) => {
    const cif = await applyTransformerRaw(plugin, ctx, ParseCif, rawDataObject(data));
    const trajectory = await applyTransformerRaw(plugin, ctx, TrajectoryFromMmCif, cif);
    return { trajectory: trajectory.data };
  },
  visuals: defaultVisuals,
};
