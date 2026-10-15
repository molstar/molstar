/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Asset } from '@molstar/core/util/assets';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { OpenFiles } from '@molstar/plugin/state/actions/file';
import type { ApiContext } from './area.js';
import { ApiArea, ApiResult } from './area.js';

// `load` (operations.md §2): the core's `files` (drag and drop, Open) now; `action` (Home's registry actions) with the
// scene areas; `entry`, `url`, `mvs`, `svs`, `snapshot` come with the `loaders` feature (plan phase 3).

/** File names a session opens from (`.molx`, `.molj`): opened by `session.open`, not loaded as data. */
const SessionExtensions = ['.molx', '.molj'];

export class LoadApi extends ApiArea {
  readonly id = 'load';

  readonly files = this.operation({
    id: 'load.files',
    title: 'Open files',
    params: {
      files: PD.FileList({ description: 'Structure, volume and other data files; the format is read from the name' }),
      visuals: PD.Boolean(true, { description: 'Show them with the default representations' }),
    },
    modules: ['files'],
    async run({ files, visuals }, ctx) {
      const list = Array.from((files as unknown as ArrayLike<File> | null) ?? []);
      if (!list.length)
        return ApiResult.fail('invalid', 'No files to open.', 'Drop files on the viewer, or pick some.');
      const sessions = list.filter((f) => SessionExtensions.some((e) => f.name.toLowerCase().endsWith(e)));
      if (sessions.length) {
        return ApiResult.fail(
          'invalid',
          `${sessions[0].name} is a session.`,
          'Open it with session.open: it replaces what is loaded.',
        );
      }
      const before = structureRefs(ctx);
      const cellsBefore = ctx.plugin.state.data.cells.size;
      await ctx.plugin.runTask(
        ctx.plugin.state.data.applyAction(OpenFiles, {
          files: list.map((f) => Asset.File(f)),
          format: { name: 'auto', params: {} },
          visuals,
        }),
      );
      const added = structureRefs(ctx).filter((r) => !before.includes(r));
      const names = list.map((f) => f.name).join(', ');
      if (ctx.plugin.state.data.cells.size === cellsBefore) {
        return ApiResult.fail(
          'failed',
          `Could not open ${names}.`,
          'Check the format: the viewer reads it from the file name.',
        );
      }
      const structures = added.length ? ` (${added.length} structure${added.length === 1 ? '' : 's'})` : '';
      return ApiResult.ok(`Opened ${names}${structures}.`, { nodes: added });
    },
  });
}

function structureRefs(ctx: ApiContext): string[] {
  return ctx.plugin.managers.structure.hierarchy.current.structures.map((s) => s.cell.transform.ref);
}
