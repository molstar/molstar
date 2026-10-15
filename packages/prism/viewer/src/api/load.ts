/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { StateTransform } from '@molstar/core/state';
import { Asset } from '@molstar/core/util/assets';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { OpenFiles } from '@molstar/plugin/state/actions/file';
import { FilesParam, filesOf } from './files.js';
import type { ApiContext } from './area.js';
import { ApiArea, ApiResult } from './area.js';

// `load` (operations.md §2): the core's `files` (drag and drop, Open) now; `action` (Home's registry actions) with the
// scene areas; `entry`, `url`, `mvs`, `svs`, `snapshot` come with the `loaders` feature (plan phase 3).

/** File names a session opens from (`.molx`, `.molj`): opened by `session.open`, not loaded as data. */
const SessionExtensions = ['.molx', '.molj'];

export class LoadApi extends ApiArea {
  readonly id = 'load';

  /**
   * A registry state action applicable to the root (Home's blocks: download a structure, open files, download density;
   * operations.md §2), by its id, with its params. Network actions ask first, from the agent (agent.md §4).
   */
  readonly action = this.operation({
    id: 'load.action',
    title: 'Run a loading action',
    surface: { kind: 'ui-only', reason: "Home's action blocks are its form; the agent uses load.entry and load.url." },
    params: {
      action: PD.Text('', { description: 'A loading action by name: download-structure, download-density, …' }),
      params: PD.Value<Record<string, unknown>>(
        {},
        { description: "The action's params; missing ones are its defaults" },
      ),
    },
    approval: 'network',
    async run({ action, params }, ctx) {
      const root = ctx.plugin.state.data.cells.get(StateTransform.RootRef)!;
      const actions = ctx.plugin.state.data.actions.fromCell(root, ctx.plugin);
      const found = actions.find((a) => actionName(a) === action || a.id === action);
      if (!found) {
        return ApiResult.fail(
          'not-found',
          `There is no loading action '${action}'.`,
          `Actions: ${actions.map(actionName).join(', ')}.`,
        );
      }
      const defs = found.definition.params?.(root.obj!, ctx.plugin) as PD.Params | undefined;
      const values = defs ? { ...PD.getDefaultValues(defs), ...params } : params;
      const before = ctx.plugin.state.data.cells.size;
      await ctx.plugin.runTask(ctx.plugin.state.data.applyAction(found, values));
      const name = found.definition.display.name;
      if (ctx.plugin.state.data.cells.size === before) {
        return ApiResult.fail('failed', `${name} made nothing.`, 'Check its values; the log says more.');
      }
      return ApiResult.ok(`${name}: done.`, { nodes: [] });
    },
  });

  readonly files = this.operation({
    id: 'load.files',
    title: 'Open files',
    params: {
      files: FilesParam({ description: 'Structure, volume and other data files; the format is read from the name' }),
      visuals: PD.Boolean(true, { description: 'Show them with the default representations' }),
    },
    modules: ['files'],
    async run({ files, visuals }, ctx) {
      const list = filesOf(files);
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
      // A lazy extension that reads one of these formats loads first (viewer-spec.md §3).
      await ctx.extensions.ensureFor({ fileExtensions: list.map((f) => extensionOf(f.name)) });
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

/**
 * A state action's name: its display name as a slug (`download-structure`). The plugin's action ids are made per session,
 * so they cannot be what a tool, a command or the run log names.
 */
export function actionName(action: { definition: { display: { name: string } } }): string {
  return action.definition.display.name
    .toLowerCase()
    .replace(/[^a-z0-9]+/g, '-')
    .replace(/^-|-$/g, '');
}

function structureRefs(ctx: ApiContext): string[] {
  return ctx.plugin.managers.structure.hierarchy.current.structures.map((s) => s.cell.transform.ref);
}

/** A file's extension, lowercase, without compression suffixes (`1abc.cif.gz` → `cif`). */
export function extensionOf(name: string): string {
  const parts = name.toLowerCase().split('.');
  while (parts.length > 2 && (parts.at(-1) === 'gz' || parts.at(-1) === 'zip')) parts.pop();
  return parts.length > 1 ? parts.at(-1)! : '';
}
