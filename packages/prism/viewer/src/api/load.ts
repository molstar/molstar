/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { StateSelection, StateTransform } from '@molstar/core/state';
import { Asset } from '@molstar/core/util/assets';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { OpenFiles } from '@molstar/plugin/state/actions/file';
import { ReadFile } from '@molstar/plugin/state/transforms/data/fetch';
import { FilesParam, filesOf } from './files.js';
import { extensionOf } from './input.js';
import type { ApiContext } from './area.js';
import { ApiArea, ApiResult } from './area.js';

// `load` (api.md §2): the core's `files` (drag and drop, Open); `action` (Home's registry actions) with the
// scene areas; `entry`, `url`, `mvs`, `svs`, `snapshot` come with the `loaders` feature.

/** File names a session opens from (`.molx`, `.molj`): opened by `session.open`, not loaded as data. */
const SessionExtensions = ['.molx', '.molj'];
/** SeqViewSpec documents: shown in the sequence by `load.svs` (the loaders feature), after the other files. */
const SvsExtensions = ['svs', 'svsj', 'svsx'];

export class LoadApi extends ApiArea {
  readonly id = 'load';

  /**
   * A registry state action applicable to the root (Home's blocks: download a structure, open files, download density;
   * api.md §2), by its id, with its params. Network actions ask first, from the agent (agent.md §4).
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
      // SeqViewSpec documents go to `load.svs` once the rest is open, so one may name a file dropped with it.
      const documents = list.filter((f) => SvsExtensions.includes(extensionOf(f.name)));
      const data = list.filter((f) => !documents.includes(f));
      if (documents.length && !ctx.api.describe('load.svs')) {
        return ApiResult.fail(
          'unavailable',
          `${documents[0].name} is a SeqViewSpec, which this viewer does not load.`,
          'Its loaders feature (load.svs) is not in the viewer.',
        );
      }
      const before = structureRefs(ctx);
      const opened: string[] = [];
      const warnings: string[] = [];
      if (data.length) {
        // A lazy extension that reads one of these formats loads first (viewer-spec.md §3).
        await ctx.extensions.ensureFor({ fileExtensions: data.map((f) => extensionOf(f.name)) });
        // One data transaction for all of them: the viewer's sources (the scene, the structures' lists, the commands)
        // sync once at its end, not once a file over every structure loaded so far (files × structures).
        await ctx.plugin.dataTransaction(async () => {
          for (const f of data) {
            const outcome = await openFile(ctx, f, visuals);
            if (outcome === 'opened') opened.push(f.name);
            else if (outcome === 'partly') {
              opened.push(f.name);
              warnings.push(`Some of ${f.name} could not be opened.`);
            } else warnings.push(`Could not open ${f.name}.`);
          }
        });
        if (!opened.length) {
          return ApiResult.fail(
            'failed',
            `Could not open ${data.map((f) => f.name).join(', ')}.`,
            'Check the format: the viewer reads it from the file name.',
          );
        }
      }
      const names = opened.join(', ');
      for (const doc of documents) {
        const r = await ctx.api.run('load.svs', { file: doc });
        if (!r.ok) return r;
      }
      const added = structureRefs(ctx).filter((r) => !before.includes(r));
      const structures = added.length ? ` (${added.length} structure${added.length === 1 ? '' : 's'})` : '';
      const shown = documents.map((f) => f.name).join(', ');
      const did = !names
        ? `The sequence shows ${shown}.`
        : `Opened ${names}${structures}${shown ? `; the sequence shows ${shown}` : ''}.`;
      return ApiResult.ok(
        did,
        { nodes: added, ...(documents.length ? { sequence: true } : {}) },
        warnings.length ? { warnings } : {},
      );
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

/**
 * Opens one file with `OpenFiles` and judges it by what it made: the plugin logs a file it cannot parse and keeps going,
 * leaving the raw data it read (a `ReadFile` node with nothing under it, or a failed node under it). Those are deleted,
 * so a failed file leaves no node behind; a zip's entries are judged one by one (`partly` when only some opened).
 */
async function openFile(ctx: ApiContext, file: File, visuals: boolean): Promise<'opened' | 'partly' | 'failed'> {
  const state = ctx.plugin.state.data;
  const rootChildren = () => state.tree.children.get(StateTransform.RootRef).toArray();
  const before = new Set(rootChildren());
  await ctx.plugin.runTask(
    state.applyAction(OpenFiles, {
      files: [Asset.File(file)],
      format: { name: 'auto', params: {} },
      visuals,
    }),
  );
  const made = rootChildren().filter((r) => !before.has(r));
  const failed = made.filter((ref) => {
    if (state.cells.get(ref)?.transform.transformer !== ReadFile) return false;
    const subtree = state.select(StateSelection.Generators.byRef(ref).subtree());
    return subtree.length < 2 || subtree.some((c) => c.status === 'error');
  });
  if (failed.length) {
    const b = state.build();
    for (const ref of failed) b.delete(ref);
    await ctx.plugin.runTask(state.updateTree(b));
  }
  if (failed.length === made.length) return 'failed';
  return failed.length ? 'partly' : 'opened';
}

function structureRefs(ctx: ApiContext): string[] {
  return ctx.plugin.managers.structure.hierarchy.current.structures.map((s) => s.cell.transform.ref);
}

export { extensionOf };
