/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import type { CommandDefinition } from '@molstar/prism-ui/core';
import { ApiResult, type OperationFn, type ViewerApi } from '../../api/area.js';
import type { ViewerCommandContext } from '../../api/managers.js';
import { operationCommand, type OperationCommandOptions, runAsCommand } from '../../commands/operations.js';
import type { ViewerExtension } from '../../spec.js';
import { SuperposeApi } from './operations.js';

declare module '../../api.js' {
  // The `superpose` area (api.md §2): superposing structures by chains, TM-align, ligands, atoms or UniProt.
  interface ExtensionApis {
    superpose: SuperposeApi;
  }
}

// The `superposition` feature (viewer-spec.md §3.1): the `superpose` area and the Superpose section of the Selection
// tab (parity.md §2). The views are a chunk of their own.

/** A list param (targets, structures) read as text in a command's form: items separated by `;`. */
function listCommand(
  api: ViewerApi,
  op: OperationFn<any, any>,
  list: 'targets' | 'structures',
  o: OperationCommandOptions,
  what: string,
): CommandDefinition<ViewerCommandContext> {
  const base = operationCommand(api, op, o);
  const split = (v: unknown) =>
    String(v ?? '')
      .split(';')
      .map((x) => x.trim())
      .filter(Boolean);
  return {
    ...base,
    fields: undefined,
    params: (ctx, values) => ({
      ...(base.params?.(ctx, values) ?? {}),
      [list]: {
        ...PD.Text('', {
          label: list === 'targets' ? 'Targets' : 'Structures',
          description: `${what}, separated by ';', the first is the reference. Empty: ${list === 'targets' ? 'the selection of each structure' : 'every structure'}.`,
        }),
        isOptional: true,
      },
    }),
    async run(args) {
      const r = await runAsCommand(api, op.definition.id, { ...args, [list]: split(args[list]) });
      return { ...r, undoable: true };
    },
  };
}

const SelectionTab = 'core.tab.selection';

export const SuperpositionFeature: ViewerExtension = {
  id: 'superposition',
  title: 'Superposition',
  builtIn: true,
  api: { superpose: (ctx) => new SuperposeApi(ctx) },
  commands: (api) => {
    const s = api.superpose;
    if (!s) return [];
    const category = 'Superposition';
    return [
      {
        id: 'superpose.open',
        title: 'Superpose structures…',
        category,
        icon: 'layers',
        alias: ['superpose', 'superposition', 'align structures', 'rmsd'],
        // The section is in the Selection panel: without the panel there is nowhere to open.
        when: () => (api.view.hasTab(SelectionTab) ? true : 'There is no Selection panel.'),
        // Opens the section: a view, not a change.
        async run() {
          const r = await api.view.open({ tab: SelectionTab, section: 'superposition.section.superpose' });
          if (ApiResult.isFailure(r)) throw new Error(r.message);
          return {
            message: 'Opened Superpose.',
            ops: 'view.open({ tab: "core.tab.selection", section: "superposition.section.superpose" })',
          };
        },
      },
      listCommand(
        api,
        s.chains,
        'targets',
        { category, icon: 'layers', alias: ['superpose chains'], action: 'Superpose' },
        'One polymer chain per structure',
      ),
      listCommand(
        api,
        s.tmAlign,
        'targets',
        { category, icon: 'layers', alias: ['tm-align', 'tm align', 'tmalign'], action: 'Superpose' },
        'One polymer chain per structure',
      ),
      listCommand(
        api,
        s.ligands,
        'targets',
        { category, icon: 'layers', alias: ['superpose ligands', 'mccs'], action: 'Superpose' },
        'One ligand per structure',
      ),
      listCommand(
        api,
        s.atoms,
        'targets',
        { category, icon: 'layers', alias: ['superpose atoms'], action: 'Superpose' },
        'The same atoms in each structure',
      ),
      listCommand(
        api,
        s.db,
        'structures',
        { category, icon: 'layers', alias: ['superpose uniprot', 'sifts'], action: 'Superpose' },
        'The structures',
      ),
      operationCommand(api, s.remove, {
        category,
        icon: 'rotate-ccw',
        title: 'Remove a superposition…',
        alias: ['remove superposition', 'undo superposition', 'reset superposition', 'unsuperpose'],
        action: 'Remove',
      }),
    ];
  },
  ui: async (ctx) => (await import('./ui.js')).superpositionUI(ctx),
};
