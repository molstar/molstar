/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { LazyViewerExtension } from '../../spec.js';

// The heavier features as lazy stubs (viewer-spec.md §3.1, review F44): each declares every part, operation and command
// it brings, and its `load` imports the feature's module, so importing this list pulls in none of their code. An app
// that wants one up front lists the feature itself, or sets `preload: true`.

const ops = (list: readonly [string, string][]) => list.map(([id, title]) => ({ id, title }));

export const LazySuperposition: LazyViewerExtension = {
  id: 'superposition',
  title: 'Superposition',
  builtIn: true,
  stub: {
    parts: [{ id: 'superposition.section.superpose', kind: 'section', slot: 'selection', title: 'Superpose' }],
    operations: ops([
      ['superpose.chains', 'Superpose by chains'],
      ['superpose.tm-align', 'Superpose by TM-align'],
      ['superpose.ligands', 'Superpose by ligands'],
      ['superpose.atoms', 'Superpose by atoms'],
      ['superpose.db', 'Superpose by UniProt (SIFTS)'],
    ]),
    commands: [
      { id: 'superpose.open', title: 'Superpose structures…', category: 'Superposition', icon: 'layers' },
      ...ops([
        ['superpose.chains', 'Superpose by chains'],
        ['superpose.tm-align', 'Superpose by TM-align'],
        ['superpose.ligands', 'Superpose by ligands'],
        ['superpose.atoms', 'Superpose by atoms'],
        ['superpose.db', 'Superpose by UniProt (SIFTS)'],
      ]).map((c) => ({ ...c, category: 'Superposition', icon: 'layers' as const })),
    ],
  },
  load: async () => (await import('../superposition/index.js')).SuperpositionFeature,
};

export const LazyStateTree: LazyViewerExtension = {
  id: 'state-tree',
  title: 'State tree',
  builtIn: true,
  stub: {
    parts: [{ id: 'state-tree.tab.tree', kind: 'tab', title: 'State tree' }],
    operations: ops([
      ['state.get', 'Read the state tree'],
      ['state.visibility', 'Show or hide a state cell'],
      ['state.update', 'Update a state cell'],
      ['state.action', 'Apply an action to a state cell'],
      ['state.remove', 'Remove a state cell'],
    ]),
    commands: [
      { id: 'state.visibility', title: 'Show or hide a state cell', category: 'State', icon: 'eye' },
      { id: 'state.remove', title: 'Remove a state cell', category: 'State', icon: 'trash-2' },
      { id: 'state-tree.show', title: 'Show the state tree', category: 'State', icon: 'list-tree' },
    ],
  },
  load: async () => (await import('../state-tree/index.js')).StateTreeFeature,
};

export const LazyRemoteStates: LazyViewerExtension = {
  id: 'remote-states',
  title: 'Remote states',
  builtIn: true,
  stub: {
    parts: [
      { id: 'remote-states.overlay.states', kind: 'overlay', title: 'Remote states' },
      { id: 'remote-states.toolbar.states', kind: 'toolbar', title: 'Remote states' },
    ],
    operations: ops([
      ['session.remote.list', 'List the remote states'],
      ['session.remote.open', 'Open a remote state'],
      ['session.remote.upload', 'Upload the state'],
      ['session.remote.remove', 'Remove a remote state'],
    ]),
    commands: [
      { id: 'remote-states.show', title: 'Show the remote states', category: 'Remote states', icon: 'cloud' },
      ...ops([
        ['session.remote.list', 'List the remote states'],
        ['session.remote.open', 'Open a remote state'],
        ['session.remote.upload', 'Upload the state'],
        ['session.remote.remove', 'Remove a remote state'],
      ]).map((c) => ({ ...c, category: 'Remote states', icon: 'cloud' as const })),
    ],
  },
  load: async () => (await import('../remote-states/index.js')).RemoteStatesFeature,
};

/** The lazy features, by id. */
export const LazyFeatures = {
  superposition: LazySuperposition,
  'state-tree': LazyStateTree,
  'remote-states': LazyRemoteStates,
} as const;
