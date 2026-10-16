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
/** What the network operations declare (their definitions' `modules`, and the remote states' `addresses`). */
const network = { modules: ['network'] as const };
const server = { addresses: ({ server }: { server?: string }) => (server?.trim() ? [server] : []) };

export const LazySuperposition: LazyViewerExtension = {
  id: 'superposition',
  title: 'Superposition',
  builtIn: true,
  stub: {
    parts: [
      {
        id: 'superposition.section.superpose',
        kind: 'section',
        slot: 'selection',
        title: 'Superpose',
        icon: 'layers',
        order: 3,
      },
    ],
    operations: [
      ...ops([
        ['superpose.chains', 'Superpose by chains'],
        ['superpose.tm-align', 'Superpose by TM-align'],
        ['superpose.ligands', 'Superpose by ligands'],
        ['superpose.atoms', 'Superpose by atoms'],
        ['superpose.db', 'Superpose by UniProt (SIFTS)'],
        ['superpose.remove', 'Remove a superposition'],
      ]),
      { id: 'superpose.list', title: 'List the superpositions', read: true },
    ],
    commands: [
      { id: 'superpose.open', title: 'Superpose structures…', category: 'Superposition', icon: 'layers' },
      ...ops([
        ['superpose.chains', 'Superpose by chains'],
        ['superpose.tm-align', 'Superpose by TM-align'],
        ['superpose.ligands', 'Superpose by ligands'],
        ['superpose.atoms', 'Superpose by atoms'],
        ['superpose.db', 'Superpose by UniProt (SIFTS)'],
      ]).map((c) => ({ ...c, category: 'Superposition', icon: 'layers' as const })),
      { id: 'superpose.remove', title: 'Remove a superposition…', category: 'Superposition', icon: 'rotate-ccw' },
    ],
  },
  load: async () => (await import('../superposition/index.js')).SuperpositionFeature,
};

export const LazyStateTree: LazyViewerExtension = {
  id: 'state-tree',
  title: 'State tree',
  builtIn: true,
  stub: {
    parts: [{ id: 'state-tree.tab.tree', kind: 'tab', title: 'State tree', icon: 'list-tree', order: 2 }],
    operations: [
      { id: 'state.get', title: 'Read the state tree', read: true },
      { id: 'state.select', title: 'Find state cells', read: true },
      { id: 'state.visibility', title: 'Show or hide state cells' },
      { id: 'state.update', title: 'Update state cells' },
      { id: 'state.apply', title: 'Apply a transform' },
      { id: 'state.action', title: 'Apply an action to a state cell', approval: 'network' },
      { id: 'state.remove', title: 'Remove state cells', approval: 'destructive' },
    ],
    commands: [
      { id: 'state.visibility', title: 'Show or hide state cells', category: 'State', icon: 'eye' },
      { id: 'state.remove', title: 'Remove state cells', category: 'State', icon: 'trash-2' },
      { id: 'state-tree.open', title: 'Show the state tree', category: 'State', icon: 'list-tree' },
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
      {
        id: 'remote-states.section.home',
        kind: 'section',
        slot: 'home.import',
        title: 'Remote states',
        icon: 'cloud',
        order: 10,
      },
      { id: 'remote-states.overlay.states', kind: 'overlay', title: 'Remote states' },
      { id: 'remote-states.toolbar.states', kind: 'toolbar', title: 'Remote states' },
    ],
    operations: [
      {
        id: 'session.remote.list',
        title: 'List the remote states',
        read: true,
        ...network,
        ...server,
        approval: 'network',
      },
      { id: 'session.remote.open', title: 'Open a remote state', ...network, ...server, approval: 'destructive' },
      {
        id: 'session.remote.upload',
        title: 'Upload the state',
        history: false,
        undoable: false,
        ...network,
        ...server,
        approval: 'upload',
      },
      {
        id: 'session.remote.remove',
        title: 'Remove a remote state',
        undoable: false,
        ...network,
        ...server,
        approval: 'destructive',
      },
    ],
    commands: [
      { id: 'remote-states.open', title: 'Show the remote states', category: 'Remote states', icon: 'cloud' },
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
