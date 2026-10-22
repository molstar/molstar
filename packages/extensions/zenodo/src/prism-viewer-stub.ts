/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { LazyViewerExtension, ViewerSkill } from '@molstar/prism-viewer/spec';

// The Zenodo import as a lazy Prism extension (viewer-spec.md §3, parity.md §6: Home's import slot; lazy): this module
// is the stub, plain data and types only. Its operations and commands exist from the start; the first use loads
// `prism-viewer.ts`, checked against this.

export const ZenodoSection = 'zenodo.section.import';

export const ZenodoCommands = {
  list: 'zenodo.list',
  load: 'zenodo.load',
  open: 'zenodo.open',
} as const;

/** A record's address on Zenodo's API: what `zenodo.list` reaches, and `zenodo.load` before the file. */
export const zenodoRecordAddress = (id: string) => `https://zenodo.org/api/records/${id}`;

/** A file's address: the record's API, by file name. */
export const zenodoFileAddress = (id: string, name: string) =>
  `${zenodoRecordAddress(id)}/files/${encodeURIComponent(name)}/content`;

/**
 * The record id in what a person writes: the id (1234567), a DOI (10.5281/zenodo.1234567) or a record's page
 * (https://zenodo.org/records/1234567); undefined when it is none of these.
 */
export function zenodoRecordId(text: string): string | undefined {
  const t = text.trim();
  if (/^\d+$/.test(t)) return t;
  return (
    t.match(/zenodo\.org\/(?:api\/)?records?\/(\d+)/i)?.[1] ??
    t.match(/^(?:https?:\/\/(?:dx\.)?doi\.org\/)?10\.5281\/zenodo\.(\d+)$/i)?.[1]
  );
}

const recordOf = (record: string) => zenodoRecordId(record) ?? encodeURIComponent(record.trim());

/**
 * The operations, with the flags their definitions have (the host compares them when the extension loads), and the
 * addresses a call reaches (the agent's policy reads them before the load: agent.md §4).
 */
export const ZenodoOperations = {
  list: {
    id: 'zenodo.list',
    title: 'List the files of a Zenodo record',
    read: true,
    approval: 'network',
    modules: ['network'],
    addresses: ({ record }: { record: string }) => [zenodoRecordAddress(recordOf(record))],
  },
  load: {
    id: 'zenodo.load',
    title: 'Load a file from Zenodo',
    approval: 'network',
    modules: ['network'],
    addresses: ({ record, file, coordinates }: { record: string; file: string; coordinates?: string }) => [
      zenodoRecordAddress(recordOf(record)),
      zenodoFileAddress(recordOf(record), file.trim()),
      ...(coordinates?.trim() ? [zenodoFileAddress(recordOf(record), coordinates.trim())] : []),
    ],
  },
} as const;

export const ZenodoSkill: ViewerSkill = {
  id: 'zenodo',
  title: 'Zenodo',
  summary: 'list and load the files of a Zenodo record (by id, DOI or link)',
  guidance: [
    'Skill zenodo: a record is given by its id (1234567), DOI (10.5281/zenodo.1234567) or link.',
    'list_zenodo_files first: it names each file, its size and whether the viewer reads it (structure, volume,',
    'topology, coordinates). load_zenodo_file loads one by name; a topology or model with coordinates gives a trajectory.',
  ].join('\n'),
  tools: [
    { id: 'list_zenodo_files', operation: 'zenodo.list', description: "Lists a Zenodo record's files." },
    { id: 'load_zenodo_file', operation: 'zenodo.load', description: 'Loads a file of a Zenodo record by its name.' },
  ],
};

export const LazyZenodo: LazyViewerExtension = {
  id: 'zenodo',
  title: 'Zenodo',
  stub: {
    parts: [{ id: ZenodoSection, kind: 'section', slot: 'home.import', title: 'Zenodo', icon: 'cloud' }],
    operations: [ZenodoOperations.list, ZenodoOperations.load],
    commands: [
      { id: ZenodoCommands.list, title: 'List the files of a Zenodo record', category: 'Load', icon: 'cloud' },
      { id: ZenodoCommands.load, title: 'Load a file from Zenodo', category: 'Load', icon: 'cloud' },
      { id: ZenodoCommands.open, title: 'Show the Zenodo import', category: 'Load', icon: 'cloud' },
    ],
    skills: [ZenodoSkill],
  },
  load: async () => (await import('./prism-viewer.js')).ZenodoExtension,
};
