/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { ParamsFormModel } from '../forms/params-form.js';
import type { HomeAction, HomeModel, RemoteState, RemoteStatesModel } from '../model/home.js';
import { observableState } from '../reactive/observable-state.js';

// The root's state actions, shaped like Mol*'s (DownloadStructure, DownloadDensity, LoadTrajectory, OpenSession): the
// source as a mapped param, each source with its own fields.
const DownloadStructure = {
  source: PD.MappedStatic(
    'pdb',
    {
      pdb: PD.Group(
        {
          id: PD.Text('1tqn', { label: 'PDB ID', isEssential: true }),
          server: PD.Select('pdbe', [
            ['pdbe', 'PDBe'],
            ['rcsb', 'RCSB PDB'],
            ['pdbj', 'PDBj'],
          ]),
        },
        { isFlat: true },
      ),
      'pdb-ihm': PD.Group({ id: PD.Text('8zzc', { label: 'PDB-IHM ID', isEssential: true }) }, { isFlat: true }),
      alphafolddb: PD.Group({ id: PD.Text('Q8W3K0', { label: 'UniProt ID', isEssential: true }) }, { isFlat: true }),
      modelarchive: PD.Group(
        { id: PD.Text('ma-bak-cepc-0003', { label: 'ModelArchive ID', isEssential: true }) },
        { isFlat: true },
      ),
      pubchem: PD.Group({ id: PD.Text('2244', { label: 'PubChem CID', isEssential: true }) }, { isFlat: true }),
      url: PD.Group(
        {
          url: PD.Text('', { label: 'URL', isEssential: true }),
          format: PD.Select('mmcif', [
            ['mmcif', 'mmCIF'],
            ['bcif', 'BinaryCIF'],
            ['pdb', 'PDB'],
            ['gro', 'GRO'],
          ]),
        },
        { isFlat: true },
      ),
    },
    {
      label: 'Source',
      isEssential: true,
      options: [
        ['pdb', 'PDB'],
        ['pdb-ihm', 'PDB-IHM'],
        ['alphafolddb', 'AlphaFold DB'],
        ['modelarchive', 'ModelArchive'],
        ['pubchem', 'PubChem'],
        ['url', 'URL'],
      ],
    },
  ),
  preset: PD.Select(
    'default',
    [
      ['default', 'Default'],
      ['all-models', 'All models'],
      ['unit-cell', 'Unit cell'],
      ['supercell', 'Supercell'],
    ],
    { label: 'Preset' },
  ),
};

const DownloadDensity = {
  source: PD.MappedStatic(
    'pdbe',
    {
      pdbe: PD.Group(
        {
          id: PD.Text('1tqn', { label: 'PDB ID', isEssential: true }),
          type: PD.Select('2fofc', [
            ['2fofc', '2Fo-Fc'],
            ['fofc', 'Fo-Fc'],
          ]),
        },
        { isFlat: true },
      ),
      emdb: PD.Group({ id: PD.Text('EMD-2984', { label: 'EMDB ID', isEssential: true }) }, { isFlat: true }),
      url: PD.Group(
        {
          url: PD.Text('', { label: 'URL', isEssential: true }),
          format: PD.Select('ccp4', [
            ['ccp4', 'CCP4/MRC'],
            ['dscif', 'DensityServer CIF'],
            ['dx', 'DX'],
            ['cube', 'Cube'],
          ]),
        },
        { isFlat: true },
      ),
    },
    {
      label: 'Source',
      isEssential: true,
      options: [
        ['pdbe', 'PDBe (X-ray)'],
        ['emdb', 'EMDB'],
        ['url', 'URL'],
      ],
    },
  ),
};

const LoadTrajectory = {
  model: PD.Text('', { label: 'Topology URL', isEssential: true }),
  coordinates: PD.Text('', { label: 'Coordinates URL', isEssential: true }),
  format: PD.Select('xtc', [
    ['xtc', 'XTC'],
    ['dcd', 'DCD'],
    ['trr', 'TRR'],
    ['nctraj', 'NetCDF'],
  ]),
};

// Mol*'s OpenFiles: the files, their format (detected by default), default visuals.
const OpenFiles = {
  files: PD.FileList({
    label: 'Files',
    accept: '.cif,.bcif,.pdb,.ent,.gro,.ccp4,.mrc,.map,.molx,.molj,.mvsj,.gz,.zip',
    multiple: true,
  }),
  format: PD.Select('auto', [
    ['auto', 'Auto'],
    ['mmcif', 'mmCIF'],
    ['bcif', 'BinaryCIF'],
    ['pdb', 'PDB'],
    ['gro', 'GRO'],
    ['ccp4', 'CCP4/MRC'],
    ['molx', 'Mol* session'],
    ['mvsj', 'MolViewSpec'],
  ]),
  visuals: PD.Boolean(true, { description: 'Add default visuals' }),
};

const OpenSession = {
  url: PD.Text('', { label: 'Session URL', isEssential: true }),
};

/** Home for the stories: the root actions with forms; applying posts what would happen. */
export function mockHome(
  post: (text: string, kind?: 'info' | 'success') => void,
  /** `remoteStates`: on, or on with this many states (a long list). */
  options: { remoteStates?: boolean | number } = {},
): HomeModel {
  const busy = observableState(false);
  const action = (
    id: string,
    title: string,
    icon: HomeAction['icon'],
    description: string,
    params: PD.Params,
    applyLabel: string,
  ): HomeAction => {
    const form = new ParamsFormModel(params);
    return {
      id,
      title,
      icon,
      description,
      form,
      applyLabel,
      busy: busy.asReadonly,
      apply: () => {
        const v = form.values.value as any;
        const what = v.source
          ? `${v.source.name}: ${v.source.params.id ?? v.source.params.url ?? ''}`
          : v.files
            ? v.files.map((f: { name: string }) => f.name).join(', ')
            : title;
        busy.set(true);
        post(`${applyLabel === 'Download' ? 'Downloading' : 'Opening'} ${what}…`);
        setTimeout(() => {
          busy.set(false);
          post(`${title}: done (${what})`, 'success');
        }, 1200);
      },
    };
  };
  const actions = observableState<readonly HomeAction[]>([
    action(
      'download-structure',
      'Download Structure',
      'download',
      'PDB, AlphaFold DB, URL, …',
      DownloadStructure,
      'Download',
    ),
    action('open-files', 'Open Files', 'folder-open', 'mmCIF, PDB, CCP4, MOLX, …', OpenFiles, 'Open'),
    action('download-density', 'Download Density', 'volume', 'PDBe, EMDB, URL', DownloadDensity, 'Download'),
    action('load-trajectory', 'Load Trajectory', 'trajectory', 'Topology and coordinates', LoadTrajectory, 'Load'),
    action('open-session', 'Open Session', 'folder-open', 'A .molx session by URL', OpenSession, 'Load'),
  ]);
  return {
    actions: actions.asReadonly,
    formats: ['mmCIF', 'BinaryCIF', 'PDB', 'CCP4', 'MRC', 'MOLX', 'MVSJ'],
    remoteStates: options.remoteStates
      ? mockRemoteStates(post, typeof options.remoteStates === 'number' ? options.remoteStates : undefined)
      : undefined,
    openFiles: (files) => post(`Opening ${files.map((f) => f.name).join(', ')}…`),
  };
}

/** A state server with a few saved sessions (Mol*'s Remote States). */
function mockRemoteStates(post: (text: string) => void, count?: number): RemoteStatesModel {
  const serverUrl = observableState('https://webchem.ncbr.muni.cz/molstar-state');
  const loading = observableState(false);
  const few: RemoteState[] = [
    { id: '1', label: 'Zika + EM', date: 'Sun, 11 Oct 2026 09:05' },
    { id: '2', label: 'P-450 superposition & validation', date: 'Fri, 9 Oct 2026 16:40' },
    { id: '3', label: 'NPC', date: 'Mon, 5 Oct 2026 11:12' },
    { id: '4', label: '1RB8 assembly symmetry', date: 'Thu, 1 Oct 2026 08:30' },
  ];
  const states = observableState<readonly RemoteState[]>(
    count === undefined
      ? few
      : Array.from({ length: count }, (_, i) => ({
          id: `${i + 1}`,
          label: `${few[i % few.length].label}${i >= few.length ? ` (${Math.floor(i / few.length) + 1})` : ''}`,
          date: `${['Sun', 'Sat', 'Fri', 'Thu', 'Wed', 'Tue', 'Mon'][i % 7]}, ${Math.max(1, 28 - (i % 28))} Sep 2026 ${String(8 + (i % 10)).padStart(2, '0')}:${String((i * 7) % 60).padStart(2, '0')}`,
        })),
  );
  return {
    serverUrl: serverUrl.asReadonly,
    setServerUrl: (url) => {
      serverUrl.set(url);
      post(`State server: ${url}`);
    },
    states: states.asReadonly,
    loading: loading.asReadonly,
    refresh: () => {
      loading.set(true);
      setTimeout(() => loading.set(false), 800);
    },
    open: (id) => post(`Opening "${states.value.find((s) => s.id === id)?.label}"…`),
    upload: (name) => {
      states.set([{ id: `${Date.now()}`, label: name, date: 'just now' }, ...states.value]);
      post(`Uploaded "${name}"`);
    },
    download: (id) => post(`Downloading "${states.value.find((s) => s.id === id)?.label}"…`),
    remove: (id) => {
      post(`Removed "${states.value.find((s) => s.id === id)?.label}"`);
      states.set(states.value.filter((s) => s.id !== id));
    },
  };
}
