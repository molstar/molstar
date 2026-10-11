/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { ParamsFormModel } from '../../forms/params-form.js';
import type { HomeAction, HomeModel, RecentEntry } from '../../model/home.js';
import { observableState } from '../../reactive/observable-state.js';

// The root's state actions, shaped like Mol*'s (DownloadStructure, DownloadDensity, LoadTrajectory, OpenSession): the
// source as a mapped param, each source with its own fields.
const DownloadStructure = {
  source: PD.MappedStatic(
    'pdb',
    {
      pdb: PD.Group({
        id: PD.Text('1tqn', { label: 'PDB ID', isEssential: true }),
        server: PD.Select('pdbe', [
          ['pdbe', 'PDBe'],
          ['rcsb', 'RCSB PDB'],
          ['pdbj', 'PDBj'],
        ]),
      }),
      'pdb-ihm': PD.Group({ id: PD.Text('8zzc', { label: 'PDB-IHM ID', isEssential: true }) }),
      alphafolddb: PD.Group({ id: PD.Text('Q8W3K0', { label: 'UniProt ID', isEssential: true }) }),
      modelarchive: PD.Group({ id: PD.Text('ma-bak-cepc-0003', { label: 'ModelArchive ID', isEssential: true }) }),
      pubchem: PD.Group({ id: PD.Text('2244', { label: 'PubChem CID', isEssential: true }) }),
      url: PD.Group({
        url: PD.Text('', { label: 'URL', isEssential: true }),
        format: PD.Select('mmcif', [
          ['mmcif', 'mmCIF'],
          ['bcif', 'BinaryCIF'],
          ['pdb', 'PDB'],
          ['gro', 'GRO'],
        ]),
      }),
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
      pdbe: PD.Group({
        id: PD.Text('1tqn', { label: 'PDB ID', isEssential: true }),
        type: PD.Select('2fofc', [
          ['2fofc', '2Fo-Fc'],
          ['fofc', 'Fo-Fc'],
        ]),
      }),
      emdb: PD.Group({ id: PD.Text('EMD-2984', { label: 'EMDB ID', isEssential: true }) }),
      url: PD.Group({
        url: PD.Text('', { label: 'URL', isEssential: true }),
        format: PD.Select('ccp4', [
          ['ccp4', 'CCP4/MRC'],
          ['dscif', 'DensityServer CIF'],
          ['dx', 'DX'],
          ['cube', 'Cube'],
        ]),
      }),
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

const OpenSession = {
  url: PD.Text('', { label: 'Session URL', isEssential: true }),
};

/** Home for the stories: the root actions with forms, a few recent entries; applying posts what would happen. */
export function mockHome(post: (text: string, kind?: 'info' | 'success') => void): HomeModel {
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
        const what = v.source ? `${v.source.name}: ${v.source.params.id ?? v.source.params.url ?? ''}` : title;
        busy.set(true);
        post(`${applyLabel === 'Load' ? 'Loading' : 'Downloading'} ${what}…`);
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
    action('download-density', 'Download Density', 'volume', 'PDBe, EMDB, URL', DownloadDensity, 'Download'),
    action('load-trajectory', 'Load Trajectory', 'trajectory', 'Topology and coordinates', LoadTrajectory, 'Load'),
    action('open-session', 'Open Session', 'folder-open', 'A .molx session by URL', OpenSession, 'Load'),
  ]);
  const recent = (id: string, label: string, description: string, icon: RecentEntry['icon']): RecentEntry => ({
    id,
    label,
    description,
    icon,
    open: () => post(`Opening ${label}…`),
  });
  return {
    actions: actions.asReadonly,
    recent: observableState<readonly RecentEntry[]>([
      recent('1tqn', '1TQN', 'PDB · Cytochrome P450 3A4 · 2 min ago', 'atom'),
      recent('emd-2984', 'EMD-2984', 'EMDB · β-galactosidase · yesterday', 'volume'),
      recent('af-q8w3k0', 'AF-Q8W3K0-F1', 'AlphaFold DB · last week', 'atom'),
    ]).asReadonly,
    formats: ['mmCIF', 'BinaryCIF', 'PDB', 'CCP4', 'MRC', 'MOLX', 'MVSJ'],
    openFiles: (files) => post(`Opening ${files.map((f) => f.name).join(', ')}…`),
  };
}
