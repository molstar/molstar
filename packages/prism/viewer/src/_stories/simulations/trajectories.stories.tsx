/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { story } from '@molstar/prism-ui/testing/index';
import type { PrismViewer } from '../../viewer.js';
import { ExampleFiles, openAt, RealViewerStory } from '../../testing/real-viewer.js';

// MD trajectories in the real viewer, one story per coordinates format: a model or topology and its coordinates
// (`load.trajectory`), or a file of frames (`load.url`), with the playback card over the frames. PDB + DCD and
// PRMTOP + NetCDF are the repository's (`data/examples/trajectory`); the rest are mdtraj's test data (one alanine
// dipeptide run in every format, LGPL, not copied here), fetched from GitHub: they reach the network.

const meta: Meta<Args> = {
  title: 'Viewer/Simulations/Trajectories',
  parameters: { layout: 'fullscreen' },
  args: { play: false },
};
export default meta;

interface Args {
  /** Play the frames, looping (`trajectory.play`). */
  play: boolean;
}

const Mdtraj = (name: string) => `https://raw.githubusercontent.com/mdtraj/mdtraj/master/tests/data/${name}`;

interface Source {
  model: string;
  modelFormat: string;
  coordinates?: string;
  coordinatesFormat?: string;
}

async function show(v: PrismViewer, s: Source, play: boolean) {
  const r = s.coordinates
    ? await v.api.load.trajectory!({
        model: s.model,
        modelFormat: s.modelFormat,
        coordinates: s.coordinates,
        coordinatesFormat: s.coordinatesFormat!,
      })
    : await v.api.load.url!({ url: s.model, format: s.modelFormat });
  if (!r.ok) console.error(`[story] ${s.coordinates ?? s.model}: ${r.message}`);
  await v.api.camera.reset({ animate: false });
  if (play) await v.api.trajectory!.play({ mode: 'loop', seconds: 5 });
  openAt(v, (n) => n.kind === 'structure');
}

const trajectory = (s: Source): StoryObj<Args> => ({
  render: ({ play }) => story(() => <RealViewerStory setup={(v) => show(v, s, play)} />),
});

/** CHARMM/NAMD DCD over a PDB model: a small protein (the repository's). */
export const PdbDcd = {
  ...trajectory({
    model: ExampleFiles.trajectoryModel,
    modelFormat: 'pdb',
    coordinates: ExampleFiles.trajectoryCoordinates,
    coordinatesFormat: 'dcd',
  }),
  name: 'PDB + DCD',
};

/** Amber NetCDF over its PRMTOP topology: the same protein (the repository's). */
export const PrmtopNetcdf = {
  ...trajectory({
    model: ExampleFiles.trajectoryTopology,
    modelFormat: 'prmtop',
    coordinates: ExampleFiles.trajectoryNetcdf,
    coordinatesFormat: 'nctraj',
  }),
  name: 'PRMTOP + NetCDF',
};

/** GROMACS XTC (compressed) over a PDB model: alanine dipeptide, 501 frames (mdtraj, network). */
export const PdbXtc = {
  ...trajectory({
    model: Mdtraj('native.pdb'),
    modelFormat: 'pdb',
    coordinates: Mdtraj('frame0.xtc'),
    coordinatesFormat: 'xtc',
  }),
  name: 'PDB + XTC',
};

/** GROMACS TRR (full precision) over a PDB model: the same run (mdtraj, network). */
export const PdbTrr = {
  ...trajectory({
    model: Mdtraj('native.pdb'),
    modelFormat: 'pdb',
    coordinates: Mdtraj('frame0.trr'),
    coordinatesFormat: 'trr',
  }),
  name: 'PDB + TRR',
};

/** LAMMPS dump (unwrapped coordinates) over a PDB model: the same run (mdtraj, network). */
export const PdbLammpstrj = {
  ...trajectory({
    model: Mdtraj('native.pdb'),
    modelFormat: 'pdb',
    coordinates: Mdtraj('frame0.lammpstrj'),
    coordinatesFormat: 'lammpstrj',
  }),
  name: 'PDB + LAMMPS dump',
};

/** A multi-frame GRO: the frames are the file's models, no coordinates file (mdtraj, network). */
export const MultiFrameGro = { ...trajectory({ model: Mdtraj('frame0.gro'), modelFormat: 'gro' }), name: 'GRO frames' };
