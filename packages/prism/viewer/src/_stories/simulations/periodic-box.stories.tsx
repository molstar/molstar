/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { story } from '@molstar/prism-ui/testing/index';
import { PolymerCenter } from '@molstar/plugin/state/helpers/periodic';
import { ModelFromTrajectory } from '@molstar/plugin/state/transforms/structure/hierarchy';
import type { PrismViewer } from '../../viewer.js';
import { ExampleFiles, openAt, RealViewerStory } from '../../testing/real-viewer.js';

// Frames from a periodic simulation: molecules split across the box by the engine, made whole by the trajectory's
// models (`ModelFromTrajectory`'s `periodic`, .v6/periodic-trajectories/plan.md), and moved into the box (`box`),
// centered on the polymers (`center: polymers`). `periodic: off` shows the frames as the file has them: bonds and the
// backbone stretched across the box. Split chain is the repository's (no network);
// AdK in water is MDAnalysis's test data (47 681 atoms, 10 frames, 3.8 MB), fetched from GitHub.

interface Args {
  periodic: 'whole' | 'box' | 'off';
  /** With `box`: the box as it is, or the polymers at its center */
  center: 'box' | 'polymers';
  /** The frame shown, from 0 */
  frame: number;
  unitCell: boolean;
}

const meta: Meta<Args> = {
  title: 'Viewer/Simulations/Periodic box',
  parameters: { layout: 'fullscreen' },
  args: { periodic: 'box', center: 'polymers', frame: 0, unitCell: true },
  argTypes: {
    periodic: { control: 'inline-radio', options: ['whole', 'box', 'off'] },
    center: { control: 'inline-radio', options: ['box', 'polymers'] },
  },
};
export default meta;

const MDAnalysis = (name: string) =>
  `https://raw.githubusercontent.com/MDAnalysis/mdanalysis/develop/testsuite/MDAnalysisTests/data/${name}`;

async function show(v: PrismViewer, load: () => Promise<{ ok: boolean; message?: string }>, args: Args) {
  const r = await load();
  if (!r.ok) console.error(`[story] ${r.message}`);
  const update = v.plugin.build();
  for (const m of v.plugin.managers.structure.hierarchy.current.models) {
    update.to(m.cell).update(ModelFromTrajectory, (p) => ({
      ...p,
      modelIndex: args.frame,
      periodic: args.periodic,
      periodicCenter:
        args.center === 'polymers'
          ? { name: 'selection' as const, params: { script: PolymerCenter } }
          : { name: 'box' as const, params: {} },
    }));
  }
  await update.commit();
  if (args.unitCell) await v.api.structure.unitCell({});
  await v.api.camera.reset({ animate: false });
  openAt(v, (n) => n.kind === 'structure');
}

/** A four-residue peptide and three waters in a rhombic dodecahedron: split across a face in most frames. */
export const SplitChain: StoryObj<Args> = {
  name: 'Split chain',
  render: (args) =>
    story(() => (
      <RealViewerStory setup={(v) => show(v, () => v.api.load.url!({ url: ExampleFiles.splitPeptide }), args)} />
    )),
};

/** Adenylate kinase in water (rhombic dodecahedron): the protein is split across the box in every frame as written. */
export const AdkInWater: StoryObj<Args> = {
  name: 'AdK in water',
  render: (args) =>
    story(() => (
      <RealViewerStory
        setup={(v) =>
          show(
            v,
            () =>
              v.api.load.trajectory!({
                model: MDAnalysis('adk_oplsaa.gro'),
                modelFormat: 'gro',
                coordinates: MDAnalysis('adk_oplsaa.xtc'),
                coordinatesFormat: 'xtc',
              }),
            args,
          )
        }
      />
    )),
};
