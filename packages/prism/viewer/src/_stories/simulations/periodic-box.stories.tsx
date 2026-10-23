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

// Frames from a periodic simulation (`trajectory.periodic`): molecules split across
// the box by the engine made whole (`whole`), moved into the box (`box`, centered on the polymers with `center`), kept
// from jumping between frames (`nojump`: play it), or as the file has them (`off`: bonds and the backbone stretched
// across the box). `images` adds the neighbouring boxes. Split chain is the repository's (no network);
// AdK in water is MDAnalysis's test data (47 681 atoms, 10 frames, 3.8 MB), fetched from GitHub.

interface Args {
  periodic: 'whole' | 'box' | 'nojump' | 'off';
  /** With `box`: the box as it is, or the polymers at its center */
  center: 'box' | 'polymers';
  /** The frame shown, from 0 */
  frame: number;
  unitCell: boolean;
  /** Play the frames, looping */
  play: boolean;
  /** The neighbouring boxes too (`structure.symmetry`, the cells -1 to 1) */
  images: boolean;
}

const meta: Meta<Args> = {
  title: 'Viewer/Simulations/Periodic box',
  parameters: { layout: 'fullscreen' },
  args: { periodic: 'box', center: 'polymers', frame: 0, unitCell: true, play: false, images: false },
  argTypes: {
    periodic: { control: 'inline-radio', options: ['whole', 'box', 'nojump', 'off'] },
    center: { control: 'inline-radio', options: ['box', 'polymers'] },
  },
};
export default meta;

const MDAnalysis = (name: string) =>
  `https://raw.githubusercontent.com/MDAnalysis/mdanalysis/develop/testsuite/MDAnalysisTests/data/${name}`;

async function show(v: PrismViewer, load: () => Promise<{ ok: boolean; message?: string }>, args: Args) {
  const r = await load();
  if (!r.ok) console.error(`[story] ${r.message}`);
  const p = await v.api.trajectory.periodic({ mode: args.periodic, center: args.center });
  if (!p.ok) console.error(`[story] ${p.message}`);
  if (args.frame) await v.api.trajectory.frame({ index: args.frame });
  if (args.unitCell) await v.api.structure.unitCell({});
  if (args.images) await v.api.structure.symmetry({ kind: 'cells', cells: 1 });
  await v.api.camera.reset({ animate: false });
  if (args.play) await v.api.trajectory.play({ mode: 'loop', seconds: 5 });
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
