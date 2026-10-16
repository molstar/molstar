/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { story } from '@molstar/prism-ui/testing/index';
import { ExampleFiles, openAt, RealViewerStory } from '../../../testing/real-viewer.js';

// The Structure panel in the real viewer: the scene tree and the inspector of what is picked (a structure, a component,
// a representation), over the repository's example files. A map's inspector is in Viewer › Volumes.

const meta: Meta = { title: 'Viewer/Structure panel', parameters: { layout: 'fullscreen' } };
export default meta;

/** Nothing picked: the tree, the inspector's hint. */
export const Tree: StoryObj = {
  render: () => story(() => <RealViewerStory load={[ExampleFiles.cbs]} setup={(v) => openAt(v)} />),
};

/** The structure picked: its contents (counts, entities), its parameters, then the extensions' sections. */
export const Structure: StoryObj = {
  render: () =>
    story(() => <RealViewerStory load={[ExampleFiles.cbs]} setup={(v) => openAt(v, (n) => n.kind === 'structure')} />),
};

/** A component picked: its actions and its representations, each to open. */
export const Component: StoryObj = {
  render: () =>
    story(() => <RealViewerStory load={[ExampleFiles.cbs]} setup={(v) => openAt(v, (n) => n.kind === 'component')} />),
};

/** A representation picked (the cartoon): its quick settings (type, color, opacity, size) above its whole form. */
export const Representation: StoryObj = {
  render: () =>
    story(() => (
      <RealViewerStory load={[ExampleFiles.cbs]} setup={(v) => openAt(v, (n) => n.kind === 'representation')} />
    )),
};

/** Several structures (1CBS, crambin, gramicidin), each its row in the tree. */
export const ManyStructures: StoryObj = {
  render: () =>
    story(() => (
      <RealViewerStory load={[ExampleFiles.cbs, ExampleFiles.crambin, ExampleFiles.grm]} setup={(v) => openAt(v)} />
    )),
};

/**
 * Density streaming turned on (1TQN from the PDB, its maps from PDBe), the ligand focused: the streaming's row picked,
 * its maps' levels in the inspector. Network.
 */
export const StreamingDensity: StoryObj = {
  render: () =>
    story(() => (
      <RealViewerStory
        setup={async (v) => {
          await v.api.load.entry?.({ source: 'pdb', id: '1TQN' });
          await v.api.volume.streaming?.({ enabled: true } as never);
          await v.api.focus.set({ target: 'ligand' });
          openAt(v, (n) => n.kind === 'extension' && n.label.value === 'Volume streaming');
        }}
      />
    )),
};
