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
import type { SceneNode } from '../../../model/scene.js';
import { ExampleFiles, RealViewerStory } from '../../../testing/real-viewer.js';
import type { PrismViewer } from '../../../viewer.js';

// The Structure tab in the real viewer: the scene tree and the inspector of what is picked.

const meta: Meta = { title: 'Views/Structure tab', parameters: { layout: 'fullscreen' } };
export default meta;

/** The first node the predicate takes, depth first. */
function find(nodes: readonly SceneNode[], is: (n: SceneNode) => boolean): SceneNode | undefined {
  for (const n of nodes) {
    if (is(n)) return n;
    const below = find(n.children.value, is);
    if (below) return below;
  }
  return undefined;
}

/** The tab open, a node picked (when given). */
function tab(v: PrismViewer, pick?: (n: SceneNode) => boolean) {
  v.views.shell?.show('structure');
  const node = pick && find(v.views.scene.roots.value, pick);
  if (node) v.views.scene.pick(node);
}

/** A representation picked: the inspector shows its parameters (essentials, then Advanced). */
export const Representation: StoryObj = {
  render: () =>
    story(() => (
      <RealViewerStory load={[ExampleFiles.cbs]} setup={(v) => tab(v, (n) => n.kind === 'representation')} />
    )),
};

/** The structure picked: facts, its parameters, then the extensions' sections. */
export const Structure: StoryObj = {
  render: () =>
    story(() => <RealViewerStory load={[ExampleFiles.cbs]} setup={(v) => tab(v, (n) => n.kind === 'structure')} />),
};

/** A component picked with nothing of its own to set: its actions in place. */
export const Component: StoryObj = {
  render: () =>
    story(() => <RealViewerStory load={[ExampleFiles.cbs]} setup={(v) => tab(v, (n) => n.kind === 'component')} />),
};

/** Density streaming turned on from the structure (1TQN from the PDB, its map from PDBe): its channels. Network. */
export const StreamingDensity: StoryObj = {
  render: () =>
    story(() => (
      <RealViewerStory
        setup={async (v) => {
          await v.api.load.entry?.({ source: 'pdb', id: '1TQN' });
          await v.api.volume.streaming?.({ enabled: true } as never);
          tab(v, (n) => n.kind === 'structure');
        }}
      />
    )),
};

/** Several structures. */
export const ManyStructures: StoryObj = {
  render: () =>
    story(() => (
      <RealViewerStory load={[ExampleFiles.cbs, ExampleFiles.crambin, ExampleFiles.grm]} setup={(v) => tab(v)} />
    )),
};

/** An EMDB map next to its fitted model; the map picked: its facts and the isosurface below it. Network. */
export const EmdbMap: StoryObj = {
  render: () =>
    story(() => (
      <RealViewerStory
        setup={async (v) => {
          await v.api.load.entry?.({ source: 'pdb', id: '5A1A' });
          await v.api.load.entry?.({ source: 'emdb', id: 'EMD-2984' });
          tab(v, (n) => n.kind === 'volume');
        }}
      />
    )),
};

/** Nothing picked: the tree, the inspector's hint. */
export const TreeOnly: StoryObj = {
  render: () => story(() => <RealViewerStory load={[ExampleFiles.cbs]} setup={(v) => tab(v)} />),
};
