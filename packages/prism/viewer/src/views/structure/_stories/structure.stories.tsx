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
import { MeshFiles } from '../../../testing/meshes.js';
import { ExampleFiles, openAt, RealViewerStory } from '../../../testing/real-viewer.js';
import type { PrismViewer } from '../../../viewer.js';

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

/**
 * A torus (OBJ) around 1CBS and an octahedron (PLY, colored corners) beside it: made in code, a mesh file each (the
 * Shapes group: a row per shape, its representations under it).
 */
async function loadMeshes(v: PrismViewer) {
  const s = v.plugin.managers.structure.hierarchy.current.structures[0]?.cell.obj?.data;
  const c = s?.boundary.sphere.center ?? [0, 0, 0];
  const r = s?.boundary.sphere.radius ?? 20;
  await v.api.load.files({
    files: [
      MeshFiles.torus({ center: [c[0], c[1], c[2]], radius: r + 4, tube: 2.5 }),
      MeshFiles.octahedron({ center: [c[0] + r + 16, c[1], c[2]], size: 7 }),
    ],
  });
}

/** Meshes beside a structure, the first shape picked: its look (sliders), its representations. */
export const Shapes: StoryObj = {
  render: () =>
    story(() => (
      <RealViewerStory
        load={[ExampleFiles.cbs]}
        setup={async (v) => {
          await loadMeshes(v);
          openAt(v, (n) => n.kind === 'shape');
        }}
      />
    )),
};

/** A shape's representation picked (the torus's mesh): its whole form (opacity, coloring, material, shading). */
export const ShapeRepresentation: StoryObj = {
  render: () =>
    story(() => (
      <RealViewerStory
        load={[ExampleFiles.cbs]}
        setup={async (v) => {
          await loadMeshes(v);
          openAt(v, (n) => n.kind === 'shape-representation');
        }}
      />
    )),
};
