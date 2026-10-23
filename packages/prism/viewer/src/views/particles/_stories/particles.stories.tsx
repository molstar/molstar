/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { story } from '@molstar/prism-ui/testing/index';
import type { PrismViewer } from '../../../viewer.js';
import { ExampleFiles, openAt, RealViewerStory } from '../../../testing/real-viewer.js';

// Particles in the real viewer, over the repository's particle files (`data/examples/particles`, real and small; their
// sources in its README): each format, the list's inspector (its facts and how it is read), targets. A trajectory's
// frames (Simularium) are under Viewer › Simulations. The particles feature is lazy: the first particle file loads it.
// No network.

const meta: Meta = { title: 'Viewer/Particles', parameters: { layout: 'fullscreen' } };
export default meta;

/** The camera to everything */
async function frameAll(v: PrismViewer) {
  await v.api.camera.reset({ animate: false });
}

async function load(v: PrismViewer, url: string) {
  const r = await v.api.load.url!({ url });
  if (!r.ok) console.error(`[story] ${url}: ${r.message}`);
  await frameAll(v);
}

/** The list picked: its inspector (count, source, orientations, targets) and how it is read (pixel size, radius). */
const pickList = (v: PrismViewer) => openAt(v, (n) => n.kind === 'particles');

/** CryoET Data Portal: 200 ribosomes, oriented points; the portal's 6.802 Å voxels set in the list's Source. */
export const CryoEtRibosomes: StoryObj = {
  render: () =>
    story(() => (
      <RealViewerStory
        setup={async (v) => {
          await load(v, ExampleFiles.cryoetRibosomes);
          await v.api.particles!.source({ params: { pixelSize: 6.802 } });
          await v.api.particles!.representation({ type: 'orientation' });
          await frameAll(v);
          pickList(v);
        }}
      />
    )),
};

/** RELION STAR: 300 HIV capsid picks with their poses (10 Å pixels from the optics block), colored by index. */
export const RelionStar: StoryObj = {
  render: () =>
    story(() => (
      <RealViewerStory
        setup={async (v) => {
          await load(v, ExampleFiles.relionSubtomo);
          await v.api.particles!.color({ theme: 'particle-index' });
          await frameAll(v);
          pickList(v);
        }}
      />
    )),
};

/** Dynamo table, two tomograms: the second only, chosen in the list's Source (`particles.source`), 10 Å pixels. */
export const DynamoOneTomogram: StoryObj = {
  render: () =>
    story(() => (
      <RealViewerStory
        setup={async (v) => {
          await load(v, ExampleFiles.dynamo);
          await v.api.particles!.source({ params: { tomos: ['2'], pixelSize: 10 } });
          await frameAll(v);
          pickList(v);
        }}
      />
    )),
};

/** Artiatomi motivelist (EM): 225 ribosomes at the demo's 13 Å pixel size. */
export const ArtiatomiMotivelist: StoryObj = {
  render: () =>
    story(() => (
      <RealViewerStory
        setup={async (v) => {
          await load(v, ExampleFiles.artiatomi);
          await v.api.particles!.source({ params: { pixelSize: 13 } });
          await frameAll(v);
          pickList(v);
        }}
      />
    )),
};

/**
 * What the particles stand for: crambin (1CRN) at each of the 200 ribosome picks (`particles.targets`), the picks'
 * spheres hidden. The picks span thousands of Å and crambin is 30 Å across: the camera zooms in on a part of the field.
 */
export const TargetsAtPicks: StoryObj = {
  render: () =>
    story(() => (
      <RealViewerStory
        setup={async (v) => {
          await load(v, ExampleFiles.cryoetRibosomes);
          await v.api.particles!.source({ params: { pixelSize: 6.802 } });
          const file = new File([await (await fetch(ExampleFiles.crambin)).blob()], '1crn.cif');
          await v.api.particles!.targets({ file });
          await v.api.particles!.visibility({ visible: false, representation: 'spacefill' });
          await frameAll(v);
          // 100×: a few picks across the view, each crambin as it sits
          await v.api.camera.zoom({ factor: 10, animate: false });
          await v.api.camera.zoom({ factor: 10, animate: false });
          pickList(v);
        }}
      />
    )),
};
