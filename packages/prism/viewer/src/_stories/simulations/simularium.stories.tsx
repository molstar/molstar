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

// A Simularium trajectory in the real viewer (`data/examples/particles/small.simularium`): a simulation's frames as a
// particle list, its unit cell. The particles feature is lazy: the file loads it. No network.

const meta: Meta = { title: 'Viewer/Simulations/Simularium', parameters: { layout: 'fullscreen' } };
export default meta;

async function load(v: PrismViewer, url: string) {
  const r = await v.api.load.url!({ url });
  if (!r.ok) console.error(`[story] ${url}: ${r.message}`);
  await v.api.camera.reset({ animate: false });
}

/** The list picked: its inspector, and its frames' row under it. */
const pickList = (v: PrismViewer) => openAt(v, (n) => n.kind === 'particles');

/** Its frames' row under the list (seek, step, play) and its unit cell as a row. */
export const Frames: StoryObj = {
  render: () =>
    story(() => (
      <RealViewerStory
        setup={async (v) => {
          await load(v, ExampleFiles.simularium);
          await v.api.particles!.unitCell({});
          await v.api.particles!.frame({ index: 2 });
          pickList(v);
        }}
      />
    )),
};

/** Playing it (`particles.play`), looping. */
export const Playing: StoryObj = {
  render: () =>
    story(() => (
      <RealViewerStory
        setup={async (v) => {
          await load(v, ExampleFiles.simularium);
          await v.api.particles!.play({ mode: 'loop', seconds: 3 });
          pickList(v);
        }}
      />
    )),
};
