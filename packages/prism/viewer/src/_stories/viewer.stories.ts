/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { DefaultPrismViewerSpec } from '../default-spec.js';
import { createPrismViewer } from '../page.js';
import type { PrismViewer } from '../viewer.js';

// The real viewer (plan 1.6): `createPrismViewer` over a real plugin, in a shadow root, as the app makes it. The mock
// stories (Shell › Viewer mode) stay beside it until each feature's real variant replaces them. Stories may use the
// network (RCSB); tests never do.

const meta: Meta = { title: 'Viewer/Real', parameters: { layout: 'fullscreen' } };
export default meta;

/** One viewer at a time: a story that renders disposes the one before (Storybook has no unmount for this renderer). */
let current: PrismViewer | undefined;

function viewerStory(setup?: (viewer: PrismViewer) => Promise<void>): HTMLElement {
  const target = document.createElement('div');
  target.style.cssText = 'position: relative; width: 100%; height: 100vh';
  void (async () => {
    await current?.dispose();
    current = undefined;
    // Wait for the element to be in the page (the shell measures it).
    while (!target.isConnected) await new Promise((r) => requestAnimationFrame(r));
    const viewer = await createPrismViewer({ target, spec: DefaultPrismViewerSpec() });
    current = viewer;
    (window as unknown as { viewer: PrismViewer }).viewer = viewer;
    await setup?.(viewer);
  })();
  return target;
}

async function loadUrl(viewer: PrismViewer, url: string, name: string) {
  const blob = await (await fetch(url)).blob();
  await viewer.api.as({ kind: 'person' }).load.files({ files: [new File([blob], name)] as never });
}

/** Nothing loaded: the canvas with the wordmark, the rail, the canvas tools hidden. */
export const Empty: StoryObj = { render: () => viewerStory() };

/** 1TQN from RCSB, opened as a dropped file is. */
export const Loaded1TQN: StoryObj = {
  render: () => viewerStory((v) => loadUrl(v, 'https://models.rcsb.org/1TQN.bcif', '1tqn.bcif')),
};

/** 1TQN on a dark background (`canvas.background`): the chrome's ink follows the canvas. */
export const DarkBackground: StoryObj = {
  render: () =>
    viewerStory(async (v) => {
      await loadUrl(v, 'https://models.rcsb.org/1TQN.bcif', '1tqn.bcif');
      await v.api.as({ kind: 'person' }).canvas.background({ color: '#15171c' });
    }),
};
