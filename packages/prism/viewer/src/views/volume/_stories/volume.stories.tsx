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
import { openAt, RealViewerStory } from '../../../testing/real-viewer.js';

// A map's inspector in the real viewer: an EMDB entry (EMD-2984) next to its fitted model (5A1A). Network.

const meta: Meta = { title: 'Viewer/Volumes', parameters: { layout: 'fullscreen' } };
export default meta;

/** The map picked: its levels (the histogram with the iso value to drag, the color, the opacity). */
export const Isosurface: StoryObj = {
  render: () =>
    story(() => (
      <RealViewerStory
        setup={async (v) => {
          await v.api.load.entry?.({ source: 'pdb', id: '5A1A' });
          await v.api.load.entry?.({ source: 'emdb', id: 'EMD-2984' });
          openAt(v, (n) => n.kind === 'volume');
        }}
      />
    )),
};

/**
 * The map as a direct volume, its representation picked: the transfer function (the opacity curve over the map's
 * histogram, the colors under it, the color stops).
 */
export const DirectVolume: StoryObj = {
  render: () =>
    story(() => (
      <RealViewerStory
        setup={async (v) => {
          await v.api.load.entry?.({ source: 'emdb', id: 'EMD-2984' });
          await v.api.volume.representation({ type: 'direct-volume' });
          openAt(v, (n) => n.kind === 'volume-representation' && /direct/i.test(n.label.value));
        }}
      />
    )),
};
