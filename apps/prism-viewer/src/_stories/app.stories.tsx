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
import { RealViewerStory } from '@molstar/prism-viewer/testing/real-viewer';
import { FullPrismViewerSpec } from '../full-spec.js';

// The viewer as the app composes it (the viewer's defaults and Mol*'s extensions: `FullPrismViewerSpec`), for what
// needs an extension. These reach the network.

const meta: Meta = { title: 'App/Viewer', parameters: { layout: 'fullscreen' } };
export default meta;

/**
 * A predicted structure with ligands: AlphaFill's model of CYP3A4 (AlphaFold's AF-P08684-F1 with ligands transplanted
 * from PDB homologs), the cartoon colored by pLDDT (the model-archive extension's annotation), the sequence's pLDDT track
 * in the same bands, linked to the 3D view. Network.
 */
export const Predicted: StoryObj = {
  render: () =>
    story(() => (
      <RealViewerStory
        spec={() => FullPrismViewerSpec()}
        setup={async (v) => {
          const loaded = await v.api.load.url!({ url: 'https://alphafill.eu/v1/aff/P08684', format: 'mmcif' });
          if (!loaded.ok) throw new Error(loaded.message);
          await v.api.annotations.apply({ annotation: 'plddt' });
        }}
      />
    )),
};
