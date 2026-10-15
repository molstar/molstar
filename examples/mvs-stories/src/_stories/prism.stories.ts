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
import { Stories } from '../stories/index.js';

// The MolViewSpec stories of this example in the Prism viewer: each built (`buildStory`) and loaded as the viewer loads a
// MolViewSpec (`load.mvs`), its snapshots in the Snapshots list and the current one's description (markdown) on the
// canvas, stepped through with the playback card. Their data comes from the network.

const meta: Meta = { title: 'Viewer/MVS stories', parameters: { layout: 'fullscreen' } };
export default meta;

function mvsStory(id: (typeof Stories)[number]['id']): StoryObj {
  return {
    render: () =>
      story(() =>
        RealViewerStory({
          setup: async (v) => {
            const entry = Stories.find((s) => s.id === id)!;
            const r = await v.api.load.mvs!({ data: JSON.stringify(entry.buildStory()) });
            if (!r.ok) throw new Error(r.message);
          },
        }),
      ),
  };
}

/** BCR-ABL: a kinase out of control. */
export const Kinase = mvsStory('kinase');
/** The TATA-binding protein and its role in transcription initiation. */
export const Tata = mvsStory('tata');
/** RCSB PDB Molecule of the Month #1. */
export const MoleculeOfTheMonth1 = mvsStory('motm1');
/** A molecular animation. */
export const Animation = mvsStory('animation-example');
/** Audio playback. */
export const Audio = mvsStory('audio-example');
/** MolQL selectors (experimental). */
export const MolQL = mvsStory('molql');
