/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { story } from '@molstar/prism-ui/testing/index';
import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { RealViewerStory } from '../../testing/real-viewer.js';
import { DpdbFixtures } from './data.js';
import { DpdbCard, dynamicPdbExtension } from './extension.js';

// Dynamic PDB (dynamicpdb.com) as an extension of the viewer, written as a host would write one (extension.ts): an area
// that loads an entry's models, a card with the entry and the models' metrics, the sequence's tracks from the API's
// per-residue data and qFit's conformers, and a density row in the Structure panel (the volume server). The data is a
// copy of the API's responses (`fixtures/`): its API answers only its own origin. The deposited model, and the density,
// come from RCSB and PDBe (the network).

const meta: Meta = { title: 'Viewer/Dynamic PDB', parameters: { layout: 'fullscreen' } };
export default meta;

interface Args {
  card: boolean;
  /** The side panel pinned open: the canvas gives it room, so the card stays beside the structure. */
  pinned: boolean;
  density: boolean;
}

/** Entry dpdb_gijqn2no: FAAH (3K7F) as deposited and as qFit's multiconformer model. */
export const Entry: StoryObj<Args> = {
  args: { card: true, pinned: true, density: false },
  render: (args) =>
    story(() =>
      RealViewerStory({
        spec: (base) => ({ ...base, extensions: [...(base.extensions ?? []), dynamicPdbExtension(DpdbFixtures)] }),
        setup: async (viewer) => {
          if (args.pinned) viewer.views.shell?.setPinned('right', true);
          const person = viewer.api.as({ kind: 'person' });
          const r = await person.run('dynamic-pdb.load', { entry: 'dpdb_gijqn2no' });
          if (!r.ok) throw new Error(r.message);
          // Without the waters: the protein, its ligand and the density read better (the deposited model has 1,100).
          const models = person['dynamic-pdb']?.loaded.value?.models ?? [];
          for (const m of models)
            if (m.ref) await person.run('component.remove', { component: 'Water', structure: m.ref });
          if (args.card) await person.view.open({ card: DpdbCard });
          const deposited = models[0]?.ref;
          if (args.density && deposited) await person.run('volume.streaming', { structure: deposited, on: true });
        },
      }),
    ),
};

/** With the deposited model's density streamed from the volume server (2Fo-Fc and Fo-Fc). */
export const WithDensity: StoryObj<Args> = { ...Entry, name: 'With density', args: { ...Entry.args, density: true } };
