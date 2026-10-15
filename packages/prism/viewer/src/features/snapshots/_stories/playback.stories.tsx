/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PluginStateSnapshotManager } from '@molstar/plugin/state/manager/snapshots';
import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { story } from '@molstar/prism-ui/testing/index';
import { ExampleFiles, RealViewerStory } from '../../../testing/real-viewer.js';
import { PlaybackCard } from '../model.js';

// Playback in the real viewer: the canvas's playback card (the snapshots, the current one's transition, a trajectory's
// models) and the Story card (a snapshot's description). The MolViewSpec stories (Viewer › MVS stories) are stories
// as they are published.

const meta: Meta = { title: 'Viewer/Playback', parameters: { layout: 'fullscreen' } };
export default meta;

/**
 * A trajectory's models, three snapshots, the last with a transition: the playback card has a row each (the snapshots,
 * the transition's scrubber, the models).
 */
export const Controls: StoryObj = {
  render: () =>
    story(() => (
      <RealViewerStory
        setup={async (v) => {
          await v.api.load.trajectory?.({
            model: ExampleFiles.trajectoryModel,
            modelFormat: 'pdb',
            modelBinary: false,
            coordinates: ExampleFiles.trajectoryCoordinates,
            coordinatesFormat: 'dcd',
            coordinatesBinary: true,
          });
          await v.api.snapshot?.add({ name: 'The whole protein' });
          await v.api.snapshot?.add({ name: 'Again, for the playback' });
          // A transition of four frames over the same state (an MVS animation makes these).
          const base = v.plugin.state.getSnapshot({ camera: false, canvas3d: false });
          const frames = [0, 1, 2, 3].map(() => ({ durationInMs: 500, data: base.data! }));
          const entry = PluginStateSnapshotManager.Entry({ ...base, transition: { frames } }, { name: 'In motion' });
          v.plugin.managers.snapshot.add(entry);
          await v.api.snapshot?.apply({ entry: entry.snapshot.id });
        }}
      />
    )),
};

/** A short story on 1CBS: the Story card open on its first snapshot, its description with a link to the next. */
export const Story: StoryObj = {
  render: () =>
    story(() => (
      <RealViewerStory
        load={[ExampleFiles.cbs]}
        setup={async (v) => {
          // A snapshot takes the camera as it is: each waits for the camera to settle (the scene's bounds, the focus's
          // animation).
          const frame = () => new Promise((r) => requestAnimationFrame(r));
          const settled = async () => {
            const canvas = v.plugin.canvas3d;
            while (!canvas || canvas.boundingSphere.radius <= 0 || canvas.camera.transition.inTransition) await frame();
          };
          await settled();
          await v.api.camera.reset({ animate: false });
          await v.api.snapshot?.add({
            name: 'CRABP II',
            key: 'overview',
            description:
              'Cellular retinoic acid-binding protein II holds **retinoic acid** inside its β-barrel.\n\nSee [the ligand](#ligand).',
          });
          await v.api.focus.set({ target: 'ligand' });
          await frame();
          await settled();
          await v.api.snapshot?.add({
            name: 'The ligand',
            key: 'ligand',
            description:
              '### Retinoic acid\n- Buried in the barrel\n- Its carboxylate toward Arg111\n\nBack to [the overview](#overview).',
          });
          await v.api.snapshot?.apply({ entry: 'overview' });
          await v.api.view.open({ card: PlaybackCard });
        }}
      />
    )),
};
