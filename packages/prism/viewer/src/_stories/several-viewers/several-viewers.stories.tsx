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
import { ExampleFiles } from '../../testing/real-viewer.js';
import { syncCameras } from '../../sources/camera.js';
import type { PrismViewer } from '../../viewer.js';
import { ViewerGrid, type ViewerCell } from './grid.js';
import { CanvasContexts } from '../../sources/canvas-presence.js';
import { onCleanup } from '@molstar/prism-ui/dom/scope';
import { fetchPoses, poseName, showPose } from './poses.js';

// Several viewers on one page (README decision 29). Each viewer has its own shadow root, plugin, query languages, kept
// records and stored state (under its `persist.key`; two live viewers with one key are refused), and takes keys only
// while it has the focus.

const meta: Meta = { title: 'Viewer/Several viewers', parameters: { layout: 'fullscreen' } };
export default meta;

async function load(viewer: PrismViewer, url: string) {
  const r = await viewer.api.load.url!({ url, format: 'mmcif' });
  if (!r.ok) throw new Error(r.message);
}

/**
 * Three independent viewers, each its own structure and theme (1CBS, crambin, gramicidin). What to try: ⌘K opens the
 * command menu of the viewer with the focus only; selecting, hovering, the context menu, the panels, Expand and Fullscreen
 * are each viewer's own; the theme and the layout you change in one are kept under its key (`persist.key`, in
 * localStorage) and are not seen by the others; a command registered in one is not in the others.
 */
export const Independent: StoryObj = {
  render: () =>
    story(() => (
      <ViewerGrid
        keep
        minWidth={460}
        cells={[
          { key: 'several-cbs', title: '1CBS', theme: 'light', setup: (v) => load(v, ExampleFiles.cbs) },
          { key: 'several-crambin', title: 'Crambin', theme: 'dark', setup: (v) => load(v, ExampleFiles.crambin) },
          { key: 'several-grm', title: 'Gramicidin', theme: 'blue', setup: (v) => load(v, ExampleFiles.grm) },
        ]}
      />
    )),
};

/** Four poses of a ligand docked into one receptor, a viewer for each. */
const PoseCount = 4;

/**
 * Docking poses side by side: the same receptor with each of four poses of the ligand (`data/examples/docking`, the poses
 * split from one SDF), a viewer each, their cameras linked (`viewer.camera`, `syncCameras`): moving one moves the
 * others, with no echo. "Link cameras" turns it off, and the viewers move apart; on again, they come to the one moved
 * last.
 */
export const PosesWithSyncedCamera: StoryObj = {
  name: 'Poses with a synced camera',
  render: () => story(() => <Poses />),
};

function Poses() {
  let files: ReturnType<typeof fetchPoses> | undefined;
  const cells: ViewerCell[] = Array.from({ length: PoseCount }, (_, i) => ({
    key: `pose-${i + 1}`,
    title: `Pose ${i + 1}`,
    setup: async (viewer) => {
      files ??= fetchPoses(ExampleFiles.dockingReceptor, ExampleFiles.dockingPoses);
      const poses = await files;
      await showPose(viewer, poses, i);
      const title = document.querySelector(`[data-viewer="pose-${i + 1}"] strong`);
      if (title) title.textContent = poseName(poses.poses[i], i);
    },
  }));

  // The viewers are made after this runs: the controls follow what `ready` hands over.
  let linked = true;
  let viewers: readonly PrismViewer[] = [];
  let lastMoved = 0;
  let unlink: (() => void) | undefined;
  const link = () => {
    unlink?.();
    unlink = undefined;
    if (!linked || !viewers.length) return;
    unlink = syncCameras(viewers, viewers[lastMoved]);
  };

  const button = document.createElement('button');
  button.type = 'button';
  const label = () => {
    button.textContent = linked ? 'Cameras linked' : 'Cameras unlinked';
    button.setAttribute('aria-pressed', String(linked));
  };
  label();
  button.addEventListener('click', () => {
    linked = !linked;
    label();
    link();
  });
  const reset = document.createElement('button');
  reset.type = 'button';
  reset.textContent = 'Reset all';
  reset.addEventListener('click', () => {
    // Each resets its own camera; linked, the others follow the first, so it is the first that frames the scene.
    for (const v of linked ? viewers.slice(0, 1) : viewers) void v.api.camera.reset();
  });
  const bar = document.createElement('div');
  bar.style.cssText = 'display: flex; gap: 8px; padding: 2px 4px; font: 13px system-ui, sans-serif';
  bar.append(button, reset);

  return (
    <ViewerGrid
      minWidth={520}
      header={bar}
      cells={cells}
      ready={(made) => {
        viewers = made;
        // The pose that moved last is the one the viewers return to when the link is made again.
        const subs = made.map((v, i) => v.camera.changed.subscribe(() => void (lastMoved = i)));
        // Looking at the first ligand; every pose is in the receptor's frame, so the one view fits them all.
        void (async () => {
          await made[0].api.camera.frame({ target: 'ligand', animate: false });
          link();
        })();
        return () => {
          for (const s of subs) s.unsubscribe();
          unlink?.();
        };
      }}
    />
  );
}

/**
 * Fifty viewers on one page, more than the browser keeps WebGL contexts for (Chrome 16). Only the viewers in view run
 * their canvas; past the page's budget, the contexts of viewers out of view are released, the longest unseen first, and
 * a viewer scrolled back into view gets its context back (canvas-presence.ts). Scroll the grid: the viewers in view
 * always draw, none stays black. The header counts the viewers, those in view and the live contexts.
 */
export const Fifty: StoryObj = {
  name: 'Fifty viewers',
  render: () =>
    story(() => {
      const files = [ExampleFiles.crambin, ExampleFiles.grm, ExampleFiles.cbs];
      const names = ['Crambin', 'Gramicidin', '1CBS'];
      const cells: ViewerCell[] = Array.from({ length: 50 }, (_, i) => ({
        key: `fifty-${i}`,
        title: `${i + 1}. ${names[i % 3]}`,
        setup: (v) => load(v, files[i % 3]),
      }));
      const readout = document.createElement('code');
      const update = () => {
        const s = CanvasContexts.stats();
        readout.textContent = `viewers ${s.viewers} · in view ${s.inView} · live WebGL contexts ${s.live} (budget ${CanvasContexts.budget}) · in view without one ${s.dark}`;
      };
      const timer = setInterval(update, 500);
      onCleanup(() => clearInterval(timer));
      update();
      return <ViewerGrid minWidth={260} cells={cells} header={readout} />;
    }),
};
