/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

/**
 * The view files that still use rxjs (its imports and `.pipe(` calls, counted) when the guard came in (2026-10-17;
 * views-without-rxjs.test.ts): an entry only shrinks, and leaves at 0.
 */
export const RxjsInViews: Record<string, number> = {
  'packages/prism/viewer/src/views/canvas/mode-switch.tsx': 3,
  'packages/prism/viewer/src/views/common/quick-controls.tsx': 6,
  'packages/prism/viewer/src/views/export/image-export.tsx': 4,
  'packages/prism/viewer/src/views/selection/selection-summary.tsx': 3,
  'packages/prism/viewer/src/views/selection/style-layers.tsx': 2,
  'packages/prism/viewer/src/views/selection/superposed.tsx': 3,
  'packages/prism/viewer/src/views/sequence/sequence-strip.tsx': 1,
  'packages/prism/viewer/src/views/state/audio-bar.tsx': 3,
  'packages/prism/viewer/src/views/state/snapshots-panel.tsx': 7,
  'packages/prism/viewer/src/views/structure/measurements-add.tsx': 6,
  'packages/prism/viewer/src/views/structure/representation-quick.tsx': 2,
  'packages/prism/viewer/src/views/structure/scene-row.tsx': 3,
  'packages/prism/viewer/src/views/structure/scene-tree.tsx': 6,
  'packages/prism/viewer/src/views/structure/structure-tab.tsx': 8,
  'packages/prism/viewer/src/views/structure/superpose.tsx': 2,
  'packages/prism/viewer/src/views/volume/volume-levels.tsx': 22,
  'packages/prism/viewer/src/views/volume/volume-transfer.tsx': 20,
};
