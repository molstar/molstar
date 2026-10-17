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
  'packages/prism/ui/src/components/action-footer/action-footer.tsx': 2,
  'packages/prism/ui/src/components/button/command-button.tsx': 4,
  'packages/prism/ui/src/components/choices/choices.tsx': 2,
  'packages/prism/ui/src/components/color/color-picker.tsx': 15,
  'packages/prism/ui/src/components/params/extra-controls.tsx': 21,
  'packages/prism/ui/src/components/params/field-binding.ts': 1,
  'packages/prism/ui/src/components/params/help.tsx': 4,
  'packages/prism/ui/src/components/params/params-form.tsx': 35,
  'packages/prism/ui/src/components/params/presets.tsx': 2,
  'packages/prism/ui/src/components/select/select.tsx': 3,
  'packages/prism/ui/src/components/virtual-list/virtual-list.tsx': 5,
  'packages/prism/ui/src/views/canvas/playback.tsx': 2,
  'packages/prism/ui/src/views/command-menu/command-menu.tsx': 22,
  'packages/prism/ui/src/views/settings/settings-tab.tsx': 2,
  'packages/prism/ui/src/views/shell/canvas-chrome.tsx': 19,
  'packages/prism/ui/src/views/shell/dock-view.tsx': 13,
  'packages/prism/ui/src/views/shell/empty-canvas.tsx': 6,
  'packages/prism/ui/src/views/shell/help-panel.tsx': 10,
  'packages/prism/ui/src/views/shell/rail.ts': 7,
  'packages/prism/ui/src/views/shell/shell-view.tsx': 22,
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
