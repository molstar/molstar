/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { InputBindingGroup } from '@molstar/prism-ui';

/**
 * Mol*'s canvas bindings, as the Help tab lists them (parity.md §5, the keyboard rows): the trackball's
 * (`DefaultTrackballBindings`), the interactivity behaviors', and the keys of the `camera-controls` and
 * `snapshot-controls` behaviors. Written out (the plugin's own `Binding.format` reads as "Press key w"), so
 * `_test/bindings.test.ts` checks every key here against the plugin's defaults; a host that changes the plugin's
 * bindings changes what the keys do, not this list.
 */
export const CanvasBindings: readonly InputBindingGroup[] = [
  {
    title: 'Camera',
    bindings: [
      { action: 'Rotate', input: ['Drag'], or: ['One-finger drag'] },
      { action: 'Move', input: ['Right drag'], or: ['⌃', 'Drag'] },
      { action: 'Zoom', input: ['Scroll'], or: ['Pinch'] },
      { action: 'Clip (move the slab)', input: ['⇧', 'Scroll'] },
      { action: 'Roll', input: ['⇧', '⌃', 'Drag'] },
      { action: 'Reset the camera', input: ['Right click', 'empty space'] },
    ],
  },
  {
    title: 'Picking',
    bindings: [
      { action: 'Focus what is clicked (Interact)', input: ['Click'] },
      { action: 'Select what is clicked (Select)', input: ['Click'] },
      { action: 'Add to the selection', input: ['⇧', 'Click'] },
      { action: 'Select everything between', input: ['⇧', 'Click', 'a second residue'] },
      { action: 'Clear the selection', input: ['Click', 'empty space'] },
      { action: 'Show what is under the pointer', input: ['Hover'], or: ['Long press'] },
    ],
  },
  {
    // The `camera-controls` behavior: keys while the pointer is over the canvas.
    title: 'View keys',
    bindings: [
      { action: 'Spin', input: ['I'] },
      { action: 'Rock', input: ['O'] },
      { action: 'Reset the view', input: ['T'] },
      { action: 'Global illumination on or off', input: ['G'] },
      { action: 'Fly mode on or off', input: ['⇧', 'Space'] },
    ],
  },
  {
    // The trackball's keys: moving through the scene (fly mode, or a close look).
    title: 'Moving with keys',
    bindings: [
      { action: 'Forward', input: ['W'] },
      { action: 'Back', input: ['S'] },
      { action: 'Left', input: ['A'] },
      { action: 'Right', input: ['D'] },
      { action: 'Up', input: ['R'] },
      { action: 'Down', input: ['F'] },
      { action: 'Roll left', input: ['Q'] },
      { action: 'Roll right', input: ['E'] },
      { action: 'Pitch up', input: ['⇧', '↑'] },
      { action: 'Pitch down', input: ['⇧', '↓'] },
      { action: 'Yaw left', input: ['⇧', '←'] },
      { action: 'Yaw right', input: ['⇧', '→'] },
      { action: 'Faster, while held', input: ['⇧'] },
      { action: 'Look with the mouse (pointer lock)', input: ['⌃', 'Space'] },
    ],
  },
  {
    // The `snapshot-controls` behavior.
    title: 'Snapshot keys',
    bindings: [
      { action: 'Next snapshot', input: ['⌃', '→'] },
      { action: 'Previous snapshot', input: ['⌃', '←'] },
      { action: 'First snapshot', input: ['⌃', '↑'] },
      { action: 'Last snapshot', input: ['⌃', '↓'] },
    ],
  },
];
