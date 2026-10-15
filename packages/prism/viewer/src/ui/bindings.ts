/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { InputBindingGroup } from '@molstar/prism-ui';

/** Mol*'s canvas bindings (the trackball and the interactivity behaviors), as the Help tab lists them. */
export const CanvasBindings: readonly InputBindingGroup[] = [
  {
    title: 'Camera',
    bindings: [
      { action: 'Rotate', input: ['Drag'], or: ['One-finger drag'] },
      { action: 'Move', input: ['Right drag'], or: ['⌃', 'Drag'] },
      { action: 'Zoom', input: ['Scroll'], or: ['Pinch'] },
      { action: 'Clip (move the slab)', input: ['⇧', 'Scroll'] },
      { action: 'Roll', input: ['⇧', 'Drag'] },
      { action: 'Reset the camera', input: ['Double-click', 'empty space'] },
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
];
