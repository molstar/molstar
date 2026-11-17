/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The viewer's codes (_spec/README.md "Bindings").

import { defineIssues } from '@molstar/core/repr/tree';

export const ViewerIssues = defineIssues('seq.viewer', {
  'seq.viewer.binding.conflict': {
    severity: 'error',
    description: 'Two bindings of one layer take the same gesture, modifiers and target; the later one is used',
    hint: 'Keep one binding for each gesture, modifiers and target',
  },
  'seq.viewer.binding.action': {
    severity: 'error',
    description:
      'A binding names an action the viewer does not have (no model method, no command defined); it is left out',
    hint: 'Name a method (hover, pan, zoom) or the id of a command the viewer defines',
  },
  'seq.viewer.behavior.missing': {
    severity: 'warning',
    description: 'A gesture needs a behavior the viewer spec does not list; the gesture does nothing',
    hint: 'List the behavior in the viewer spec (the default spec has it)',
  },
  'seq.viewer.binding.key': {
    severity: 'warning',
    description: 'A key gesture in the pointer bindings: keys are bound through the keymap',
    hint: 'Give the command its key in the keymap',
  },
});
