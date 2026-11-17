/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The resolved view's codes (_spec/README.md "Resolved view", Issues).

import { defineIssues } from '@molstar/core/repr/tree';

export const ViewIssues = defineIssues('seq.view', {
  'seq.view.layer.data': {
    severity: 'error',
    description: 'A layer’s data names nothing the static model holds: no sequence, no annotation',
    hint: 'Name a sequence or an annotation of the spec, or a query that starts at one',
  },
  'seq.view.layer.unplaced': {
    severity: 'warning',
    description: 'A layer’s data is on a system no static path joins to its section’s; it is not drawn',
    hint: 'Add a map between the systems, or show the data in a section of its own system',
  },
  'seq.view.follow.unplaced': {
    severity: 'warning',
    description: 'A section follows one whose system no static path joins to its own; it does not follow',
    hint: 'Add a map between the sections’ systems',
  },
});
