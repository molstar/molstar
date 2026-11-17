/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The render codes (_spec/README.md "Render", Issues).

import { defineIssues } from '@molstar/core/repr/tree';

export const RenderIssues = defineIssues('seq.render', {
  'seq.render.view': {
    severity: 'error',
    description: 'The spec has no view at the address asked for (or no view at all): nothing is drawn',
    hint: 'Name a view of the spec, or add one under `views`',
  },
  'seq.render.color': {
    severity: 'warning',
    description: 'A layer’s color is neither a hex color (`#rrggbb`, `#rrggbbaa`) nor a token (`token:accent`)',
    hint: 'Write a hex color or a token name; the layer’s default color is used meanwhile',
  },
  'seq.render.layer.type': {
    severity: 'warning',
    description: 'No type in a layer’s fallback chain is one this target draws: the layer is not drawn',
    hint: 'Give the layer a fallback this target draws',
  },
  'seq.render.layer.letters': {
    severity: 'warning',
    description: 'A letters layer shows another system’s letters through a path: letters are drawn on their own system',
    hint: 'Show the letters in a section of their own system',
  },
  'seq.render.axis.pending': {
    severity: 'info',
    description: 'A section’s system has no known length yet (its source is not open): its tracks wait',
    hint: 'The section draws once its source gives the length',
  },
  'seq.render.source': {
    severity: 'error',
    description: 'A source could not be opened for the render: its layers stay pending',
    hint: 'Check the source’s location and the assets the render may read',
  },
  'seq.render.unsettled': {
    severity: 'warning',
    description: 'The export did not settle before its timeout: what was pending is drawn as pending',
    hint: 'Give the export a longer timeout, or check the sources that did not answer',
  },
});
