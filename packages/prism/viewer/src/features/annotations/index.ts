/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { AnnotationsApi } from '../../api/annotations.js';
import { operationCommand } from '../../commands/operations.js';
import type { ViewerExtension } from '../../spec.js';

// The `annotations` feature (viewer-spec.md §3.1; README decision 27): what a structure carries and the viewer can
// show (B-factors, secondary structure, and the providers extensions register: pLDDT, validation), as
// `annotations.list` and `annotations.apply`. Its providers color representations: it needs `structure`.

export const AnnotationsFeature: ViewerExtension = {
  id: 'annotations',
  title: 'Annotations',
  builtIn: true,
  requires: ['structure'],
  api: { annotations: (c) => new AnnotationsApi(c) },
  commands: (api) => [
    operationCommand(api, api.annotations.apply, {
      category: 'Style',
      icon: 'tag',
      title: 'Show an annotation…',
      alias: ['annotation', 'b-factor', 'secondary structure', 'color by annotation'],
      action: 'Show',
    }),
  ],
};
