/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { getCellBoundingSphere } from '../manager/focus-camera/focus-object.js';
import type { MarkdownExtension } from '../manager/markdown-extensions.js';
import { findRepresentations, parseArray } from './helpers.js';

export const CenterCameraMarkdownExtension: MarkdownExtension = {
  name: 'center-camera',
  execute: ({ event, args, manager }) => {
    if (event !== 'click') return;
    if ('center-camera' in args) {
      manager.plugin.managers.camera.reset();
    }
  },
};

export const FocusRefsMarkdownExtension: MarkdownExtension = {
  name: 'focus-refs',
  execute: ({ event, args, manager }) => {
    if (event !== 'click') return;
    const refs = parseArray(args['focus-refs']);
    if (!refs?.length) return;

    const cells = manager.findCells(refs);
    if (!cells.length) return;

    const reprs = findRepresentations(manager.plugin, cells);
    if (!reprs.length) return;

    const spheres = reprs.map((c) => getCellBoundingSphere(manager.plugin, c.transform.ref)).filter((s) => !!s);
    if (!spheres.length) return;
    manager.plugin.managers.camera.focusSpheres(spheres, (s) => s, { extraRadius: 3 });
  },
};
