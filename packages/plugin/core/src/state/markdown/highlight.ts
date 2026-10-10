/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { MarkdownExtension } from '../manager/markdown-extensions.js';
import { findRepresentations, parseArray } from './helpers.js';

export const HighlightRefsMarkdownExtension: MarkdownExtension = {
  name: 'highlight-refs',
  execute: ({ event, args, manager }) => {
    const refs = parseArray(args['highlight-refs']);
    if (!refs?.length) return;

    if (event === 'mouse-leave' && refs.length) {
      manager.plugin.managers.interactivity.lociHighlights.clearHighlights();
      return;
    } else if (event === 'mouse-enter') {
      const cells = manager.findCells(refs);
      for (const cell of findRepresentations(manager.plugin, cells)) {
        if (!cell.obj?.data) continue;
        const { repr } = cell.obj.data;
        for (const loci of repr.getAllLoci()) {
          manager.plugin.managers.interactivity.lociHighlights.highlight({ loci, repr }, false);
        }
      }
    }
  },
};
