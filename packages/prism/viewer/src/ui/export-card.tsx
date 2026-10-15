/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type CanvasToolbarItem, domView, SectionList, type ShellModel, type ShellOverlay } from '@molstar/prism-ui';
import { map } from 'rxjs';

// The Export card (the core's home for the `export` slot, viewer-spec.md §3.1): image and session with their
// features, geometry, models, video with Mol*'s extensions. Its button shows once something exports.

/** The Export card and its button: there once a section is in the `export` slot. */
export function exportCard(shell: ShellModel): { overlay: ShellOverlay; toolbar: CanvasToolbarItem } {
  return {
    overlay: {
      id: 'export',
      title: 'Export',
      // Folded until opened, but the image (`open`): the card's first use.
      render: domView(({ ctx }) => <SectionList shell={ctx.shell} slot="export" focus={ctx.focus} collapsed />),
    },
    toolbar: {
      id: 'export',
      // Last in the scene's group, under Canvas, Animation and the Story (David, 2026-10-15).
      group: 'scene',
      order: 100,
      icon: 'download',
      title: 'Export…',
      overlay: 'export',
      visible: shell.sections('export').pipe(map((list) => list.length > 0)),
    },
  };
}
