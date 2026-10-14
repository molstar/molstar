/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The controls the Seq* stories that show a viewer share: its navigator, the alignment rows' color theme and the UI
// (React parts or the `<seq-viewer>` element). A story file spreads `viewerControls` into its meta.

import { NAVIGATORS, type NavigatorChoice, ROWS_COLORS, type RowsColor, type UiChoice } from './seqstar-viewer.js';

export interface Args {
  /** The navigator, switchable in the Controls panel to compare them on each case. */
  navigator: NavigatorChoice;
  /** The alignment rows' color theme, to compare the themes on the families (cases 3, 9, 10). */
  rowsColor: RowsColor;
  /** The default UI by React (`/react`) or by the `<seq-viewer>` element (`/dom`): the same parts and skin. */
  ui: UiChoice;
}

export const viewerControls = {
  // The stories' source as written, not the rendered JSX with whole documents as props, sent over the dev server's channel.
  parameters: { layout: 'fullscreen', docs: { source: { type: 'code' } } },
  args: { navigator: 'rulerRail', rowsColor: 'document', ui: 'React' } satisfies Args,
  argTypes: {
    ui: {
      control: 'inline-radio',
      options: ['React', 'DOM'],
      description: 'The default UI: React parts (`@molstar/seq-viewer/react`) or the `<seq-viewer>` element (`/dom`)',
    },
    rowsColor: {
      control: 'select',
      options: Object.keys(ROWS_COLORS),
      description: "The rows layers' color theme: the document's, or one of the themes set on every rows layer",
    },
    navigator: {
      control: 'inline-radio',
      options: Object.keys(NAVIGATORS),
      description:
        "The linear presets' navigator: rulerRail (their default), rail, strip, band, band with labels, none",
    },
  },
} as const;
