/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { P53Views } from '@molstar/seqx-examples';
import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { hoverCardStory } from './hover-card-root.js';

// The hover card and the keys panel: the viewer's own chrome, on every view with the default spec (the examples show
// it with no code of theirs). The card is drawn through a port: a host draws the same data its own way.

const meta: Meta = { title: 'Features/Hover card and keys' };
export default meta;

const view = (id: string) => P53Views.find((v) => v.id === id)!;

/**
 * The hover card: rest the pointer on a residue, a domain or a variant, and a card beside it says what is there (the
 * position, the letter, the items) in the words the live region says; it hides while a drag pans, on Esc and when the
 * pointer leaves. Arrow keys move the keyboard cursor: the card follows it.
 */
export const HoverCard: StoryObj = {
  render: () =>
    hoverCardStory({
      view: view('tracks'),
      hint: 'Rest the pointer on a residue, a domain or a variant: the card says what is there. Drag to pan (the card goes), Esc hides it, the arrow keys move the keyboard cursor and the card with it.',
    }),
};

/**
 * A host's own card: the same data (a title, rows of a label and a value, the announcer's text) and the same place
 * (`placeCard`), drawn by the host's adapter through the hover card port, ink on paper in a monospace.
 */
export const HostCard: StoryObj = {
  render: () =>
    hoverCardStory({
      view: view('tracks'),
      host: true,
      hint: 'The same card through a host’s adapter: its own element and styles, the viewer’s data and place.',
    }),
};

/**
 * The keys panel: `?` (or the toolbar's keyboard button, "Keys and gestures (?)") opens it over the view, listing what
 * the view answers to, generated from its commands and gesture bindings, grouped; Esc closes it. Beside it, the host's
 * own list from the same data (`describeKeys`).
 */
export const KeysPanel: StoryObj = {
  render: () =>
    hoverCardStory({
      view: view('tracks'),
      keys: true,
      hint: 'The keys panel is open: Esc closes it, ? or the toolbar’s keyboard button opens it again. Below, a host’s own list from the same data.',
    }),
};
