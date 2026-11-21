/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ExampleViews, P53Views } from '@molstar/seqx-examples';
import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { structureLinkStory } from '../structure-link/structure-link-root.js';
import { trackHeadersStory } from './track-headers-root.js';

// The track headers (Design/Track headers): actions as registry entries in the header's fixed slots, states as chips,
// the menu, the host's own parts. Every story is a view of p53 (or the plasmid) with something in its headers.

const meta: Meta = { title: 'Features/Track headers' };
export default meta;

const view = (id: string) => ExampleViews.find((v) => v.id === id)!;
const headers = view('headers');

/**
 * The built-ins on p53, drawn on the canvas: point at a track's label for its hide slot and the "⋯", which opens
 * the menu with the rest (color, layers, move). The column holds one slot, so the domains track's two states ("1
 * filtered", "1 layer hidden") are marks, their words in the card. The two hidden tracks are the application's: the view
 * draws nothing for them (no line, no chip), the toolbar's eye button lists them. A click on a
 * state undoes it in one step (the story's Cmd/Ctrl+Z, Shift for redo, puts it back; the toolbar's Back is the view's own history, not edits).
 */
export const BuiltIns: StoryObj = {
  render: () =>
    trackHeadersStory({
      view: headers,
      hint: 'Point at a track’s label: its hide slot and the “⋯” menu appear (Color and the rest are in the menu). The marks on Domains are “1 filtered” and “1 layer hidden” (hover one for its words); click a mark: the state is undone, and Cmd/Ctrl+Z (Shift+Cmd/Ctrl+Z redoes) brings it back; the toolbar’s Back is for jumps of the view, not edits. Nothing in the view says two tracks are hidden: the toolbar’s eye button (“Hidden tracks (2)”) lists them, and the “⋯” menu of the section’s tracks has “Show hidden”.',
    }),
};

/**
 * The same entries as `@molstar/ui` buttons over the headers' boxes (`headerPartsDom()`): the canvas draws nothing
 * in them. The keyboard cursor (← at a track's first position, ↑ ↓ between headers) moves the focus onto the header's
 * parts; ← and → move between them, Enter runs one, Esc goes back to the canvas; Tab is never captured.
 */
export const BuiltInsAsDomParts: StoryObj = {
  render: () =>
    trackHeadersStory({
      view: headers,
      parts: 'dom',
      hint: 'The same headers with their slots and chips as DOM buttons over the same boxes. Click in the canvas, then press ← at the first position of a track (or ↑ ↓ between headers): the focus goes onto that header’s parts; ← and → move between them, Enter runs one, Esc returns to the canvas. The figure exported from either is the same.',
    }),
};

/**
 * A host adds an entry ("Copy the track’s address", a registry entry and a command of its own) and leaves a built-in
 * out (Color): its slot is on every track, the menu lists it, the command’s message is the host’s.
 */
export const HostEntryAndALeftOutBuiltIn: StoryObj = {
  render: () =>
    trackHeadersStory({
      view: headers,
      host: true,
      labelWidth: 170,
      without: ['seq:color'],
      hint: 'The host’s “Copy the track’s address” sits beside Hide; Color is not offered, in the slots or in the menu. Click the host’s slot: its command runs and the log below says what it copied.',
    }),
};

/**
 * A section whose tracks are all hidden collapses with them, and nothing takes its place: no header, no chip. The
 * application lists the hidden tracks (the toolbar's eye button, `mountHiddenTracks`) and restores one or all.
 */
export const HiddenTracksList: StoryObj = {
  render: () =>
    trackHeadersStory({
      view: view('headers-hidden'),
      hint: 'The first section has no track left to draw and shows nothing for it: hiding never adds chrome or moves the view. The toolbar’s eye button lists the hidden tracks (“Domains”, “Conservation”, then “Show all”); a row brings one back, Show all both, in one step. Cmd/Ctrl+Z hides them again.',
    }),
};

/**
 * Tracks of 14 px have no slots and no room for a chip's words: each state is a small mark (its words in the chip's
 * card), and the track's actions are in its menu (secondary click, or the menu key).
 */
export const ShortTracks: StoryObj = {
  render: () =>
    trackHeadersStory({
      view: view('headers-short'),
      hint: 'The 14 px tracks show their states as marks (“1 filtered”, “1 layer hidden”) and have no slots: secondary-click one for its menu. The tall track below has the full header. Hover a mark: its card says the state in words.',
    }),
};

/**
 * A ring track has no header line to point at; its menu is reached by its secondary click, the menu key, or a host's
 * command. The buttons are a host's way to the same menu (`openMenu`), the entries the ring's tracks take: no Reverse,
 * no Collapse (they are for linear axes).
 */
export const RingThroughItsMenu: StoryObj = {
  render: () =>
    trackHeadersStory({
      view: view('ring'),
      menus: [
        {
          label: 'Menu of the forward strand',
          target: { section: '/views/main/ring', track: '/views/main/ring/plus' },
        },
        { label: 'Menu of the sites', target: { section: '/views/main/ring', track: '/views/main/ring/sites' } },
      ],
      hint: 'The ring’s tracks have no header line. Secondary-click a ring track, or use a button: the menu lists Hide, Color, Move; no Reverse or Collapse, which a ring does not have. The detail below has headers as usual.',
    }),
};

/**
 * A host's own DOM widget placed from `viewer.headers`: a numbered badge at the left of each track's line, at the box
 * the viewer publishes; it follows scrolling, resizing, a track hidden or moved. A click opens the track's menu.
 */
export const HostWidgetFromTheGeometry: StoryObj = {
  render: () =>
    trackHeadersStory({
      view: headers,
      widgets: true,
      hint: 'The numbered buttons at the left of each track are the page’s own, not Seq*’s: the page places each from the box the viewer publishes for that track’s header, so they stay with their tracks. Click one to open that track’s menu (Hide, Color, Move up or down: choose “Move down” or “Hide track” and the buttons follow). The menu opens over the data to the right of the button.',
    }),
};

/**
 * "Color by this track" is a host's entry, not Seq*’s: it lives in Features/Structure link with the fake structure
 * viewer it needs. The same story composed here for the header’s side: the box icon in a track’s header.
 */
export const ColorByThisTrack: StoryObj = {
  render: () =>
    structureLinkStory({
      view: P53Views.find((v) => v.id === 'structure')!,
      dataset: 'p53',
      colorBy: true,
      hint: 'The box icon in a track’s header is the host’s entry (it is in the track’s menu too): click it and the fake structure viewer below is colored by that track.',
    }),
};
