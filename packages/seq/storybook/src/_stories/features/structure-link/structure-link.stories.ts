/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ComplexViews, P53Views } from '@molstar/seqx-examples';
import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { structureLinkStory } from './structure-link-root.js';

// The structure link (Design/Structure syncing): a Seq* view beside a fake structure viewer, linked by `linkStructure`
// through the host's adapter. The fake draws each chain instance as a row of residue cells (the gaps are what the entry
// has no atoms for), has a grammar of its own and speaks MolViewSpec; what crosses is listed under it.

const meta: Meta = { title: 'Features/Structure link' };
export default meta;

const p53 = (id: string) => P53Views.find((v) => v.id === id)!;
const complex = (id: string) => ComplexViews.find((v) => v.id === id)!;

/**
 * Hover R175 on P04637: the same residue lights on every chain of the entity in the fake, and back. Hover a cell in the
 * fake: it lights on P04637 and on that chain's block, never on the sibling copies; the gaps (the ends the crystal has
 * no atoms for) light nothing either way. Hover the zinc chip: its chain's zinc lights.
 */
export const HoverBothWays: StoryObj = {
  render: () =>
    structureLinkStory({
      view: p53('structure'),
      dataset: 'p53',
      hint: 'Hover a residue of P04637: it lights on every chain of the structure. Hover a cell of a chain (or the zinc chip): it lights on P04637 and on that chain. The dashed gaps are unmodeled ends: they light nothing.',
    }),
};

/**
 * The fake answers 600 ms late, as a real viewer may: move the pointer from a residue of P04637 to a cell of a chain
 * and away; what the fake tells arrives after the pointer has moved on, and the late message never clears or replaces
 * the hover the pointer has now.
 */
export const HoverWithLatency: StoryObj = {
  render: () =>
    structureLinkStory({
      view: p53('structure'),
      dataset: 'p53',
      latency: 600,
      hint: 'The fake tells Seq* what it hovers 600 ms late. Sweep the pointer across the fake and back to P04637: the stale messages do not take the hover from the pointer.',
    }),
};

/**
 * Select in either view: a click or a drag on P04637 selects the same residues on every chain, a click on a cell of a
 * chain selects it on P04637 (revealed, the zoom kept). Cmd/Ctrl+click toggles, Esc clears. "Select chain" and "Within
 * 5 Å" make selections only the structure has: they arrive on P04637 whole.
 */
export const SelectionBothWays: StoryObj = {
  render: () =>
    structureLinkStory({
      view: p53('structure-split'),
      dataset: 'p53',
      hint: 'Click or drag on P04637 to select; click a cell in the fake to select there (Cmd/Ctrl+click toggles, Esc clears). "Select chain" and "Within 5 Å" select in the fake alone: Seq* receives them and reveals them.',
    }),
};

/**
 * Double-click a cell of the fake: its camera frame moves over 400 ms and Seq* focuses the residue. Focus a residue in
 * Seq* (double-click, or Enter on a selection): the fake's frame follows over the same time; Alt+Left goes back, and
 * the fake's frame with it.
 */
export const FocusAndBack: StoryObj = {
  render: () =>
    structureLinkStory({
      view: complex('units'),
      dataset: 'complex',
      hint: '1FIN: CDK2 and cyclin A, two copies, an ATP on each CDK2. Double-click a residue in either view to focus it; Alt+Left goes back. The dashed frame in the fake is its camera; it moves over the time the jump took.',
    }),
};

/**
 * P04637 reaches the entity by two maps: SIFTS's `structure` map (ranked first) and an alignment shifted by one residue.
 * The first section of the complex follows the ranked map, the second goes `via` the alignment's; a residue hovered on
 * P04637 lights two neighbors. The link goes through the ranked map unless the toggle says the alignment: the fake
 * lights the residue the chosen map gives.
 */
export const ChainMaps: StoryObj = {
  render: () =>
    structureLinkStory({
      view: p53('structure-compare'),
      dataset: 'p53',
      chooseMap: true,
      hint: 'Hover a residue of P04637: the two sections of the complex light neighbors (one follows SIFTS, the other the alignment, shifted by one). The fake lights what the link’s map says; the toggle changes the map the link goes through.',
    }),
};

/**
 * "Color the structure by this track": a host's entry on every track's header (the box icon). It writes the track's colors
 * as the fake's MolViewSpec state (a snapshot, not following the track): try P04637's letters (amino acid colors),
 * the domains, the secondary structure.
 */
export const ColorByThisTrack: StoryObj = {
  render: () =>
    structureLinkStory({
      view: p53('structure'),
      dataset: 'p53',
      colorBy: true,
      hint: 'Click the box icon in a track’s header (or open its menu): the structure is colored by that track, on every chain it reaches. The gaps stay gray: the structure has no atoms there.',
    }),
};
