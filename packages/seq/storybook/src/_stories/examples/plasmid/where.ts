/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { storyId } from '../../story-id.js';

// Where each Plasmid view's story is (examples-and-stories): the dataset's figures are Examples; the index lists them
// and the Features that show the same data.

/** A story: its title (under `Seq* (next)`) and its export. */
export interface Where {
  readonly title: string;
  readonly export: string;
}

/** The Plasmid views' figures: the stories under Examples. */
export const Figures: Readonly<Record<string, Where>> = {
  ring: { title: 'Examples/Plasmid/Ring', export: 'Ring' },
  crossing: { title: 'Examples/Plasmid/Features across the origin', export: 'FeaturesAcrossTheOrigin' },
  'rotated-linear': { title: 'Examples/Plasmid/Rotated linear view', export: 'RotatedLinearView' },
  scale: { title: 'Examples/Plasmid/Ring at scale', export: 'RingAtScale' },
  'detail-strip': { title: 'Examples/Plasmid/Ring', export: 'DetailWithAStrip' },
};

/** The editing views, and the Features that act on the same data. */
export const Tries: readonly { readonly where: Where; readonly title: string; readonly summary: string }[] = [
  {
    where: { title: 'Examples/Plasmid/Editing', export: 'EditingTheRing' },
    title: 'Editing the ring',
    summary:
      'Hide a strand’s track, show the sites of some enzymes only, bookmark bla from its start to its stop; undo.',
  },
  {
    where: { title: 'Examples/Plasmid/Editing', export: 'EditingTheStart' },
    title: 'Editing the linear view’s start',
    summary: 'The rotated linear view’s start set as a patch, undone.',
  },
  {
    where: { title: 'Features/Files', export: 'PlasmidGff3' },
    title: 'The features as GFF3 (Features)',
    summary:
      'The rotated variant’s features read from a GFF3 file, bla and a primer site across the origin as one item each.',
  },
];

/** The id of a story of `where`. */
export const idOf = (where: Where): string => storyId(where.title, where.export);
