/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { DefaultNavigators, type PresetEntry } from '@molstar/seq-viewer';

// The sequence model's presets (prism-sequence.md §4): the tracks version without Seq*'s header column (the `strip`
// preset's look: `headerWidth` 0), its rail the person's or the host's choice: which navigator, how tall. Registered by
// the model as a host preset, so Seq* needs no mode of its own.

/** The tracks version's rail (the zoomable axis' navigator): the ruler's own rail, or one of Seq*'s navigators. */
export interface TracksRail {
  /** `ruler`: a hairline on the ruler track (no row of its own); the others are Seq*'s navigators of that name. */
  readonly navigator: 'ruler' | 'strip' | 'band' | 'rail' | 'none';
  /** The navigator's row, px; absent: its own (the ruler's rail has none). */
  readonly height?: number;
}
export const DefaultTracksRail: TracksRail = { navigator: 'ruler' };

/** The preset the Text version uses: Seq*'s `compact` (its render props are the model's, the lines centered). */
export const TextPreset = 'compact';
/** The preset the Tracks version uses. */
export const TracksPreset = 'prism_tracks';

const NAVIGATORS = { ruler: 'rulerRail', strip: 'strip', band: 'band', rail: 'rail', none: 'none' } as const;

/** The tracks' left gutter where a coloring icon sits per track (`TrackGutter`), px: the icon's cell and the air after it. */
export const TrackGutterPx = 28;

/** `gutter`: the tracks leave `TrackGutterPx` at their left for the track icons (the host offers colorings). */
export function tracksPreset(rail: TracksRail, gutter = false): PresetEntry {
  const name = NAVIGATORS[rail.navigator];
  const own = DefaultNavigators.find((n) => n.name === name);
  const navigator =
    rail.height !== undefined && own && !own.parts.onRuler ? { ...own.parts, height: rail.height } : name;
  return {
    name: TracksPreset,
    geometry: 'linear',
    render: { headerWidth: 0, gutterLeft: 8 + (gutter ? TrackGutterPx : 0), gutterRight: 8 },
    navigator,
  };
}
