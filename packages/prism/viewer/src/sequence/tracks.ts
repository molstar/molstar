/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Structure } from '@molstar/model/model/structure';
import { SecondaryStructureProvider as ComputedSecondaryStructure } from '@molstar/model/props/computed/secondary-structure';
import type { LegendItem, ReadonlyState } from '@molstar/prism-ui';
import type { DocNode } from '@molstar/seq-core/format';
import type { SequenceTrackMenuItem } from './behaviors.js';
import { PrismSecondaryColors } from './palette.js';

// Track providers (sequence.md §3): what tracks a structure's sequence can have. A provider says what its track is
// and when it applies; the viewer puts it in and takes it out (per structure, sequence.md §7.1), lists it in the
// structure tree (§8.1) and offers it under Add track. Plain data and functions over the structure: no Seq* module, so
// the core spec lists them (a spec without a sequence still bundles none of Seq*).

/** What a provider builds its track over: the structure, and the document's refs it may name. */
export interface SequenceTrackContext {
  /** The structure the 3D view draws (its model carries the custom properties). */
  readonly structure: Structure;
  /** The complex's ref: a computed data node names it in `options.of`. */
  readonly complex: string;
  /** The letters' color as the person chose it (a layer's `color`): the letters track's. */
  readonly letters: Readonly<Record<string, unknown>>;
}

export interface SequenceTrackProvider {
  /** `<owner>.sequence.<name>` (viewer-spec.md §4), one namespace with behaviors. */
  readonly id: string;
  /** Its name for people (the tree's row, Add track, the band's menus). */
  readonly label: string;
  /** Where its track goes among the others: lower first. */
  readonly order: number;
  /** Its track's ref in the document (a provider's refs start with its id; the viewer's own keep their short ones). */
  readonly track: string;
  /** In by default; false: only once added (Add track). Absent: true. */
  readonly shown?: boolean;
  /** Part of the built document: a loaded SeqViewSpec brings its own, so the provider is out there unless added. */
  readonly base?: boolean;
  /** The structure carries its data; absent: always. Reads only the structure and what the plugin holds. */
  applies?(structure: Structure): boolean;
  /** Told when `applies` may answer otherwise without a new structure (a property computed later). */
  readonly whenChanged?: ReadonlyState<unknown>;
  /** Its data nodes (under the document's data node) and its one track (in the section), refs starting with its id. */
  build(ctx: SequenceTrackContext): { readonly data?: readonly DocNode[]; readonly track: DocNode };
  /** Its track's menu: colorings of the 3D view by a color theme (the first is the gutter's icon). */
  readonly menu?: readonly SequenceTrackMenuItem[];
  /**
   * What it draws, in a few words: the second line of its row in the structure tree ("Letters · Zappo", as a
   * representation's "Cartoon · Chain ID"); `letters` is the letters' color scheme by name.
   */
  describe?(o: { readonly letters: string }): string;
  /** Its legend blocks, shown while its track is in. */
  readonly legend?: readonly { readonly title: string; readonly items: readonly LegendItem[] }[];
}

/** The viewer's own tracks' refs (their menus and saved edits use them). */
export const RulerTrack = 'ruler';
export const LettersTrack = 'letters';
export const SecondaryData = 'ss';
export const SecondaryTrack = `${SecondaryData}-track`;

const track = (ref: string, layers: readonly DocNode[]): DocNode => ({
  kind: 'track',
  ref,
  params: { label: '' },
  children: [...layers],
});

/** The ruler: positions along the axis (Tracks draws it; Text draws none). */
export const RulerTrackProvider: SequenceTrackProvider = {
  id: 'prism-viewer.sequence.ruler',
  label: 'Ruler',
  order: 0,
  track: RulerTrack,
  base: true,
  describe: () => 'Positions',
  build: () => ({ track: track(RulerTrack, [{ kind: 'layer', ref: 'sequence.ruler', params: { type: 'ruler' } }]) }),
};

/** The residues' letters, in the person's colors (sequence.md §10). */
export const LettersTrackProvider: SequenceTrackProvider = {
  id: 'prism-viewer.sequence.letters',
  label: 'Residues',
  order: 10,
  track: LettersTrack,
  base: true,
  describe: (o) => `Letters · ${o.letters}`,
  build: (ctx) => ({
    track: track(LettersTrack, [
      { kind: 'layer', ref: 'sequence.letters', params: { color: { ...ctx.letters }, type: 'letters' } },
    ]),
  }),
  menu: [{ id: 'color-default', label: 'Apply default coloring', icon: 'paintbrush', theme: 'default' }],
};

/**
 * The settings of the secondary structure a structure's 3D view draws (Mol*'s computed property: Automatic, Model, DSSP
 * and its params): the track's data names them, so a change of them (the person's, in the 3D view) recomputes it.
 */
export function secondaryAssignmentOf(structure: Structure): string {
  return JSON.stringify(ComputedSecondaryStructure.props(structure.root).type);
}

/** The helices, strands and turns the 3D view draws (blocks and arrows over a coil line; Text: underlines). */
export const SecondaryTrackProvider: SequenceTrackProvider = {
  id: 'prism-viewer.sequence.secondary',
  label: 'Secondary structure',
  order: 20,
  track: SecondaryTrack,
  base: true,
  describe: () => 'Helices, strands, turns',
  build: (ctx) => ({
    data: [
      {
        kind: 'computed',
        ref: SecondaryData,
        params: {
          provider: 'molstar:secondary_structure',
          options: { of: ctx.complex, assignment: secondaryAssignmentOf(ctx.structure) },
        },
      },
    ],
    track: track(SecondaryTrack, [
      {
        kind: 'layer',
        ref: 'sequence.secondary',
        params: { data: SecondaryData, color: PrismSecondaryColors, type: 'secondary' },
      },
    ]),
  }),
  menu: [
    { id: 'color-secondary', label: 'Color by secondary structure', icon: 'paintbrush', theme: 'secondary-structure' },
  ],
};

/** The viewer's own providers (`spec.sequence.tracks`). */
export const DefaultSequenceTracks: readonly SequenceTrackProvider[] = [
  RulerTrackProvider,
  LettersTrackProvider,
  SecondaryTrackProvider,
];

/** A structure's changes over its tracks and chains (sequence.md §7.1), each a change from the document's own value. */
export interface StructureSequenceEdits {
  /** By track ref: put in or taken out (a provider's, Add track and Remove), hidden or shown (the eye). */
  readonly tracks?: Readonly<Record<string, TrackEdit>>;
  /** The chain filter; absent: every chain, copies merged. */
  readonly chains?: ChainFilter;
}

export interface TrackEdit {
  readonly in?: boolean;
  readonly hidden?: boolean;
}

/**
 * What the section shows of the structure (sequence.md §3): its polymers only, the chains of one molecule each a block
 * of their own, or the chains kept (unit names: `label_asym_id`, the deposited model's: an assembly's copies are their
 * chain once).
 */
export interface ChainFilter {
  readonly polymer?: boolean;
  readonly each?: boolean;
  readonly chains?: readonly string[];
}

/** Every structure's edits, by its state ref. */
export type SequenceEdits = Readonly<Record<string, StructureSequenceEdits>>;

/** Whether a provider's track is in for a structure: the person's choice, else its default (out under a SeqViewSpec when base). */
export function trackIn(p: SequenceTrackProvider, edits: StructureSequenceEdits | undefined, svs: boolean): boolean {
  const chosen = edits?.tracks?.[p.track]?.in;
  if (chosen !== undefined) return chosen;
  return svs && p.base ? false : p.shown !== false;
}

/** The default `in` of a provider's track, for dropping an edit that only restates it. */
export const defaultIn = (p: SequenceTrackProvider, svs: boolean) => !(svs && p.base) && p.shown !== false;

/** The providers that apply to a structure, in order. */
export function applying(providers: readonly SequenceTrackProvider[], structure: Structure): SequenceTrackProvider[] {
  return providers.filter((p) => !p.applies || p.applies(structure)).sort((a, b) => a.order - b.order);
}

/** An empty filter is none: every chain, copies merged. */
export function isFiltered(f: ChainFilter | undefined): f is ChainFilter {
  return !!f && (!!f.polymer || !!f.each || f.chains !== undefined);
}

/** One structure's edits with one track's entry replaced (undefined: dropped), empty parts dropped. */
export function withTrackEdit(
  edits: StructureSequenceEdits | undefined,
  ref: string,
  entry: TrackEdit | undefined,
): StructureSequenceEdits {
  const tracks = { ...edits?.tracks };
  if (entry && (entry.in !== undefined || entry.hidden !== undefined)) tracks[ref] = entry;
  else delete tracks[ref];
  const { tracks: _, ...rest } = edits ?? {};
  return Object.keys(tracks).length ? { ...rest, tracks } : rest;
}

/** Every structure's edits with one structure's replaced; an empty one dropped. */
export function withStructureEdits(all: SequenceEdits, ref: string, edits: StructureSequenceEdits): SequenceEdits {
  const next = { ...all };
  if (edits.tracks || isFiltered(edits.chains)) next[ref] = edits;
  else delete next[ref];
  return next;
}
