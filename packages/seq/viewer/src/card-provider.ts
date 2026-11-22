/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { RangesList } from '@molstar/seqx-core/data';
import type { Location } from '@molstar/seqx-core/query';
import type { ItemDescription, LayerLayout, SectionLayout } from '@molstar/seqx-core/render';
import type { StaticModel } from '@molstar/seqx-core/resolve';
import type { ResolvedView } from '@molstar/seqx-core/view';
import type { StoreModel } from '@molstar/seqx-viewer/models';
import type { Place } from './describe.js';
import type { HoverCardRow } from './hover-card.js';

// The hover card composed by providers (_spec/README.md "Hover card", Providers): registry entries of kind
// `card-provider` in the viewer spec, beside `header-action`, ordered as header entries are (by `order`, then name).
// Each is given the place (the hit, its section and layer, the floored positions) and the parts so far, and adds its
// own: rows of the card, phrases of the live region's sentence. The spec says which words attach to which data; the
// application says which providers speak and in which order. Types only: `describe.ts` runs them.

/** What a description collects (`out`): the card's title, rows and the sentence's phrases. */
export interface CardParts {
  /** The card's title: a provider sets it when it is empty (the first phrase when none does). */
  title: string;
  readonly rows: HoverCardRow[];
  /** The live region's phrases, joined by `, `: never one a row (a provider writes several rows and one phrase). */
  readonly words: string[];
  /**
   * Set by a provider that says all there is (an item, a region, a pair said whole): the place's providers then say
   * nothing; a provider about the data under the mark (the spec's tooltips) reads it as it chooses.
   */
  whole: boolean;
  /** Adds a row; `item`: a row naming an item at the place (what `hover.cardItems` counts). */
  row(label: string, value: string, item?: boolean): void;
  /** Matches past the cap not made into rows (counted by an index): the card's `+k more` adds them. */
  moreCount: number;
}

/** What a provider reads beside the place: the viewer's lookups, made once a frame. */
export interface CardReads {
  /** A system's name as the spec gives it. */
  systemName(system: string): string;
  /** A layer's items as drawn in the last frame. */
  itemsOf(l: LayerLayout): readonly ItemDescription[];
  /** The items of a section's tracks under a position (each with its extent in words). */
  itemsAt(section: SectionLayout, position: number, only?: LayerLayout): { label: string; extent: string }[];
  /** The letter a letters layer of the section shows at a position ('' while none is loaded). */
  letterAt(section: SectionLayout, position: number): string;
  /** A system's letter at a position ('' while none is loaded). */
  letterOn(system: string, position: number): string;
  /** Each value track's value at a position, by its track's label. */
  valuesAt(section: SectionLayout, position: number): [string, string][];
  /** A grid's value at a cell, NaN when none. */
  gridValue(l: LayerLayout, system: string, i: number, ySystem: string, j: number): number;
  /** The position's labels in the system's declared numberings, with their ids. */
  numberedLabels(section: SectionLayout, system: string, position: number): [string, string][];
  /** An item's strand (`+`, `−`, ''). */
  strandOf(l: LayerLayout, key: string | number, start: number, row?: { t: RangesList; i: number } | null): string;
  /** An item's row in its source's loaded list tile. */
  sourceRowOf(l: LayerLayout, key: string | number, start: number): { t: RangesList; i: number } | undefined;
  /** An item's category ('' when none). */
  categoryOf(l: LayerLayout, key: string | number): string;
  /** An address as positions in the spec shown; undefined until it resolves (its data not arrived, or nothing). */
  locate(address: string): Location | undefined;
}

/** What a provider is given (`at`): the place, its section and layer in the last frame, the models. */
export interface CardAt {
  readonly words: CardReads;
  readonly model: StaticModel;
  /** The view resolved (the spec's tooltips are its), and the store its data is read from. */
  readonly view: ResolvedView | undefined;
  readonly store: StoreModel;
  readonly section: SectionLayout;
  readonly layer: LayerLayout | undefined;
  readonly place: Place;
  /** The position, floored; NaN when none. */
  readonly position: number;
  readonly yPosition: number;
  /** The `hover.cardItems` cap: a provider listing matches stops past it and counts the rest (`moreCount`). */
  readonly cardItems: number;
}

/**
 * A card provider (`{ kind: 'card-provider', name, order, describe }`): namespaced (`seq:item`, `lab:assay`); a later
 * entry under a name replaces the earlier one. `describe` adds to `out` and changes nothing else: never the hit, never
 * the emphasis. A provider that throws has its rows dropped, said once (`seq.viewer.card.provider`).
 */
export interface CardProviderEntry {
  readonly kind: 'card-provider';
  readonly name: string;
  /** Providers run by order, then name (the header entries' rule); absent: after every ordered one. */
  readonly order?: number;
  describe(at: CardAt, out: CardParts): void;
}

/** The providers a registry lists, composed: by name (a later one replacing), then sorted by order, then name. */
export function composeCardProviders(entries: readonly CardProviderEntry[]): readonly CardProviderEntry[] {
  const byName = new Map<string, CardProviderEntry>();
  for (const e of entries) {
    byName.delete(e.name);
    byName.set(e.name, e);
  }
  return [...byName.values()].sort(
    (a, b) =>
      (a.order ?? Number.POSITIVE_INFINITY) - (b.order ?? Number.POSITIVE_INFINITY) ||
      (a.name < b.name ? -1 : a.name > b.name ? 1 : 0),
  );
}
