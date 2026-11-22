/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { StaticModel } from '@molstar/seqx-core/resolve';
import type { ResolvedView } from '@molstar/seqx-core/view';
import type { ItemTitle } from './accessibility.js';
import type { StoreReader } from './layer-type.js';
import { TooltipRender } from './tooltip-hooks.js';

// An export's item titles (_spec/README.md "Render", Export, Titles): `name` (the default: a figure's size as before)
// or `card`, each exported item's named box titled by its composed card as plain text: a viewer passes its providers'
// (`cardTitle`), a headless export the spec's tooltips alone (`tooltipsAt` over the item's extent). A provider call an
// item, once an export, never a frame's.

/** What the titles read of the export's options. */
export interface TitleOptions {
  /** An item's `<title>`: its `name` or its `card` as text; the view's `export.titles` when absent, else `name`. */
  readonly titles?: 'name' | 'card';
  /** With `card` titles, an item's card as text (a viewer's composed providers); absent: the spec's tooltips alone. */
  readonly cardTitle?: ItemTitle;
}

/** The item titles an export asks for; undefined: the items' names. */
export function itemTitles(
  options: TitleOptions,
  view: ResolvedView,
  model: StaticModel,
  data: StoreReader,
): ItemTitle | undefined {
  const titles = options.titles ?? (view.export?.titles as TitleOptions['titles']) ?? 'name';
  if (titles !== 'card') return undefined;
  if (options.cardTitle) return options.cardTitle;
  return view.tooltips ? TooltipRender.chunk?.tooltipTitles(view, model, data) : undefined;
}
