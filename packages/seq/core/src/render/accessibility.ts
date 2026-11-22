/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Location } from '@molstar/seqx-core/query';
import type { StaticModel } from '@molstar/seqx-core/resolve';
import type { ItemDescription } from './layer-type.js';
import type { SvgAccessible } from './svg.js';
import type { EmphasisHighlight } from './layer-type.js';
import type { LayerLayout, SectionLayout, ViewLayout } from './tracks.js';

/** An item's title in an export with `card` titles (the composed card's rows as text); undefined: its name. */
export type ItemTitle = (
  section: SectionLayout,
  layer: LayerLayout,
  item: ItemDescription,
  name: string,
) => string | undefined;

// Words for what a figure shows (_spec/README.md "Render", Accessibility): the titles of an exported SVG and the
// sentences a viewer's live region says come from the same functions, so a reader of the export and a reader of the
// viewer hear the same names. Positions are 1-based and inclusive, as every label in Seq* says them.

const count = (n: number) => Math.round(n).toLocaleString('en-US');

/** Positions `[start, end)` as `12–20` (`12` for one). */
export function spanText(start: number, end: number): string {
  return end - start <= 1 ? count(start + 1) : `${count(start + 1)}–${count(end)}`;
}

/** The ranges of a location part as `1–10, 20` (at most `max` of them, then how many more). */
export function rangesText(ranges: ArrayLike<number>, max = 3): string {
  const n = ranges.length >> 1;
  const out: string[] = [];
  for (let i = 0; i < Math.min(n, max); i++) out.push(spanText(ranges[2 * i], ranges[2 * i + 1]));
  if (n > max) out.push(`${n - max} more`);
  return out.join(', ');
}

/** How many positions a location holds on its first system, and where, as `3 positions, p53 12–14`. */
export function locationText(model: StaticModel, loc: Location | undefined): string {
  if (!loc || !loc.parts.length) return 'nothing';
  const part = loc.parts[0];
  let n = 0;
  for (let i = 0; i + 1 < part.ranges.length; i += 2) n += part.ranges[i + 1] - part.ranges[i];
  const name = model.systems.get(part.system)?.name ?? part.system;
  return `${count(n)} position${n === 1 ? '' : 's'}, ${name} ${rangesText(part.ranges)}`;
}

/**
 * An item's extent: one run (`12–20`), or across a cyclic origin its two runs from its start (`2,487–2,686, 1–175`).
 */
export function extentText(item: Pick<ItemDescription, 'start' | 'end' | 'runs'>): string {
  return item.runs ? rangesText(item.runs) : spanText(item.start, item.end);
}

/** An item as a title: its label and its extent on its system (`Domain A, p53 12–20`). */
export function itemText(system: string, item: ItemDescription): string {
  return `${item.label}, ${system} ${extentText(item)}`;
}

/**
 * The most items an export names in all (a figure of thousands names its tracks, and says how many more, per track
 * and in the description) and the most a track names (so one track does not use the figure's whole share).
 */
export const MaxNamedItems = 2000;
export const MaxNamedItemsPerTrack = 200;

/**
 * What an export's SVG says of the figure: its title (the spec's), a description of its sections and tracks, each
 * track's band and its items' places with their titles, and the highlights'. A layout of lines (a text) names its
 * tracks' first line.
 */
export function accessibleOf(
  model: StaticModel,
  layout: ViewLayout,
  title: string,
  highlights: readonly EmphasisHighlight[] = [],
  /** `card` titles: an item's composed card as text (absent: its name, the default). */
  titleOf?: ItemTitle,
): SvgAccessible {
  const tracks: SvgAccessible['tracks'][number][] = [];
  const sections: string[] = [];
  let named = 0;
  let unnamed = 0;
  for (const s of layout.sections) {
    const system = model.systems.get(s.section.axis.system)?.name ?? s.section.axis.system;
    sections.push(`${s.id} (${system}, ${s.tracks.length} track${s.tracks.length === 1 ? '' : 's'})`);
    for (const t of s.tracks) {
      const items: SvgAccessible['tracks'][number]['items'][number][] = [];
      let total = 0;
      // its own layers, then its `section` layers (pairs: laid out once a section, their arcs named on their track)
      for (const l of [...t.layers, ...s.sectionLays.filter((x) => t.track.layers.includes(x.layer))]) {
        const described = l.type.items?.(l.lay);
        if (!described) continue;
        const dataSystem = (l.layer.data?.system && model.systems.get(l.layer.data.system)?.name) || system;
        for (const it of described) {
          total++;
          if (items.length >= MaxNamedItemsPerTrack || named >= MaxNamedItems) continue;
          named++;
          const name = itemText((it.system && model.systems.get(it.system)?.name) || dataSystem, it);
          items.push({
            title: (titleOf && titleOf(s, l, it, name)) || name,
            x: it.x0,
            y: it.y ?? t.box.y,
            width: Math.max(0, it.x1 - it.x0),
            height: it.height ?? t.box.height,
          });
        }
      }
      unnamed += total - items.length;
      const more = total > items.length ? `, ${total - items.length} more not named` : '';
      tracks.push({
        title: `${t.track.label}${total ? `: ${total} item${total === 1 ? '' : 's'}${more}` : ''} (${s.id})`,
        x: t.box.x,
        y: t.box.y,
        width: t.box.width,
        height: t.box.height,
        items,
      });
    }
  }
  return {
    title,
    description:
      (sections.length
        ? `${sections.length} section${sections.length === 1 ? '' : 's'}: ${sections.join('; ')}`
        : 'An empty view') + (unnamed ? `. ${count(unnamed)} more items not named` : ''),
    tracks,
    highlights: highlights.map((h, i) => `Highlight ${i + 1}: ${locationText(model, h.location)}`),
  };
}
