/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The boxes that scroll under the one rail (_spec/README.md "Input", The one rail): a rows band (display rows), a lanes
// band (lanes), a section past its cap (px: a tracks section's viewport `top`, a text's anchor moved by px). One
// interface the rail's gesture, the wheel and the keys drive, so each answers alike: its extent, its window, a set
// clamped to its ends (the same floor and ceiling whatever moved it), a page, and what Esc puts back. Internal: not in
// the behaviors' entry.

import { type RangesLay, type RowsLay, type SectionLayout, textTop, type ViewLayout } from '@molstar/seqx-core/render';
import type { ViewModel } from '@molstar/seqx-viewer/models';
import { lanesBand, pageLanes, pageRows } from './rows-band.js';
import { trackOfLayer } from './two-axes.js';

/** What a scroller is: rows and lanes a band's, a section's its own. */
export type ScrollerKind = 'rows' | 'lanes' | 'section';

export interface Scroller {
  readonly kind: ScrollerKind;
  /** The track's (rows, lanes) or the section's address. */
  readonly key: string;
  /** Its extent along (units: rows, lanes, px) and the units its window shows. */
  readonly total: number;
  readonly shown: number;
  /** The units a page moves (Page Up, Page Down): what it shows less a little, at least one. */
  readonly page: number;
  /** A move is session state, never a history step (lanes, sections: the wheel's); rows' jumps are steps. */
  readonly session: boolean;
  /** Its top now, from the stored state (moves before a frame add up). */
  top(): number;
  /** Its top at `top` clamped to `[0, total − shown]`; false when it was there. */
  set(top: number): boolean;
  /** What Esc puts back: taken at a press. */
  save(): unknown;
  restore(saved: unknown): void;
}

const clamp = (v: number, lo: number, hi: number): number => Math.max(lo, Math.min(hi, v));

/** A rows band's scroller, by its track's address. */
export function rowsScroller(view: ViewModel, l: RowsLay, track: string): Scroller {
  const area = Math.max(0, l.box.height - l.pinnedHeight);
  const total = Math.max(0, l.display.count - l.display.pinned);
  const max = l.maxTop(l.rowScale, l.display.count);
  return {
    kind: 'rows',
    key: track,
    total,
    shown: area / Math.max(1e-9, l.rowScale),
    page: pageRows(l),
    session: false,
    top: () => view.band(track)?.bandTop ?? l.top,
    set(top) {
      const was = view.band(track)?.bandTop;
      const next = clamp(top, 0, Math.max(0, max));
      if (next === (was ?? l.top)) return false;
      view.setBand(track, { bandTop: next });
      return true;
    },
    save: () => view.band(track),
    restore: (s) => view.restoreBand(track, s as ReturnType<ViewModel['band']>),
  };
}

/** A lanes band's scroller, by its track's address. */
export function lanesScroller(view: ViewModel, l: RangesLay, track: string): Scroller {
  return {
    kind: 'lanes',
    key: track,
    total: l.scrollMax + l.room,
    shown: l.room,
    page: pageLanes(l),
    session: true,
    top: () => view.band(track)?.bandTop ?? l.top,
    set(top) {
      const now = view.band(track)?.bandTop ?? l.top;
      const next = clamp(top, 0, Math.max(0, l.scrollMax));
      if (next === now) return false;
      view.setBand(track, { bandTop: next });
      return true;
    },
    save: () => view.band(track),
    restore: (s) => view.restoreBand(track, s as ReturnType<ViewModel['band']>),
  };
}

/** The px a section's page keeps in view: a page is the box less it (the last lines seen stay). */
const PageOverlapPx = 40;

/**
 * A section's own scroller while its content overflows its box (its `sectionRail`): px of its scrolled content; a
 * tracks section's top is its viewport's `top`, a text's is moved by px through its anchor (its layout's pan).
 */
export function sectionScroller(view: ViewModel, s: SectionLayout): Scroller | undefined {
  const rail = s.sectionRail;
  if (!rail) return undefined;
  const section = s.section.address;
  const p = s.placement;
  const shown = rail.end - rail.start;
  const max = Math.max(0, rail.total - shown);
  const text = s.layout.name === 'text';
  // a text's from its stored anchor through its placed lines: moves before a frame add up
  const top = (): number => {
    const v = view.get(section);
    if (text) return v ? textTop(p, v.start) : (p.top ?? 0);
    return v?.top ?? p.top ?? 0;
  };
  return {
    kind: 'section',
    key: section,
    total: rail.total,
    shown,
    page: Math.max(1, Math.max(shown / 2, shown - PageOverlapPx)),
    session: true,
    top,
    set(to) {
      const now = top();
      const next = clamp(to, 0, max);
      if (Math.abs(next - now) < 1e-6) return false;
      if (text) return view.panPx(section, 0, next - now);
      view.setTop(section, next);
      return true;
    },
    save: () => view.get(section),
    restore(saved) {
      const v = saved as ReturnType<ViewModel['get']>;
      if (v) view.setViewport(section, v);
    },
  };
}

/** The scroller a rail names: its kind and key (a rows or lanes band's track, a section's address). */
export function scrollerOf(
  view: ViewModel,
  layout: ViewLayout | undefined,
  kind: ScrollerKind,
  key: string,
): Scroller | undefined {
  for (const s of layout?.sections ?? []) {
    if (kind === 'section') {
      if (s.section.address === key) return sectionScroller(view, s);
      continue;
    }
    for (const band of s.lays)
      for (const l of band ?? []) {
        if (trackOfLayer(l.key) !== key) continue;
        if (kind === 'rows' && l.type.name === 'rows') return rowsScroller(view, l.lay as RowsLay, key);
        if (kind === 'lanes' && l.type.name === 'ranges') {
          const r = l.lay as RangesLay;
          return r.scrollMax > 0 ? lanesScroller(view, r, key) : undefined;
        }
      }
  }
  return undefined;
}

/**
 * The innermost scroller under the pointer that can move: an overflowing lanes band (its layer), a rows band (its
 * layer), else the section (`section`) while its content overflows its box. The wheel's chain and the keys read it.
 */
export function scrollersUnder(
  view: ViewModel,
  layout: ViewLayout | undefined,
  at: { readonly section?: string; readonly layer?: string },
): Scroller[] {
  const out: Scroller[] = [];
  if (at.layer !== undefined) {
    const b = lanesBand(layout, at.layer);
    if (b && b.lay.scrollMax > 0) out.push(lanesScroller(view, b.lay, b.track));
    else {
      // a layer's track, or the track itself (a rail's hit names it)
      const rows =
        scrollerOf(view, layout, 'rows', trackOfLayer(at.layer)) ?? scrollerOf(view, layout, 'rows', at.layer);
      if (rows && rows.total > rows.shown) out.push(rows);
    }
  }
  if (at.section !== undefined) {
    const s = scrollerOf(view, layout, 'section', at.section);
    if (s) out.push(s);
  }
  return out;
}

/** The first of the scrollers that moves by `by` units (a page, a wheel's px as units); false when none moved. */
export function scrollChain(chain: readonly Scroller[], by: (s: Scroller) => number): Scroller | undefined {
  for (const s of chain) if (s.set(s.top() + by(s))) return s;
  return undefined;
}
