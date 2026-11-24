/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type ReadonlyState, State } from '@molstar/core/reactive';
import type { BandViewport } from '@molstar/seqx-core/render';
import type { Axis } from '@molstar/seqx-core/view';

// The view's bands (_spec/README.md "ViewModel", Bands): every band's own viewport by its track's address (a rows
// band's scroll, row zoom, folds and opened tree; session state), one map replaced on each change; a section's bands
// as a frame reads them, kept while the map is the same; what names a band's rows (the viewer's rows layers).

/**
 * Every band's own viewport by its track's address (_spec/README.md "ViewModel", Bands): a rows band's scroll
 * (`bandTop`, display rows below its pinned band, fractional), row zoom (`rowScale`, px a row) and folded clades
 * (`folds`, session state). One object, replaced on each change.
 */
export type Bands = ReadonlyMap<string, BandViewport>;

export const NoBands: Bands = new Map();

/**
 * What names a band's rows (the rows layer's last layout, given by the viewer): the key of the row at a band's top
 * and the fraction into it, and the top (display rows) a key stands at under a band's folds and row zoom, clamped (a
 * folded key at its summary row, a hidden one at the nearest row shown).
 */
export interface BandRows {
  keyAt(track: string, band: BandViewport): { readonly key: string; readonly offset: number } | undefined;
  topOf(track: string, band: BandViewport, key: string, offset: number): number | undefined;
  /** A rows system as an axis in its rows layer's display order (what Shift+click extends along). */
  axisOf(system: string): Axis | undefined;
}

export class ViewBands {
  private readonly _bands = new State<Bands>(NoBands);
  readonly bands: ReadonlyState<Bands> = this._bands.asReadonly;
  /** What names the bands' rows for the history (the viewer's rows layers); none: a band kept as its `bandTop`. */
  rows: BandRows | undefined;
  private readonly _kept = new Map<string, [Bands, Readonly<Record<string, BandViewport>> | undefined]>();

  /** A band's own viewport (its track's address), undefined when it has none yet. */
  get(track: string): BandViewport | undefined {
    return this._bands.value.get(track);
  }

  /** Every band now. */
  get value(): Bands {
    return this._bands.value;
  }

  /** Replaces every band (a step of the history restored); none: `NoBands`. */
  replace(map: Bands): void {
    this._bands.set(map.size ? map : NoBands);
  }

  /** Puts a band back as it was (a gesture cancelled): the viewport it held, or none at all when it had none. */
  restore(track: string, band: BandViewport | undefined): void {
    const prev = this._bands.value.get(track);
    if (prev === band || (prev && band && JSON.stringify(prev) === JSON.stringify(band))) return;
    const map = new Map(this._bands.value);
    if (band) map.set(track, band);
    else map.delete(track);
    this.replace(map);
  }

  /** Sets a band's own viewport, merged over what it holds (`ViewModel.setBand`). */
  set(track: string, band: Partial<BandViewport>): void {
    const prev = this._bands.value.get(track);
    const next: { bandTop: number; rowScale?: number; folds?: readonly string[]; tree?: boolean } = {
      bandTop: Math.max(0, band.bandTop ?? prev?.bandTop ?? 0),
    };
    const scale = 'rowScale' in band ? band.rowScale : prev?.rowScale;
    const folds = 'folds' in band ? band.folds : prev?.folds;
    // the tree a narrow frame collapsed, opened by the reader (`rows.tree`)
    const tree = 'tree' in band ? band.tree : prev?.tree;
    if (scale !== undefined) next.rowScale = scale;
    if (folds !== undefined) next.folds = folds;
    if (tree) next.tree = true;
    if (
      prev &&
      prev.bandTop === next.bandTop &&
      prev.rowScale === next.rowScale &&
      prev.folds === next.folds &&
      prev.tree === next.tree
    )
      return;
    const map = new Map(this._bands.value);
    map.set(track, next);
    this._bands.set(map);
  }

  /** The bands of a section's tracks (by track address), as a frame reads them; undefined when none has one. */
  of(section: string): Readonly<Record<string, BandViewport>> | undefined {
    const all = this._bands.value;
    if (!all.size) return undefined;
    const kept = this._kept.get(section);
    if (kept && kept[0] === all) return kept[1];
    let out: Record<string, BandViewport> | undefined;
    for (const [track, b] of all)
      if (track.startsWith(`${section}/`)) {
        out ??= {};
        out[track] = b;
      }
    this._kept.set(section, [all, out]);
    return out;
  }
}
