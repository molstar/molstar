/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { BandViewport, Viewport as RenderViewport } from '@molstar/seqx-core/render';
import type { StepViewport } from '@molstar/seqx-core/story';
import type { ResolvedView } from '@molstar/seqx-core/view';
import type { ViewModel } from '@molstar/seqx-viewer/models';
import type { ViewerMotion } from './viewer-motion.js';

// The viewports a frame reads (_spec/README.md "Composition"): the renderer's window and every leading section's
// viewport with its bands, as the motion's frame slot shows them (`ViewerMotion.shownViewport`), kept while equal so
// the frame's key holds and a frame makes none.

export class ViewerViewports {
  constructor(
    private readonly view: ViewModel,
    private readonly motion: ViewerMotion,
  ) {}

  /** The window the renderer draws (reused), the viewports kept while equal and how many they are, the bands merged. */
  private readonly _windowOut = { start: 0, end: 0 };
  private _kept: Readonly<Record<string, RenderViewport>> | undefined;
  private _keptCount = 0;
  private readonly _bandsMerged = new Map<
    string,
    [RenderViewport['bands'] | undefined, NonNullable<StepViewport['bands']>, RenderViewport['bands']]
  >();

  /**
   * The window the renderer draws: the first section that leads (the renderer takes one window for all), as the frame
   * slot shows it during a motion. Undefined when it is not placed. A reused record.
   */
  window(resolved: ResolvedView, eye = true): RenderViewport | undefined {
    for (const s of resolved.sections) {
      if (this.view.isFollowing(s.address)) continue;
      const w = this.motion.shownViewport(s.address, eye);
      if (!w) return undefined;
      const out = this._windowOut;
      // written only when moved: a store of a large fractional number into the record makes a box (a genome's units)
      if (out.start !== w.start) out.start = w.start;
      if (out.end !== w.end) out.end = w.end;
      return out;
    }
    return undefined;
  }

  /**
   * Every leading section's viewport (the frame's `viewports`, by address), as its layout reads it (a window; a text's
   * anchor and letter size; a top, a band's scroll), as the frame slot shows it during a motion: each section's own,
   * field by field. Computed once a frame and kept while equal, so the frame's key holds; undefined when no section is
   * placed.
   */
  all(resolved: ResolvedView, eye = true): Readonly<Record<string, RenderViewport>> | undefined {
    const kept = this._kept;
    let n = 0;
    let same = kept !== undefined;
    for (const s of resolved.sections) {
      // a grid's own y window: its entry even when x follows (the frame reads its `y`)
      const y = this.view.get(s.address)?.y;
      if (this.view.isFollowing(s.address) && !y) continue;
      const w = this.motion.shownViewport(s.address, eye);
      if (!w) continue;
      const k = kept?.[s.address];
      if (
        !k ||
        k.start !== w.start ||
        k.end !== w.end ||
        k.scale !== w.scale ||
        k.top !== (w.top || undefined) ||
        k.y !== y ||
        k.groups !== this.view.get(s.address)?.groups ||
        k.bands !== this._bandsFor(s.address, w)
      )
        same = false;
      n++;
    }
    if (!n) return (this._kept = undefined);
    if (same && this._keptCount === n) return kept;
    const out: Record<string, RenderViewport> = {};
    for (const s of resolved.sections) {
      const y = this.view.get(s.address)?.y;
      if (this.view.isFollowing(s.address) && !y) continue;
      const w = this.motion.shownViewport(s.address, eye);
      if (!w) continue;
      const bands = this._bandsFor(s.address, w);
      const groups = this.view.get(s.address)?.groups;
      out[s.address] = {
        start: w.start,
        end: w.end,
        ...(w.scale !== undefined ? { scale: w.scale } : {}),
        // `top`: the frame's during a motion, else the model's (when set)
        ...(w.top ? { top: w.top } : {}),
        // `bands`: the ViewModel's (scroll, row zoom, folds), a shown frame's `bandTop` and `rowScale` over them
        ...(bands ? { bands } : {}),
        ...(y ? { y } : {}),
        // `groups`: the section's folded track groups, the session's set (absent: each group's `collapsed`)
        ...(groups ? { groups } : {}),
      };
    }
    this._kept = out;
    this._keptCount = n;
    return out;
  }

  /**
   * A section's bands as a frame reads them (_spec/README.md "ViewModel", Bands): the ViewModel's (a rows band's scroll,
   * row zoom and folds, session state), a shown story frame's `bandTop` and `rowScale` merged over them, the session's
   * folds kept. Kept while both are the same objects, so the frame's key holds.
   */
  private _bandsFor(section: string, w: StepViewport): RenderViewport['bands'] | undefined {
    const own = this.view.bandsOf(section);
    const framed = this.motion.shown ? w.bands : undefined;
    if (!framed) return own;
    const kept = this._bandsMerged.get(section);
    if (kept && kept[0] === own && kept[1] === framed) return kept[2];
    const out: Record<string, BandViewport> = { ...own };
    for (const track of Object.keys(framed)) {
      const f = framed[track];
      const o = own?.[track];
      out[track] = {
        bandTop: f.bandTop ?? o?.bandTop ?? 0,
        ...((f.rowScale ?? o?.rowScale) !== undefined ? { rowScale: f.rowScale ?? o?.rowScale } : {}),
        ...(o?.folds !== undefined ? { folds: o.folds } : {}),
        ...(o?.tree ? { tree: true } : {}),
      };
    }
    this._bandsMerged.set(section, [own, framed, out]);
    return out;
  }
}
