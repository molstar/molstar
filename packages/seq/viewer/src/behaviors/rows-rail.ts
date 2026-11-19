/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The rows rail's gestures (_spec/README.md "Input", Navigators, The rows rail): a rows band's default navigator, the
// ruler rail upright, in display rows below the pinned band. A drag on its window scrolls the rows (`bandTop`); a grip
// zooms them about the other end, which stays (`rowScale`, within the rows behavior's bounds); a press outside the
// window is a click that centres the rows there (a jump); Esc puts the band back where the press found it. A drag is
// one step of the view's history when it went farther than the rows in view or zoomed (`navigator.dragStep`).

import { MaxRowScale, MinRowScale, type WindowPart } from '@molstar/seqx-core/render';
import { DefaultViewerConfig, type ViewerConfig } from '@molstar/seqx-viewer/config';
import type { ViewModel, ViewSnapshot } from '@molstar/seqx-viewer/models';

/** The rows' zoom bounds (px a row over the band's base): core's, one source (the rows rail's idle rule reads them). */
export { MaxRowScale, MinRowScale };

/** A rows rail's window as the last frame drew it: display rows below the pinned band, and the rail's length in px. */
export interface RowsRailWindow {
  readonly start: number;
  readonly end: number;
  readonly total: number;
  /** The band's scrolled area, px (the rail's length). */
  readonly width: number;
}

export class RowsRailGesture {
  private _press:
    | {
        track: string;
        part: WindowPart;
        at0: number;
        start: number;
        end: number;
        total: number;
        area: number;
        band: { bandTop: number; rowScale?: number };
        before: ViewSnapshot;
      }
    | undefined;

  constructor(
    private readonly view: ViewModel,
    private readonly config: ViewerConfig = DefaultViewerConfig,
  ) {}

  get pressed(): boolean {
    return this._press !== undefined;
  }

  /** A press on a part at display row `at`: a drag of the window (`body`, `start`, `end`); outside it a click's. */
  press(track: string, part: WindowPart, at: number, w: RowsRailWindow): boolean {
    if (!(w.total > 0) || !(w.width > 0) || !Number.isFinite(at)) return false;
    const b = this.view.band(track);
    this._press = {
      track,
      part,
      at0: at,
      start: w.start,
      end: w.end,
      total: w.total,
      area: w.width,
      band: { bandTop: b?.bandTop ?? w.start, ...(b?.rowScale !== undefined ? { rowScale: b.rowScale } : {}) },
      before: this.view.snapshot(),
    };
    return true;
  }

  /**
   * The pointer at display row `at`: the body scrolls by the rows it travelled; a grip moves its end by them, the other
   * end kept, the rows zoomed to fit (a floored window maps by travel too, so it reaches every row).
   */
  dragTo(at: number): boolean {
    const p = this._press;
    if (!p || p.part === 'outside' || !Number.isFinite(at)) return false;
    const d = at - p.at0;
    const shown = p.end - p.start;
    if (p.part === 'body') {
      this.view.setBand(p.track, { bandTop: clamp(p.start + d, 0, Math.max(0, p.total - shown)) });
      return true;
    }
    // the rows a grip leaves in view, within the zoom's bounds and the rows there are
    const least = p.area / MaxRowScale;
    const most = p.area / MinRowScale;
    if (p.part === 'start') {
      const rows = clamp(shown - d, least, Math.min(most, p.end));
      this.view.setBand(p.track, { rowScale: p.area / rows, bandTop: Math.max(0, p.end - rows) });
    } else {
      const rows = clamp(shown + d, least, Math.min(most, p.total - p.start));
      this.view.setBand(p.track, { rowScale: p.area / rows, bandTop: p.start });
    }
    return true;
  }

  /** The release: a drag farther than the rows in view, or a zoom, pushes the viewports at its press. */
  release(): void {
    const p = this._press;
    this._press = undefined;
    if (!p || p.part === 'outside') return;
    const now = this.view.band(p.track);
    const far = Math.abs((now?.bandTop ?? p.start) - p.start) > p.end - p.start;
    const zoomed = now?.rowScale !== p.band.rowScale;
    const step = this.config.navigator.dragStep;
    if (step === 'always' || (step === 'far' && (far || zoomed))) this.view.record(p.before);
  }

  /** Esc: the band back where the press found it. */
  cancel(): void {
    const p = this._press;
    this._press = undefined;
    if (p && p.part !== 'outside') this.view.setBand(p.track, p.band);
  }

  /**
   * A click on the rail (a press that never moved; `detail` its count, the platform's): outside the window it centres
   * the rows there, a jump with its own step; on the window, and a double-click's second, nothing.
   */
  click(track: string, part: WindowPart, at: number, w: RowsRailWindow, detail = 1): boolean {
    this._press = undefined;
    if (detail > 1 || part !== 'outside' || !Number.isFinite(at)) return false;
    const shown = w.end - w.start;
    return (
      this.view.jump(() => {
        const before = this.view.band(track)?.bandTop;
        this.view.setBand(track, { bandTop: clamp(at - shown / 2, 0, Math.max(0, w.total - shown)) });
        return this.view.band(track)?.bandTop !== before;
      }) === 'moved'
    );
  }
}

function clamp(v: number, lo: number, hi: number): number {
  return Math.max(lo, Math.min(hi, v));
}
