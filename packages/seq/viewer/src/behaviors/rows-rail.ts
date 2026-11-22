/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The rows rail's gestures (_spec/README.md "Input", Navigators, The rows rail): a rows band's default navigator, the
// ruler rail upright, in display rows below the pinned band. A drag on its window scrolls the rows (`bandTop`); a grip
// zooms them about the other end, which stays (`rowScale`, within the rows behavior's bounds); a press outside the
// window is a click that centres the rows there (a scroll, no step); Esc puts the band back where the press found it. A drag is
// one step of the view's history when it went farther than the rows in view or zoomed (`navigator.dragStep`). A
// double-click on its window resets the rows' zoom (`resetZoom`); the vertical overview's grips are the same gesture.
// The one rail (_spec/README.md "Input", The one rail): the same gesture drives a lanes band's rail, pan only (a
// section's is the viewer's own, `section-scroll.ts`); a press beside the window is a scroll, never a history step.

import { type BandViewport, LanesRail, MaxRowScale, MinRowScale, type WindowPart } from '@molstar/seqx-core/render';
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
  /** A lanes band's rail (in lanes): pan only, session state as the wheel's (no history step, no zoom). */
  readonly lanes?: boolean;
  /** A rail itself is a window: its hit says which rail it is (`LanesRail`). */
  readonly hit?: { readonly navigator?: string };
  /** The least px a row the layout placed it with (the larger of the host's `minPx` and the source's floor). */
  readonly rowMin?: number;
  /** The host turned the rows' zoom off (`rows.zoom.enabled`): its grips are body. */
  readonly pan?: boolean;
}

const isLanes = (w: RowsRailWindow): boolean => w.lanes ?? w.hit?.navigator === LanesRail;

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
        lanes: boolean;
        band: BandViewport | undefined;
        /** The rail as the frames draw it: its total changes when the lanes or the rows do. */
        live: RowsRailWindow;
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
      lanes: isLanes(w),
      band: b,
      live: w,
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
    // the lanes or the rows grew or shrank under the press (a rail in the live frames): the press starts again where the
    // pointer is, the band's own top and the rail's new window
    if (p.live.total !== p.total) {
      const top = this.view.band(p.track)?.bandTop ?? p.start;
      const shown = p.live.end - p.live.start;
      p.total = p.live.total;
      p.area = p.live.width;
      p.start = clamp(top, 0, Math.max(0, p.total - shown));
      p.end = p.start + shown;
      p.at0 = at;
    }
    const d = at - p.at0;
    const shown = p.end - p.start;
    if (p.part === 'body' || p.lanes || p.live.pan === true || !this.config.rows.zoom.enabled) {
      this.scrollTo(p.track, p.start + d, p.total - shown);
      return true;
    }
    // the rows a grip leaves in view, within the zoom's bounds and the rows there are
    const zoom = this.config.rows.zoom;
    const least = p.area / zoom.maxPx;
    const most = p.area / Math.max(zoom.minPx, p.live.rowMin ?? 0);
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
    if (!p || p.part === 'outside' || p.lanes) return;
    const now = this.view.band(p.track);
    const far = Math.abs((now?.bandTop ?? p.start) - p.start) > p.end - p.start;
    const zoomed = now?.rowScale !== p.band?.rowScale;
    const step = this.config.navigator.dragStep;
    if (step === 'always' || (step === 'far' && (far || zoomed))) this.view.record(p.before);
  }

  /**
   * A double-click on a rows window (the rows rail's, a vertical overview's): the rows' zoom back to the band's default
   * (`rowScale` unset), the row at the band's top kept there; a jump, its own step. False when it had none set.
   */
  resetZoom(track: string): boolean {
    this._press = undefined;
    if (this.view.band(track)?.rowScale === undefined) return false;
    return (
      this.view.jump(() => {
        this.view.setBand(track, { rowScale: undefined });
        return true;
      }) === 'moved'
    );
  }

  /** Esc: the band back where the press found it. */
  cancel(): void {
    const p = this._press;
    this._press = undefined;
    if (p && p.part !== 'outside') this.view.restoreBand(p.track, p.band);
  }

  /**
   * A click on the rail (a press that never moved; `detail` its count, the platform's): outside the window it centres
   * the rows or lanes there, a scroll as the wheel's (session state: no history step, as every rail's); on the window,
   * and a double-click's second, nothing.
   */
  click(track: string, part: WindowPart, at: number, w: RowsRailWindow, detail = 1): boolean {
    this._press = undefined;
    if (detail > 1 || part !== 'outside' || !Number.isFinite(at)) return false;
    const shown = w.end - w.start;
    return this.scrollTo(track, at - shown / 2, w.total - shown);
  }

  /** The band's top at `top`, within `[0, max]` (the rows or lanes there are); false when it was there already. */
  private scrollTo(track: string, top: number, max: number): boolean {
    const was = this.view.band(track)?.bandTop;
    const next = clamp(top, 0, Math.max(0, max));
    if (next === was) return false;
    this.view.setBand(track, { bandTop: next });
    return was !== undefined || next !== 0;
  }
}

function clamp(v: number, lo: number, hi: number): number {
  return Math.max(lo, Math.min(hi, v));
}
