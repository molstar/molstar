/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type ReadonlyState, State } from '@molstar/core/reactive';
import { cancelFrameCallback, requestFrameCallback } from '@molstar/ui/dom/frame';
import { DomHost } from '@molstar/ui/dom/host';
import type { SeqViewer } from './viewer.js';

// The viewer's frame clock (_spec/README.md "Composition"): one frame asked for at a time, whatever asked (`@molstar/ui/dom`'s
// frame callbacks, on `DomHost`); the view settled once it rests `canvas.settleMs` after its last move; the picture at
// rest (`rested`) once the settled frame is drawn, or `canvas.settleMs` after a change while the view was still.

export class FrameClock {
  private _frameId: number | undefined;
  /** A re-hit runs (the frame's): its own hover change is drawn by the frame that made it, none asked for. */
  rehitting = false;
  /** The view is still: a frame laid out now is settled. */
  quiet = true;
  /** The view came to rest: `rested` moves once the frame that lays it out settled is drawn. */
  restDue = false;
  private _settleTimer: unknown;
  private _restTimer: unknown;
  private readonly _rested = new State(0);
  readonly rested: ReadonlyState<number> = this._rested.asReadonly;

  constructor(private readonly viewer: SeqViewer) {}

  /**
   * Asks for a frame: one a frame, whatever asked; `dispose` cancels it. The callback is made once: a closure here made
   * V8 allocate this method's context (40 bytes) on every call, an early return's too, once optimized code had asked.
   */
  schedule(): void {
    if (this._frameId !== undefined || this.viewer.retired || this.rehitting) return;
    this._frameId = requestFrameCallback(this._frame);
  }

  private readonly _frame = (): void => {
    this._frameId = undefined;
    if (!this.viewer.retired) this.viewer.draw();
  };

  /** The view moved: it is not settled until it rests for `canvas.settleMs`. */
  moved(): void {
    if (this.viewer.retired) return;
    this.quiet = false;
    this.restDue = false;
    DomHost.clearTimeout(this._settleTimer);
    this._settleTimer = DomHost.setTimeout(() => {
      this._settleTimer = undefined;
      this.quiet = true;
      this.restDue = true;
      this.schedule();
    }, this.viewer.config.canvas.settleMs);
  }

  /** The picture changed while the view is quiet (data arrived, the height followed): it rests `canvas.settleMs` later. */
  changed(): void {
    if (this.viewer.retired || !this.quiet) return;
    DomHost.clearTimeout(this._restTimer);
    this._restTimer = DomHost.setTimeout(() => {
      this._restTimer = undefined;
      if (this.quiet && !this.viewer.retired) this._rested.set(this._rested.value + 1);
    }, this.viewer.config.canvas.settleMs);
  }

  /** A frame was drawn: when it is the settled picture, it rests now, once (a layout of it waits for no second rest). */
  drawn(): void {
    if (!this.restDue) return;
    this.restDue = false;
    DomHost.clearTimeout(this._restTimer);
    this._restTimer = undefined;
    this._rested.set(this._rested.value + 1);
  }

  /** Cancels the frame asked for and the timers (the viewer retired). */
  dispose(): void {
    if (this._frameId !== undefined) cancelFrameCallback(this._frameId);
    this._frameId = undefined;
    DomHost.clearTimeout(this._settleTimer);
    this._settleTimer = undefined;
    DomHost.clearTimeout(this._restTimer);
    this._restTimer = undefined;
  }
}
