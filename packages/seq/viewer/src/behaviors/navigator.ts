/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The `navigator` behavior (_spec/README.md "Behaviors", Navigator): a navigator's gestures on the `window` behavior,
// in the bound section's display units (the navigator model carries its px there and back). A drag on a window's body
// pans by hand, an end zooms about the other (no animation); a press outside the window is a click that centres there
// (`navigator.center`, a jump: animated) or, dragged, sweeps a new window and zooms to it on release
// (`navigator.zoom-to`, a jump; `navigator.dragOutside: pan` centres there and then pans); a double-click fits, or
// fitted zooms so the window is `navigator.zoomPx` navigator px wide; Esc puts the view back where the press found it.
// One gesture is one step of the view's history: a click centres at once and keeps the viewports it left until a
// double-click in time pushes them as its own step (or the time runs out); a drag pushes the viewports at its press
// on release only when the view went farther than one window width (`navigator.dragStep: far`).

import { State, type ReadonlyState } from '@molstar/core/reactive';
import type { CommandDefinition, CommandResult } from '@molstar/ui/commands/types';
import type { WindowPart } from '@molstar/seqx-core/render';
import { DefaultViewerConfig, type ViewerConfig } from '@molstar/seqx-viewer/config';
import type { ViewModel, ViewSnapshot } from '@molstar/seqx-viewer/models';
import { DomHost } from '@molstar/ui/dom';
import { windowBehavior, type WindowBehavior } from './window.js';

/** What the navigator's commands act on: the section, and the display units it centres on or zooms to. */
export interface NavigatorContext {
  readonly section: string | undefined;
  readonly at?: number;
  readonly start?: number;
  readonly end?: number;
}

/** How long a click waits for a second click before its step is pushed, ms (a double-click's time). */
export const DoubleClickMs = 400;

export interface NavigatorBehavior {
  /** The window drags (its own: a navigator's drag is not the canvas's). */
  readonly window: WindowBehavior;
  /** The section's view centred on display unit `at`, its width kept: a jump (animated), pushed when `record`. */
  center(section: string, at: number, record?: boolean): boolean;
  /** The section's view on `[start, end)`: a jump. */
  zoomTo(section: string, start: number, end: number, record?: boolean): boolean;
  /**
   * A press on a part at display unit `at` (`drawn`: the window as drawn when floored): a drag of the window (`body`,
   * `start`, `end`), or outside it a sweep (or, `dragOutside: pan`, a centre then a pan). True when a gesture began.
   */
  press(section: string, part: WindowPart, at: number, drawn?: readonly [number, number]): boolean;
  /** The pointer at display unit `at` while pressed; `moved`: past the drag slop (a press that never moves is a click). */
  dragTo(at: number): boolean;
  /** The release: a sweep zooms to what it swept; a far drag pushes its step. */
  release(): void;
  /** Esc: the gesture is cancelled, the view back where the press found it. */
  cancel(): void;
  /** A click outside the window at `at`: the view centres there, its step deferred for a double-click. */
  click(section: string, at: number): boolean;
  /** A double-click at `at`: fit, or fitted, a window `units` wide centred there; one step with its first click's. */
  doubleClick(section: string, at: number, units: number): boolean;
  /** The swept region while a sweep runs (display units, ordered). */
  readonly sweep: ReadonlyState<readonly [number, number] | undefined>;
  readonly pressed: boolean;
  readonly commands: readonly CommandDefinition<NavigatorContext>[];
  dispose(): void;
}

export function navigatorBehavior(view: ViewModel, config: ViewerConfig = DefaultViewerConfig): NavigatorBehavior {
  const win = windowBehavior(view, config);
  const sweep = new State<readonly [number, number] | undefined>(undefined);
  /** A press: its section, part, where it began, the viewports then, and whether it pans from outside (`pan`). */
  let press:
    | {
        section: string;
        part: WindowPart;
        at0: number;
        before: ViewSnapshot;
        width: number;
        center: number;
        sweeping: boolean;
      }
    | undefined;
  /** A click's step, deferred until a double-click or its time. */
  let deferred: { snapshot: ViewSnapshot; timer: unknown } | undefined;

  const flush = () => {
    const d = deferred;
    deferred = undefined;
    if (!d) return;
    DomHost.clearTimeout(d.timer as never);
    view.record(d.snapshot);
  };

  const width = (section: string) => {
    const v = view.get(section);
    return v ? v.end - v.start : 0;
  };

  const center = (section: string, at: number, record = true) => {
    const w = width(section);
    if (!(w > 0) || !Number.isFinite(at)) return false;
    return view.jump(() => view.setRange(section, at - w / 2, at + w / 2), record) === 'moved';
  };

  const zoomTo = (section: string, start: number, end: number, record = true) => {
    if (!Number.isFinite(start) || !Number.isFinite(end) || start === end) return false;
    return view.jump(() => view.setRange(section, Math.min(start, end), Math.max(start, end)), record) === 'moved';
  };

  const commands: CommandDefinition<NavigatorContext>[] = [
    {
      id: 'navigator.center',
      title: 'Centre there',
      category: 'View',
      when: (ctx) =>
        ctx.section !== undefined && ctx.at !== undefined && view.canMove(ctx.section) ? true : 'No window',
      run: (_, ctx): CommandResult => ({
        message: ctx.section !== undefined && ctx.at !== undefined && center(ctx.section, ctx.at) ? 'Centred' : '',
      }),
    },
    {
      id: 'navigator.zoom-to',
      title: 'Zoom to the region',
      category: 'View',
      when: (ctx) =>
        ctx.section !== undefined && ctx.start !== undefined && ctx.end !== undefined && view.canMove(ctx.section)
          ? true
          : 'No region',
      run: (_, ctx): CommandResult => ({
        message:
          ctx.section !== undefined &&
          ctx.start !== undefined &&
          ctx.end !== undefined &&
          zoomTo(ctx.section, ctx.start, ctx.end)
            ? 'Zoomed to the region'
            : '',
      }),
    },
  ];

  return {
    window: win,
    center,
    zoomTo,
    commands,
    sweep: sweep.asReadonly,
    get pressed() {
      return press !== undefined;
    },
    press(section, part, at, drawn) {
      const v = view.get(section);
      if (!v || v.pending || !view.canMove(section) || !Number.isFinite(at)) return false;
      // a gesture after a click: the click's step is its own
      flush();
      const before = view.snapshot();
      if (part === 'outside') {
        press = {
          section,
          part,
          at0: at,
          before,
          width: v.end - v.start,
          center: (v.start + v.end) / 2,
          sweeping: false,
        };
        return true;
      }
      if (!win.begin(section, part, at, drawn)) return false;
      press = {
        section,
        part,
        at0: at,
        before,
        width: v.end - v.start,
        center: (v.start + v.end) / 2,
        sweeping: false,
      };
      return true;
    },
    dragTo(at) {
      const p = press;
      if (!p || !Number.isFinite(at)) return false;
      if (p.part !== 'outside') return win.dragTo(at);
      if (config.navigator.dragOutside === 'pan') {
        // jump there, then the drag pans the window as its body
        if (!p.sweeping) {
          p.sweeping = true;
          const w = width(p.section);
          view.setRange(p.section, p.at0 - w / 2, p.at0 + w / 2);
          if (!win.begin(p.section, 'body', p.at0)) return false;
        }
        return win.dragTo(at);
      }
      p.sweeping = true;
      sweep.set([Math.min(p.at0, at), Math.max(p.at0, at)]);
      return true;
    },
    release() {
      const p = press;
      press = undefined;
      if (!p) return;
      if (p.part === 'outside' && config.navigator.dragOutside === 'sweep') {
        const s = sweep.value;
        sweep.set(undefined);
        // a sweep of some width zooms to it: one jump, its step the press's
        if (p.sweeping && s && s[1] > s[0]) {
          if (zoomTo(p.section, s[0], s[1], false)) view.record(p.before);
        }
        return;
      }
      win.end();
      const v = view.get(p.section);
      if (!v) return;
      const far = Math.abs((v.start + v.end) / 2 - p.center) > p.width;
      const step = config.navigator.dragStep;
      if (step === 'always' || (step === 'far' && far)) view.record(p.before);
    },
    cancel() {
      const p = press;
      press = undefined;
      sweep.set(undefined);
      if (!p) return;
      if (p.part !== 'outside') win.cancel();
      else if (p.sweeping && config.navigator.dragOutside === 'pan') {
        win.cancel();
        view.setRange(p.section, p.center - p.width / 2, p.center + p.width / 2);
      }
    },
    click(section, at) {
      // the click's step waits: a double-click in time pushes it as its own
      const snapshot = deferred?.snapshot ?? view.snapshot();
      if (deferred) DomHost.clearTimeout(deferred.timer as never);
      const moved = center(section, at, false);
      deferred = { snapshot, timer: DomHost.setTimeout(flush, DoubleClickMs) };
      return moved;
    },
    doubleClick(section, at, units) {
      const mode = config.navigator.doubleClick;
      if (mode === 'none') return false;
      const snapshot = deferred?.snapshot ?? view.snapshot();
      if (deferred) DomHost.clearTimeout(deferred.timer as never);
      deferred = undefined;
      const v = view.get(section);
      const axisFitted = v !== undefined && isFitted(view, section);
      const moved =
        mode === 'fit' || !axisFitted
          ? view.jump(() => view.fit(section), false) === 'moved'
          : zoomTo(section, at - units / 2, at + units / 2, false);
      view.record(snapshot);
      return moved;
    },
    dispose() {
      if (deferred) DomHost.clearTimeout(deferred.timer as never);
      deferred = undefined;
    },
  };
}

/** Whether a section shows its whole axis. */
function isFitted(view: ViewModel, section: string): boolean {
  const v = view.get(section);
  const axis = view.resolved.value?.sections.find((s) => s.address === section)?.axis;
  const r = axis?.ranges;
  if (!v || !r) return false;
  let total = 0;
  for (let k = 0; k + 1 < r.length; k += 2) total += r[k + 1] - r[k];
  return v.start <= 1e-9 && v.end >= total - 1e-9;
}
