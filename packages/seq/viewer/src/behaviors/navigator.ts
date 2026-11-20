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
// One gesture is one step of the view's history: a click centres with its own step, and the double-click it begins
// (the platform's, by the click's count) moves within that step; a drag pushes the viewports at its press on release
// only when the view went farther than one window width (`navigator.dragStep: far`).

import { State, type ReadonlyState } from '@molstar/core/reactive';
import type { CommandDefinition, CommandResult } from '@molstar/ui/commands/types';
import type { WindowPart } from '@molstar/seqx-core/render';
import { axisTotal } from '@molstar/seqx-core/view';
import { DefaultViewerConfig, type ViewerConfig } from '@molstar/seqx-viewer/config';
import type { ViewModel, ViewSnapshot } from '@molstar/seqx-viewer/models';
import { windowBehavior, type WindowBehavior } from './window.js';

/** What the navigator's commands act on: the section, and the display units it centres on or zooms to. */
export interface NavigatorContext {
  readonly section: string | undefined;
  readonly at?: number;
  readonly start?: number;
  readonly end?: number;
}

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
  /**
   * A click (a press that never moved; `detail` its count, the platform's): outside the window (`at` given) the view
   * centres there, a jump with its own step; on the window (`at` undefined) nothing. A second click (`detail` 2) is a
   * double-click's: nothing, the double-click follows.
   */
  click(section: string, at: number | undefined, detail?: number): boolean;
  /**
   * A double-click at `at`: fit, or fitted, a window `units` wide centred there; within its first click's step when
   * that click centred the view (one gesture, one step), else a step of its own.
   */
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
  /** The section the last click centred with a step of its own: a double-click it begins moves within that step. */
  let clicked: string | undefined;

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
      // a drag is a gesture of its own: a double-click no longer folds into the click before it
      clicked = undefined;
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
    click(section, at, detail = 1) {
      if (detail > 1) return false;
      clicked = undefined;
      if (at === undefined) return false;
      const moved = center(section, at);
      if (moved) clicked = section;
      return moved;
    },
    doubleClick(section, at, units) {
      const mode = config.navigator.doubleClick;
      const folded = clicked === section;
      clicked = undefined;
      if (mode === 'none') return false;
      const fitted = mode !== 'fit' && isFitted(view, section);
      return fitted
        ? zoomTo(section, at - units / 2, at + units / 2, !folded)
        : view.jump(() => view.fit(section), !folded, { kind: 'fit' }) === 'moved';
    },
    dispose() {
      clicked = undefined;
    },
  };
}

/** Whether a section shows its whole axis. */
function isFitted(view: ViewModel, section: string): boolean {
  const v = view.get(section);
  const axis = view.resolved.value?.sections.find((s) => s.address === section)?.axis;
  const total = axis && axisTotal(axis);
  return v !== undefined && total !== undefined && v.start <= 1e-9 && v.end >= total - 1e-9;
}
