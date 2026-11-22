/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The `window` behavior (_spec/README.md "Behaviors", Window): a window's drags and its slider, one for every window (a
// ring's, a navigator's). The window is its section's viewport: a drag on its `body` pans it, a drag on its `start` or
// `end` resizes it about the other end, never narrower than the floor (an inward drag on a floored window keeps it);
// Esc during a drag puts the view back where the press found it. The pointer comes in as the display unit the
// placement has under it (the layout's inverse: an angle on a ring), so the behavior knows no layout. As a slider:
// commands pan it a step and grow or shrink it about its middle (the focused window's arrows and Shift+arrows), and its
// `aria-valuetext` is said when its viewport settles (a drag's end, a step), never a frame.

import type { WindowPart } from '@molstar/seqx-core/render';
import { State, type ReadonlyState } from '@molstar/core/reactive';
import type { CommandDefinition, CommandResult } from '@molstar/ui/commands/types';
import { DefaultViewerConfig, type ViewerConfig } from '@molstar/seqx-viewer/config';
import { axisTotal, positionAt, type ViewModel } from '@molstar/seqx-viewer/models';

/** The section a window command moves. */
export interface WindowContext {
  readonly section: string | undefined;
}

/** A window as a slider says it: its section, its first and last positions shown (1-based), its system's length. */
export interface WindowValue {
  readonly section: string;
  readonly text: string;
  readonly start: number;
  readonly end: number;
  readonly total: number;
}

export interface WindowBehavior {
  /**
   * A press on a window's part at display unit `at`: true when a drag begins (the window can move). `drawn`: the window
   * as drawn, in display units, when the floor drew it wider than the truth (a navigator's): an end dragged inward then
   * keeps the view, outward it zooms out from the drawn end.
   */
  begin(section: string, part: WindowPart, at: number, drawn?: readonly [number, number]): boolean;
  /** The drag at display unit `at`: the body pans, an end resizes about the other. */
  dragTo(at: number): boolean;
  /** The drag ends: the slider says where the window settled. */
  end(): void;
  /** Esc: the drag is cancelled, the window back where the press found it. */
  cancel(): void;
  readonly dragging: boolean;
  /** What the slider says, set when a window settles (`aria-valuetext`), never a frame. */
  readonly value: ReadonlyState<WindowValue | undefined>;
  readonly commands: readonly CommandDefinition<WindowContext>[];
}

/** The least share of the axis a window keeps (with the most px a position, the floor in display units). */
const MinUnits = 1;

export function windowBehavior(view: ViewModel, config: ViewerConfig = DefaultViewerConfig): WindowBehavior {
  const value = new State<WindowValue | undefined>(
    undefined,
    (a, b) => a?.text === b?.text && a?.section === b?.section,
  );
  /** A drag: where it began, the pointer's last place and the units moved since (unwrapped on a cyclic axis's total). */
  let drag:
    | {
        section: string;
        part: WindowPart;
        at0: number;
        start: number;
        end: number;
        last: number;
        moved: number;
        total: number | undefined;
        drawn: readonly [number, number] | undefined;
      }
    | undefined;

  /** The window's width at least: its floor in display units (`viewport.maxPxPerPosition`), at least one. */
  const floorOf = (section: string) => {
    const v = view.get(section);
    return Math.max(MinUnits, v && v.width > 0 ? v.width / config.viewport.maxPxPerPosition : MinUnits);
  };

  const say = (section: string) => {
    const v = view.get(section);
    const axis = view.resolved.value?.sections.find((s) => s.address === section)?.axis;
    const total = axis && axisTotal(axis);
    if (!v || v.pending || !axis || total === undefined) return;
    const a = positionAt(axis, v.start, 'start');
    // a window over a cyclic axis's end runs on past the total: its end is before its start (an arc, never the rest)
    const over = axis.cyclic && v.end > total;
    const b = positionAt(axis, over ? v.end - total : Math.min(v.end, total), 'end');
    if (a === undefined || b === undefined) return;
    const arc = over || (axis.cyclic && a > b);
    const start = (arc ? a : Math.min(a, b)) + 1;
    const end = arc ? b : Math.max(a, b);
    value.set({
      section,
      start,
      end,
      total,
      text: `${start.toLocaleString('en-US')} to ${end.toLocaleString('en-US')} of ${total.toLocaleString('en-US')}`,
    });
  };

  const move = (section: string, start: number, end: number) => {
    const w = floorOf(section);
    if (end - start < w) {
      // the floor: grown about the end held (a resize) or the middle
      if (drag?.part === 'start') start = end - w;
      else if (drag?.part === 'end') end = start + w;
      else {
        const c = (start + end) / 2;
        start = c - w / 2;
        end = c + w / 2;
      }
    }
    return view.setRange(section, start, end);
  };

  const step = (section: string, f: (s: number, e: number) => [number, number]) => {
    const v = view.get(section);
    if (!v || v.pending) return false;
    const [s, e] = f(v.start, v.end);
    const moved = move(section, s, e);
    if (moved) say(section);
    return moved;
  };

  const { keyPanFraction, wheelZoomStep } = config.navigate;
  const command = (
    id: string,
    title: string,
    f: (s: number, e: number) => [number, number],
  ): CommandDefinition<WindowContext> => ({
    id,
    title,
    category: 'View',
    when: (ctx) => (ctx.section !== undefined && view.canMove(ctx.section) ? true : 'No window to move'),
    run: (_, ctx): CommandResult => ({ message: ctx.section !== undefined && step(ctx.section, f) ? title : '' }),
  });

  const commands = [
    command('window.pan-back', 'Window back', (s, e) => [s - (e - s) * keyPanFraction, e - (e - s) * keyPanFraction]),
    command('window.pan-on', 'Window on', (s, e) => [s + (e - s) * keyPanFraction, e + (e - s) * keyPanFraction]),
    command('window.grow', 'Window wider', (s, e) => {
      const c = (s + e) / 2;
      const h = ((e - s) * wheelZoomStep) / 2;
      return [c - h, c + h];
    }),
    command('window.shrink', 'Window narrower', (s, e) => {
      const c = (s + e) / 2;
      const h = (e - s) / wheelZoomStep / 2;
      return [c - h, c + h];
    }),
  ];

  return {
    commands,
    value: value.asReadonly,
    get dragging() {
      return drag !== undefined;
    },
    begin(section, part, at, drawn) {
      const v = view.get(section);
      if (part === 'outside' || !v || v.pending || !view.canMove(section) || !Number.isFinite(at)) return false;
      // unwrapped over the origin where the section's window wraps it (`view.wraps`: a cyclic axis in tracks and rings)
      const axis = view.resolved.value?.sections.find((s) => s.address === section)?.axis;
      const total = axis && view.wraps(section) ? axisTotal(axis) : undefined;
      const floored = drawn && drawn[1] - drawn[0] > v.end - v.start ? drawn : undefined;
      drag = { section, part, at0: at, start: v.start, end: v.end, last: at, moved: 0, total, drawn: floored };
      return true;
    },
    dragTo(at) {
      const d = drag;
      if (!d || !Number.isFinite(at)) return false;
      if (d.part === 'body') {
        // on a cyclic axis the pointer is unwrapped against its last place: a step over the origin is a short one on
        // (a unit after the last is the first), so the window keeps going past it
        if (d.total !== undefined) {
          const t = d.total;
          const step = ((((at - d.last + t / 2) % t) + t) % t) - t / 2;
          d.last = at;
          d.moved += step;
        } else d.moved = at - d.at0;
        return move(d.section, d.start + d.moved, d.end + d.moved);
      }
      // an end: the other held; an inward drag on a floored window keeps it (the floor holds its far end). On a cyclic
      // axis the pointer is unwrapped against its last place as the body's is (a ring's grip dragged over the origin),
      // from the grip's own turn
      if (d.total !== undefined) {
        const t = d.total;
        if (d.moved === 0 && d.last === d.at0) {
          const grip = d.part === 'start' ? d.start : d.end;
          d.last = d.at0 + Math.round((grip - d.at0) / t) * t;
        }
        const step = ((((at - d.last + t / 2) % t) + t) % t) - t / 2;
        d.last += step;
        d.moved = 1;
        at = d.last;
      }
      // a floored window (drawn wider than the truth): inward keeps the view, outward zooms out from the drawn end
      if (d.drawn && d.total === undefined) {
        const out = d.part === 'start' ? d.at0 - at : at - d.at0;
        if (d.part === 'start') return move(d.section, out > 0 ? d.drawn[0] - out : d.start, d.end);
        return move(d.section, d.start, out > 0 ? d.drawn[1] + out : d.end);
      }
      const held = d.part === 'start' ? d.end : d.start;
      const s = d.part === 'start' ? Math.min(at, held) : held;
      const e = d.part === 'start' ? held : Math.max(at, held);
      return move(d.section, s, e);
    },
    end() {
      const d = drag;
      drag = undefined;
      if (d) say(d.section);
    },
    cancel() {
      const d = drag;
      drag = undefined;
      if (d) view.setRange(d.section, d.start, d.end);
    },
  };
}

/**
 * Whether a section's window has nowhere to go along its positions (_spec/README.md "Behaviors", Window): it cannot move
 * (a follower under `ignore`, a full-state step), or it shows its whole axis and zoom-in is at its limit
 * (`viewport.maxPxPerPosition`). Idle is drawn quiet and is no target; it wakes when a zoom becomes possible.
 */
export function windowIdle(view: ViewModel, section: string, config: ViewerConfig = DefaultViewerConfig): boolean {
  const v = view.get(section);
  if (!v || v.pending || v.width <= 0 || !view.canMove(section)) return true;
  // (a loop, not `find`: a frame asks this per rail and allocates no closure)
  const sections = view.resolved.value?.sections;
  let total: number | undefined;
  if (sections)
    for (let i = 0; i < sections.length; i++) {
      const s = sections[i];
      if (s.address !== section) continue;
      total = s.axis && axisTotal(s.axis);
      break;
    }
  if (total === undefined) return true;
  const whole = v.start <= 0 && v.end >= total;
  return whole && total <= v.width / config.viewport.maxPxPerPosition;
}
