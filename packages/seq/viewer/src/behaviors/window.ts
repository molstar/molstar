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
  /** A press on a window's part at display unit `at`: true when a drag begins (the window can move). */
  begin(section: string, part: WindowPart, at: number): boolean;
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
  let drag: { section: string; part: WindowPart; at0: number; start: number; end: number } | undefined;

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
    const b = positionAt(axis, Math.min(v.end, total), 'end');
    if (a === undefined || b === undefined) return;
    const start = Math.min(a, b) + 1;
    const end = Math.max(a, b);
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
    begin(section, part, at) {
      const v = view.get(section);
      if (!v || v.pending || !view.canMove(section) || !Number.isFinite(at)) return false;
      drag = { section, part, at0: at, start: v.start, end: v.end };
      return true;
    },
    dragTo(at) {
      const d = drag;
      if (!d || !Number.isFinite(at)) return false;
      if (d.part === 'body') {
        const delta = at - d.at0;
        return move(d.section, d.start + delta, d.end + delta);
      }
      // an end: the other held; an inward drag on a floored window keeps it (the floor holds its far end)
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
