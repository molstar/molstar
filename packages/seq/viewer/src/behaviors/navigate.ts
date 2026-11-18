/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The `navigate` behavior (_spec/README.md "Behaviors"): drags, wheels and pinches to viewports through the bindings
// (`pan`, `zoom`), and the view's moves as commands for the keymap and the menu. Wheel deltas within a frame are
// applied once, at the frame (`DomHost`), when `navigate.coalesceWheel`. A move asked of a follower is the view model's
// to carry to its leader, ignore or detach (`follow.followerMoves`).

import type { CommandDefinition, CommandResult } from '@molstar/ui/commands/types';
import { cancelFrameCallback, requestFrameCallback } from '@molstar/ui/dom';
import type { BindingTable, ModifierState } from '@molstar/seqx-viewer/bindings';
import { DefaultViewerConfig, type ViewerConfig } from '@molstar/seqx-viewer/config';
import type { HitTarget, ViewModel } from '@molstar/seqx-viewer/models';

/** Which section a command acts on (the one under the pointer or with focus), and what to focus (offsets). */
export interface NavigateContext {
  readonly section: string | undefined;
  readonly focus?: { readonly start: number; readonly end: number };
}

export interface NavigateBehavior {
  /** A drag by `dx` pixels on a section; true when it moved the viewport. */
  drag(section: string, dx: number, modifiers: ModifierState, target: HitTarget): boolean;
  /** A wheel's deltas (pixels) at `x` (pixels from the section's left edge). */
  wheel(section: string, dx: number, dy: number, x: number, modifiers: ModifierState, target: HitTarget): void;
  /** A pinch scaling by `scale` about `x`. */
  pinch(section: string, scale: number, x: number, modifiers: ModifierState, target: HitTarget): boolean;
  /** The view's moves as commands (keys through the keymap; every move reachable by a command). */
  readonly commands: readonly CommandDefinition<NavigateContext>[];
  /** Drops pending wheel deltas; later frames do nothing. */
  dispose(): void;
}

interface Pending {
  pan: number;
  zoom: number;
  x: number;
  /** The frame callback asked for, until it runs. */
  frame: number | undefined;
}

export function navigateBehavior(
  view: ViewModel,
  bindings: BindingTable,
  config: ViewerConfig = DefaultViewerConfig,
): NavigateBehavior {
  const { wheelZoomStep, keyPanFraction, coalesceWheel } = config.navigate;
  // the wheel's pending deltas, one record per section, applied at the next frame
  const pending = new Map<string, Pending>();
  let disposed = false;

  const offsetAt = (section: string, x: number) => {
    const v = view.get(section);
    return v && v.width > 0 ? v.start + (x * (v.end - v.start)) / v.width : undefined;
  };

  const apply = (section: string) => {
    const p = pending.get(section);
    if (disposed || !p) return;
    p.frame = undefined;
    // a section gone since: its record goes too
    if (!view.get(section)) {
      pending.delete(section);
      return;
    }
    if (p.pan) view.panPx(section, p.pan);
    if (p.zoom) view.zoom(section, wheelZoomStep ** p.zoom, offsetAt(section, p.x));
    p.pan = p.zoom = 0;
  };

  const movable = (ctx: NavigateContext): true | string =>
    ctx.section === undefined
      ? 'No section to move'
      : view.canMove(ctx.section)
        ? true
        : 'The section follows another and does not move';

  const command = (
    id: string,
    title: string,
    move: (section: string, ctx: NavigateContext) => boolean,
    done: string,
    shortcut: CommandDefinition<NavigateContext>['shortcut'],
    when: (ctx: NavigateContext) => true | string = movable,
  ): CommandDefinition<NavigateContext> => ({
    id,
    title,
    category: 'View',
    shortcut,
    // a single key works where the viewer has focus, never globally
    ...(shortcut && !shortcut.alt && !shortcut.mod && !shortcut.shift ? { scope: 'canvas' } : {}),
    when,
    run: (_, ctx): CommandResult => ({ message: ctx.section !== undefined && move(ctx.section, ctx) ? done : '' }),
  });

  const panBy = (section: string, fraction: number) => {
    const v = view.get(section);
    return !!v && view.pan(section, fraction * (v.end - v.start));
  };

  const commands: CommandDefinition<NavigateContext>[] = [
    command(
      'view.fit',
      'Fit the view',
      (s) => view.fit(s),
      'Fitted',
      { key: '0' },
      (ctx) => (ctx.section === undefined ? 'No section to fit' : true),
    ),
    command(
      'view.focus',
      'Focus',
      (s, ctx) => view.setRange(s, ctx.focus!.start, ctx.focus!.end),
      'Focused',
      undefined,
      (ctx) => (ctx.focus ? movable(ctx) : 'Nothing to focus'),
    ),
    command('view.zoom-in', 'Zoom in', (s) => view.zoom(s, wheelZoomStep), 'Zoomed in', { key: '=' }),
    command('view.zoom-out', 'Zoom out', (s) => view.zoom(s, 1 / wheelZoomStep), 'Zoomed out', { key: '-' }),
    command('view.pan-left', 'Pan left', (s) => panBy(s, -keyPanFraction), 'Panned left', {
      key: 'ArrowLeft',
      alt: true,
    }),
    command('view.pan-right', 'Pan right', (s) => panBy(s, keyPanFraction), 'Panned right', {
      key: 'ArrowRight',
      alt: true,
    }),
  ];

  return {
    commands,
    drag(section, dx, modifiers, target) {
      // dragging content right shows what is left of it
      return !disposed && bindings.resolve('drag', modifiers, target) === 'pan' && view.panPx(section, -dx);
    },
    wheel(section, dx, dy, x, modifiers, target) {
      if (disposed) return;
      const action = bindings.resolve('wheel', modifiers, target);
      if (action !== 'pan' && action !== 'zoom') return;
      let p = pending.get(section);
      if (!p) pending.set(section, (p = { pan: 0, zoom: 0, x, frame: undefined }));
      // the larger delta pans; with Shift a vertical wheel pans sideways too
      if (action === 'pan') p.pan += Math.abs(dx) > Math.abs(dy) ? dx : dy;
      else {
        p.zoom += -dy / 100;
        p.x = x;
      }
      if (!coalesceWheel) apply(section);
      else if (p.frame === undefined) p.frame = requestFrameCallback(() => apply(section));
    },
    pinch(section, scale, x, modifiers, target) {
      return (
        !disposed &&
        bindings.resolve('pinch', modifiers, target) === 'zoom' &&
        view.zoom(section, scale, offsetAt(section, x))
      );
    },
    dispose() {
      disposed = true;
      for (const p of pending.values()) if (p.frame !== undefined) cancelFrameCallback(p.frame);
      pending.clear();
    },
  };
}
