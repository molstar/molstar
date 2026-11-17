/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The `navigate` behavior (_spec/README.md "Behaviors"): drags, wheels and pinches to viewports through the bindings
// (`pan`, `zoom`), and the view's moves as commands for the keymap and the menu. Wheel deltas within a frame are
// applied once, at the frame (`DomHost`), when `navigate.coalesceWheel`. Followers never lead: a move asked of one is
// the view model's to ignore or detach (EC-114).

import type { CommandDefinition } from '@molstar/ui/commands/types';
import { DomHost } from '@molstar/ui/dom';
import type { BindingTable, ModifierState } from '../bindings.js';
import { DefaultViewerConfig, type ViewerConfig } from '../config.js';
import type { HitTarget } from '../models/hover.js';
import type { ViewModel } from '../models/view.js';

/** Which section a command acts on: the one under the pointer or with focus. */
export interface NavigateContext {
  readonly section: string | undefined;
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
}

const NoModifiers: ModifierState = {};

export function navigateBehavior(
  view: ViewModel,
  bindings: BindingTable,
  config: ViewerConfig = DefaultViewerConfig,
): NavigateBehavior {
  const { wheelZoomStep, keyPanFraction, coalesceWheel } = config.navigate;
  // the wheel's pending deltas, one record per section, applied at the next frame
  const pending = new Map<string, { pan: number; zoom: number; x: number; scheduled: boolean }>();

  const offsetAt = (section: string, x: number) => {
    const v = view.get(section);
    return v && v.width > 0 ? v.start + (x * (v.end - v.start)) / v.width : undefined;
  };

  const apply = (section: string) => {
    const p = pending.get(section);
    if (!p) return;
    p.scheduled = false;
    if (p.pan) view.panPx(section, p.pan);
    if (p.zoom) view.zoom(section, wheelZoomStep ** p.zoom, offsetAt(section, p.x));
    p.pan = p.zoom = 0;
  };

  const step = (section: string | undefined, fn: (s: string) => boolean, done: string) => {
    if (section === undefined) throw new Error('No section to move');
    if (!fn(section)) throw new Error('The section cannot move');
    return { message: done };
  };

  const commands: CommandDefinition<NavigateContext>[] = [
    {
      id: 'view.fit',
      title: 'Fit the view',
      category: 'View',
      shortcut: { key: '0' },
      run: (_, ctx) => step(ctx.section, (s) => view.fit(s), 'Fitted'),
    },
    {
      id: 'view.zoom-in',
      title: 'Zoom in',
      category: 'View',
      shortcut: { key: '=' },
      run: (_, ctx) => step(ctx.section, (s) => view.zoom(s, wheelZoomStep), 'Zoomed in'),
    },
    {
      id: 'view.zoom-out',
      title: 'Zoom out',
      category: 'View',
      shortcut: { key: '-' },
      run: (_, ctx) => step(ctx.section, (s) => view.zoom(s, 1 / wheelZoomStep), 'Zoomed out'),
    },
    {
      id: 'view.pan-left',
      title: 'Pan left',
      category: 'View',
      shortcut: { key: 'ArrowLeft', alt: true },
      run: (_, ctx) => step(ctx.section, (s) => panBy(s, -keyPanFraction), 'Panned left'),
    },
    {
      id: 'view.pan-right',
      title: 'Pan right',
      category: 'View',
      shortcut: { key: 'ArrowRight', alt: true },
      run: (_, ctx) => step(ctx.section, (s) => panBy(s, keyPanFraction), 'Panned right'),
    },
  ];

  function panBy(section: string, fraction: number): boolean {
    const v = view.get(section);
    return !!v && view.pan(section, fraction * (v.end - v.start));
  }

  return {
    commands,
    drag(section, dx, modifiers, target) {
      // dragging content right shows what is left of it
      return bindings.resolve('drag', modifiers, target) === 'pan' && view.panPx(section, -dx);
    },
    wheel(section, dx, dy, x, modifiers = NoModifiers, target) {
      const action = bindings.resolve('wheel', modifiers, target);
      if (action !== 'pan' && action !== 'zoom') return;
      let p = pending.get(section);
      if (!p) pending.set(section, (p = { pan: 0, zoom: 0, x, scheduled: false }));
      if (action === 'pan') p.pan += Math.abs(dx) > Math.abs(dy) ? dx : dy;
      else {
        p.zoom += -dy / 100;
        p.x = x;
      }
      if (!coalesceWheel) apply(section);
      else if (!p.scheduled) {
        p.scheduled = true;
        DomHost.requestAnimationFrame(() => apply(section));
      }
    },
    pinch(section, scale, x, modifiers, target) {
      return bindings.resolve('pinch', modifiers, target) === 'zoom' && view.zoom(section, scale, offsetAt(section, x));
    },
  };
}
