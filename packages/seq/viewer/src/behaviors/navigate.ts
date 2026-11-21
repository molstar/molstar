/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The `navigate` behavior (_spec/README.md "Behaviors"): drags, wheels and pinches to viewports through the bindings
// (`pan`, `zoom`), and the view's moves as commands for the keymap and the menu. Wheel deltas within a frame are
// applied once, at the frame (`DomHost`), when `navigate.coalesceWheel`. A move asked of a follower is the view model's
// to carry to its leader, ignore or detach (`follow.followerMoves`). A focus and a fit are jumps (the view's history:
// `view.back`, `view.forward`); a double-click on what is focused goes back; a step zoom is about the keyboard cursor,
// the pointer or the centre.

import type { CommandDefinition, CommandResult } from '@molstar/ui/commands/types';
import { cancelFrameCallback, requestFrameCallback } from '@molstar/ui/dom/frame';
import { type BindingTable, type KeyChord, type ModifierState, parseChord } from '@molstar/seqx-viewer/bindings';
import { DefaultViewerConfig, type ViewerConfig } from '@molstar/seqx-viewer/config';
import type { Location } from '@molstar/seqx-core/query';
import { blockSpan } from '@molstar/seqx-core/view';
import { blockOfUnit, blockSwitchRange } from '@molstar/seqx-viewer/block-switch';
import type { HitTarget, ViewModel } from '@molstar/seqx-viewer/models';

/**
 * Which section a command acts on (the one under the pointer or with focus), what to focus (offsets), and for a step
 * zoom the offset it is about (`at`: the pointer's; `atY`, down a grid's y) and its factor (`step`: a double-click's).
 * `cell`: the focus is a grid's cell, which zooms a step about `at` (and `atY`) instead of a jump to its positions.
 */
export interface NavigateContext {
  readonly section: string | undefined;
  readonly focus?: { readonly start: number; readonly end: number };
  /** What the focus names (an element's extent, a position), told to the view's last jump: a host reads it. */
  readonly location?: Location;
  readonly at?: number;
  readonly atY?: number;
  readonly cell?: boolean;
  readonly step?: number;
}

/** Where a step zoom with no pointer is about: the keyboard cursor's offset on a section, when it has one there. */
export type CursorOffset = (section: string) => number | undefined;

export interface NavigateBehavior {
  /**
   * A drag by `dx`, `dy` pixels on a section (a window pans across, a text scrolls); with `units`, by display units
   * the placement read under the pointer (round a ring: the window follows the pointer's angle). True when it could move.
   */
  drag(section: string, dx: number, modifiers: ModifierState, target: HitTarget, dy?: number, units?: number): boolean;
  /**
   * A wheel's deltas (pixels) at `x` (pixels from the section's left edge); `at`, the display unit under the pointer
   * from the placement (what a zoom keeps where it is), else read from the window at `x`.
   */
  wheel(
    section: string,
    dx: number,
    dy: number,
    x: number,
    modifiers: ModifierState,
    target: HitTarget,
    at?: number,
    atY?: number,
  ): void;
  /** A pinch scaling by `scale` about `x` (or the unit `at`; on a grid's y, `atY`). */
  pinch(
    section: string,
    scale: number,
    x: number,
    modifiers: ModifierState,
    target: HitTarget,
    at?: number,
    atY?: number,
  ): boolean;
  /** The view's moves as commands (keys through the keymap; every move reachable by a command). */
  readonly commands: readonly CommandDefinition<NavigateContext>[];
  /** Drops pending wheel deltas; later frames do nothing. */
  dispose(): void;
}

interface Pending {
  pan: number;
  /** A grid's vertical pan (its y moves apart from x). */
  panY: number;
  zoom: number;
  x: number;
  at: number | undefined;
  atY: number | undefined;
  /** The frame callback asked for, until it runs. */
  frame: number | undefined;
}

export function navigateBehavior(
  view: ViewModel,
  bindings: BindingTable,
  config: ViewerConfig = DefaultViewerConfig,
  cursorAt?: CursorOffset,
): NavigateBehavior {
  const { wheelZoomStep, keyPanFraction, coalesceWheel } = config.navigate;
  /** A command's keys: the configured chords, the first its shortcut. */
  const keysOf = (id: string): KeyChord | undefined => parseChord(config.view.keys[id]?.[0] ?? '');
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
    if (p.pan || p.panY) view.panPx(section, p.pan, p.panY, p.x);
    if (p.zoom) view.zoom(section, wheelZoomStep ** p.zoom, p.at ?? offsetAt(section, p.x), p.atY);
    p.pan = p.panY = p.zoom = 0;
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
    move: (section: string, ctx: NavigateContext) => boolean | string,
    done: string,
    when: (ctx: NavigateContext) => true | string = movable,
  ): CommandDefinition<NavigateContext> => {
    const shortcut = keysOf(id);
    return {
      id,
      title,
      category: 'View',
      ...(shortcut ? { shortcut: { ...shortcut } } : {}),
      // a single key works where the viewer has focus, never globally
      ...(shortcut && !shortcut.alt && !shortcut.mod && !shortcut.shift ? { scope: 'canvas' } : {}),
      when,
      run: (_, ctx): CommandResult => {
        const r = ctx.section !== undefined ? move(ctx.section, ctx) : false;
        return { message: typeof r === 'string' ? r : r ? done : '' };
      },
    };
  };

  /**
   * A step zoom about the context's offset (a double-click's pointer), else the keyboard cursor's on the section when
   * it is in view, else the centre: animated, never a step of the history (ux: zoom by a step animates).
   */
  const stepZoom = (section: string, ctx: NavigateContext, factor: number) => {
    const v = view.get(section);
    let anchor = ctx.at;
    if (anchor === undefined && v) {
      const c = cursorAt?.(section);
      if (c !== undefined && c >= v.start && c <= v.end) anchor = c;
    }
    return view.jump(() => view.zoom(section, factor, anchor), false, 'zoom') !== 'none';
  };

  /**
   * A double-click's focus: a jump; on what the view already shows focused, back (`view.backOnFocused`). On a grid's
   * cell a step zoom about it on both axes, their ratio kept, a step of the history each (again: further; back: undone)
   * (_spec/README.md "Input"): a jump to both positions would re-centre, never zoom.
   */
  const focus = (section: string, ctx: NavigateContext): boolean | string => {
    if (ctx.cell && ctx.at !== undefined)
      return view.jump(() => view.zoom(section, ctx.step ?? wheelZoomStep, ctx.at, ctx.atY), true, 'zoom') !== 'none';
    const r = view.jump(() => view.setRange(section, ctx.focus!.start, ctx.focus!.end), true, 'focus', ctx.location);
    if (r === 'same' && config.view.backOnFocused && view.canBack) return view.back() ? 'Back' : false;
    return r !== 'none';
  };

  const panBy = (section: string, fraction: number) => {
    const v = view.get(section);
    return !!v && view.pan(section, fraction * (v.end - v.start));
  };

  /** Scrolls by a line: the anchor moved by the layout's step up or down from it (up from the first line: its top). */
  const lineBy = (section: string, dir: 'up' | 'down') => {
    const v = view.get(section);
    if (!v || v.pending) return false;
    const from = Math.floor(v.start);
    const to = view.step(section, from, dir);
    if (to < 0) return dir === 'up' && v.start > 0 && view.pan(section, -v.start);
    return view.pan(section, to - from);
  };

  /** A section's blocks in a row (a genome's chromosomes): none for a stack or a section of one system. */
  const rowBlocks = (section: string | undefined) => {
    const s = section === undefined ? undefined : view.resolved.value?.sections.find((x) => x.address === section);
    return s?.blocks && s.blocks.length > 1 && s.params.arrange !== 'stack' ? s.blocks : undefined;
  };

  /**
   * The next or the previous block of a row (`block.next`, `block.previous`; genome coordinates): from the block holding
   * the window's centre, a jump to that block by the switching rule (`blockSwitchRange`: from a block shown whole the
   * next whole, zoomed in the same width at the same relative place; a step of the history); none past the last or
   * before the first.
   */
  const blockStep = (section: string, dir: 1 | -1): boolean | string => {
    const blocks = rowBlocks(section);
    const v = view.get(section);
    if (!blocks || !v) return false;
    const b = blockOfUnit(blocks, (v.start + v.end) / 2);
    let n = b + dir;
    while (n >= 0 && n < blocks.length && !(blockSpan(blocks[n]) > 0)) n += dir;
    if (n < 0 || n >= blocks.length) return dir > 0 ? 'At the last block' : 'At the first block';
    const r = blockSwitchRange(blocks, v, n, v.width > 0 ? (v.end - v.start) / v.width : 0);
    return view.jump(() => view.setRange(section, r[0], r[1]), true, 'move') !== 'none';
  };
  const hasBlocks = (ctx: NavigateContext): true | string => {
    const m = movable(ctx);
    return m !== true ? m : rowBlocks(ctx.section) ? true : 'No row of blocks here';
  };

  const commands: CommandDefinition<NavigateContext>[] = [
    command(
      'view.fit',
      'Fit the view',
      (s) => view.jump(() => view.fit(s), true, 'fit') !== 'none',
      'Fitted',
      (ctx) => (ctx.section === undefined ? 'No section to fit' : true),
    ),
    command('view.focus', 'Focus', focus, 'Focused', (ctx) => (ctx.focus ? movable(ctx) : 'Nothing to focus')),
    // back and forward restore every leading section together: they need no section
    command(
      'view.back',
      'Back',
      () => view.back(),
      'Back',
      () => (view.canBack ? true : 'No view to go back to'),
    ),
    command(
      'view.forward',
      'Forward',
      () => view.forward(),
      'Forward',
      () => (view.canForward ? true : 'No view to go forward to'),
    ),
    command('view.zoom-in', 'Zoom in', (s, ctx) => stepZoom(s, ctx, ctx.step ?? wheelZoomStep), 'Zoomed in'),
    command('view.zoom-out', 'Zoom out', (s, ctx) => stepZoom(s, ctx, 1 / (ctx.step ?? wheelZoomStep)), 'Zoomed out'),
    command('view.pan-left', 'Pan left', (s) => panBy(s, -keyPanFraction), 'Panned left'),
    command('view.pan-right', 'Pan right', (s) => panBy(s, keyPanFraction), 'Panned right'),
    // a line in a layout of lines (a text scrolls a line); a window has one line, so nothing
    command('view.pan-up', 'Scroll a line up', (s) => lineBy(s, 'up'), 'Scrolled up'),
    command('view.pan-down', 'Scroll a line down', (s) => lineBy(s, 'down'), 'Scrolled down'),
    command('block.next', 'Next block', (s) => blockStep(s, 1), 'Next block', hasBlocks),
    command('block.previous', 'Previous block', (s) => blockStep(s, -1), 'Previous block', hasBlocks),
  ];

  return {
    commands,
    drag(section, dx, modifiers, target, dy = 0, units) {
      if (disposed || bindings.resolve('drag', modifiers, target) !== 'pan') return false;
      // round a ring the window follows the pointer's angle (its units); elsewhere dragging content right (down) shows
      // what is left of (above) it
      return units !== undefined ? view.pan(section, units) : view.panPx(section, -dx, -dy);
    },
    wheel(section, dx, dy, x, modifiers, target, at, atY) {
      if (disposed) return;
      const action = bindings.resolve('wheel', modifiers, target);
      if (action !== 'pan' && action !== 'zoom') return;
      let p = pending.get(section);
      if (!p) pending.set(section, (p = { pan: 0, panY: 0, zoom: 0, x, at, atY, frame: undefined }));
      // the larger delta pans; with Shift a vertical wheel pans sideways too; a grid's y takes the vertical delta
      // the pointer's px for the pan too (the latest of the frame's)
      if (action === 'pan') p.x = x;
      if (action === 'pan' && view.movesY(section) && !modifiers.shift) {
        p.pan += dx;
        p.panY += dy;
      } else if (action === 'pan') p.pan += Math.abs(dx) > Math.abs(dy) ? dx : dy;
      else {
        p.zoom += -dy / 100;
        p.x = x;
        p.at = at;
        p.atY = atY;
      }
      if (!coalesceWheel) apply(section);
      else if (p.frame === undefined) p.frame = requestFrameCallback(() => apply(section));
    },
    pinch(section, scale, x, modifiers, target, at, atY) {
      return (
        !disposed &&
        bindings.resolve('pinch', modifiers, target) === 'zoom' &&
        view.zoom(section, scale, at ?? offsetAt(section, x), atY)
      );
    },
    dispose() {
      disposed = true;
      for (const p of pending.values()) if (p.frame !== undefined) cancelFrameCallback(p.frame);
      pending.clear();
    },
  };
}
