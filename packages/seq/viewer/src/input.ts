/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Hit as RenderHit } from '@molstar/seqx-core/render';
import { type ModifierState, modifiersOf } from '@molstar/seqx-viewer/bindings';
import { changeKindOf } from './commands.js';
import {
  type HoverBehavior,
  HoverProvider,
  type NavigateBehavior,
  NavigateProvider,
  type SelectBehavior,
  SelectProvider,
} from '@molstar/seqx-viewer/behaviors';
import { type Hit, type HitTarget, hitTarget, offsetOf } from '@molstar/seqx-viewer/models';
import type { SeqViewer } from './viewer.js';

// The viewer's input (_spec/README.md "Input"): pointer, wheel and key events in canvas coordinates (CSS px from the
// canvas's top left), DOM-free, so the mounted UI only translates events and a test drives the same calls. Every
// gesture goes through the binding table: what a drag, a wheel or a pinch does is the table's, the behaviors do it.
// Per pointer move nothing is allocated: one hit record and one modifiers record are reused.

/** The modifier keys of an event, and where it is. */
export interface KeyState {
  readonly shiftKey?: boolean;
  readonly metaKey?: boolean;
  readonly ctrlKey?: boolean;
  readonly altKey?: boolean;
}

/** A click's state: its count (`detail`: 2 for the second click of a double-click). */
export interface ClickState extends KeyState {
  readonly detail?: number;
}

/** How far a pointer moves before a press is a drag, not a click (CSS px). */
const DragSlop = 3;

export interface WheelState extends KeyState {
  readonly deltaX: number;
  readonly deltaY: number;
  /** 0: pixels, 1: lines, 2: pages. */
  readonly deltaMode?: number;
}

const NoModifiers: ModifierState = {};

function platformIsMac(): boolean {
  return typeof navigator !== 'undefined' && /Mac|iPhone|iPad/.test(navigator.platform ?? '');
}

export class ViewerInput {
  /** Cmd is `mod` on macOS, Ctrl elsewhere. */
  mac = platformIsMac();
  private _space = false;
  private readonly _mods: ModifierState = {};
  private readonly _hit: { -readonly [K in keyof Hit]: Hit[K] } = {} as never;
  /** The pointer's last place while it rests over the canvas (EC-143's re-hit reads it). */
  private _x = Number.NaN;
  private _y = Number.NaN;
  private _drag:
    | { x: number; x0: number; section: string; target: HitTarget; select: boolean; moved: boolean }
    | undefined;
  /** A drag moved before the button went up: the click the platform sends after it is not one. */
  private _suppressClick = false;
  /** Space was held while a drag began: its release is no tap (no toggle). */
  private _spaceUsed = false;

  constructor(private readonly viewer: SeqViewer) {}

  private get hoverBehavior(): HoverBehavior | undefined {
    return this.viewer.behavior(HoverProvider);
  }

  private get navigateBehavior(): NavigateBehavior | undefined {
    return this.viewer.behavior(NavigateProvider);
  }

  private get selectBehavior(): SelectBehavior | undefined {
    return this.viewer.behavior(SelectProvider);
  }

  /** Whether the pointer is over the canvas. */
  get resting(): boolean {
    return !Number.isNaN(this._x);
  }

  /** The section of the layout at a y (a pointer over its label, its tracks or the gaps within), by address. */
  sectionAt(y: number): string | undefined {
    const layout = this.viewer.frame?.layout;
    if (!layout) return undefined;
    for (const s of layout.sections) if (y >= s.y && y < s.y + s.height) return s.section.address;
    return undefined;
  }

  /** The hit under a point as the hover model takes it (the section is read from the layer), in the reused record. */
  private _hitAt(x: number, y: number): Hit | undefined {
    const h: RenderHit | undefined = this.viewer.hit(x, y);
    if (!h) return undefined;
    const out = this._hit as Record<string, unknown>;
    out.layer = h.layer;
    out.system = h.system;
    out.position = h.position;
    out.element = h.element;
    out.extent = h.extent;
    return this._hit;
  }

  private _modifiers(e: KeyState): ModifierState {
    return modifiersOf(e, this.mac, this._space, this._mods);
  }

  /** The pointer moved to (x, y): a drag pans, and the hover follows through the bindings. */
  move(x: number, y: number, e: KeyState = {}): void {
    if (this.viewer.retired) return;
    this._x = x;
    this._y = y;
    const mods = this._modifiers(e);
    const d = this._drag;
    const hit = this._hitAt(x, y);
    if (d) {
      if (!d.moved && Math.abs(x - d.x0) >= DragSlop) d.moved = true;
      if (d.select) {
        // the range from where it began to the pointer, read from the view even past the data
        if (d.moved) this.selectBehavior?.dragTo(d.section, x - this.viewer.dataX, hit);
      } else {
        const dx = x - d.x;
        d.x = x;
        this.navigateBehavior?.drag(d.section, dx, mods, d.target);
      }
    }
    this.hoverBehavior?.move(hit, mods);
  }

  /** The pointer left the canvas: the hover it owned goes, and no re-hit follows. */
  leave(): void {
    this._x = this._y = Number.NaN;
    this.hoverBehavior?.leave();
  }

  /**
   * A resting pointer hit again (EC-143): the view moved under it, so what it is over is not what it was. Once a frame
   * (the viewer's draw asks after a layout); true when the hover changed.
   */
  rehit(): boolean {
    if (this.viewer.retired || !this.resting || this._drag) return false;
    const before = this.viewer.hover.version.value;
    this.hoverBehavior?.move(this._hitAt(this._x, this._y), this._mods);
    return this.viewer.hover.version.value !== before;
  }

  /** The primary button went down at (x, y): a drag begins when the bindings give one here. */
  down(x: number, y: number, e: KeyState = {}): boolean {
    if (this.viewer.retired) return false;
    // a new press: a click after it is its own
    this._suppressClick = false;
    const section = this.sectionAt(y);
    if (section === undefined) return false;
    const mods = this._modifiers(e);
    const hit = this._hitAt(x, y);
    const target = hitTarget(hit);
    const action = this.viewer.bindings.resolve('drag', mods, target);
    if (action === undefined) return false;
    // a drag-select begins on what was hit; it selects once the pointer moves past the slop (else it is a click)
    const select = action === 'drag-select' && !!this.selectBehavior?.beginDrag(section, hit, mods, target);
    if (action === 'drag-select' && !select) return false;
    if (this._space) this._spaceUsed = true;
    this._drag = { x, x0: x, section, target, select, moved: false };
    return true;
  }

  /**
   * The button went up (or the gesture was cancelled). The platform calls it more than once for one press (`pointerup`,
   * then `lostpointercapture`): only the call that ends a drag touches what follows it.
   */
  up(): void {
    const d = this._drag;
    if (!d) return;
    if (d.select) this.selectBehavior?.endDrag();
    this._suppressClick = d.moved;
    this._drag = undefined;
  }

  /**
   * A click at (x, y) through the bindings (select, extend, toggle, clear): what was hit is the granularity. The click
   * that ends a drag, and the second click of a double-click, are not clicks. True when a binding took it.
   */
  click(x: number, y: number, e: ClickState = {}): boolean {
    if (this.viewer.retired) return false;
    if (this._suppressClick) {
      this._suppressClick = false;
      return false;
    }
    if ((e.detail ?? 1) > 1) return false;
    const select = this.selectBehavior;
    if (!select) return false;
    const section = this.sectionAt(y);
    const hit = this._hitAt(x, y);
    // the behavior keeps nothing of the reused record but what it copies
    return select.click(
      section,
      hit ? { ...hit, ...(section !== undefined ? { section } : {}) } : undefined,
      this._modifiers(e),
    );
  }

  get dragging(): boolean {
    return this._drag !== undefined;
  }

  /**
   * A wheel at (x, y). A trackpad's pinch arrives as a wheel with Ctrl held and small pixel deltas: it is the `pinch`
   * gesture (`pointer.pinchFromCtrlWheel`), its zoom factor `exp(-deltaY · pinchSpeed)`; any other wheel is the `wheel`
   * gesture. True when a binding took it (the page should not scroll).
   */
  wheel(x: number, y: number, e: WheelState): boolean {
    if (this.viewer.retired) return false;
    const section = this.sectionAt(y);
    const nav = this.navigateBehavior;
    if (section === undefined || !nav) return false;
    const { pointer } = this.viewer.config;
    // lines are `pointer.lineHeight` px; a page is the viewport's own extent on each axis
    const ux = e.deltaMode === 1 ? pointer.lineHeight : e.deltaMode === 2 ? this.viewer.width : 1;
    const uy = e.deltaMode === 1 ? pointer.lineHeight : e.deltaMode === 2 ? this.viewer.height : 1;
    const dx = e.deltaX * ux;
    const dy = e.deltaY * uy;
    const target = hitTarget(this._hitAt(x, y));
    // over the label column the zoom is about the section's first visible position, never before it
    const width = this.viewer.view.get(section)?.width ?? 0;
    const local = Math.min(Math.max(x - this.viewer.dataX, 0), width);
    // on macOS Ctrl is no zoom modifier, so a Ctrl+wheel is a pinch whatever its delta; elsewhere a large delta is a mouse
    const pinch =
      e.ctrlKey && pointer.pinchFromCtrlWheel && (this.mac || (!e.deltaMode && Math.abs(dy) <= pointer.pinchMaxDelta));
    if (pinch) {
      // the pinch's own modifiers: none (the Ctrl is the gesture's, not a key held)
      const taken = this.viewer.bindings.resolve('pinch', NoModifiers, target) !== undefined;
      if (taken) nav.pinch(section, Math.exp(-dy * pointer.pinchSpeed), local, NoModifiers, target);
      this.viewer.schedule();
      return taken;
    }
    const mods = this._modifiers(e);
    const taken = this.viewer.bindings.resolve('wheel', mods, target) !== undefined;
    if (taken) nav.wheel(section, dx, dy, local, mods, target);
    // the frame asked for after the behavior's own, so a wheel's deltas apply in the frame that draws them
    this.viewer.schedule();
    return taken;
  }

  /** A double-click: the bound command, with the element or position under the pointer as what to focus. */
  doubleClick(x: number, y: number, e: KeyState = {}): boolean {
    if (this.viewer.retired) return false;
    const section = this.sectionAt(y);
    if (section === undefined) return false;
    const mods = this._modifiers(e);
    const hit = this._hitAt(x, y);
    const action = this.viewer.bindings.resolve('double-click', mods, hitTarget(hit));
    if (action === undefined) return false;
    const focus = hit ? this._focusOf(section, hit) : undefined;
    this.viewer.runCommand(action, { section, ...(focus ? { focus } : {}) });
    return true;
  }

  /** What to focus for a hit: the element's extent, else the position's cell, in the section's offsets. */
  private _focusOf(section: string, hit: Hit): { start: number; end: number } | undefined {
    const axis = this.viewer.view.resolved.value?.sections.find((s) => s.address === section)?.axis;
    if (!axis) return undefined;
    const part = hit.extent?.parts.find((p) => p.system === hit.system);
    if (part && part.ranges.length >= 2) {
      const start = offsetOf(axis, part.ranges[0]);
      const end = offsetOf(axis, part.ranges[part.ranges.length - 1]);
      return start !== undefined && end !== undefined ? { start, end } : undefined;
    }
    if (hit.position === undefined) return undefined;
    const start = offsetOf(axis, Math.floor(hit.position));
    return start === undefined ? undefined : { start, end: start + 1 };
  }

  /** Whether a binding takes Space as a modifier: the page's own Space (a scroll) is then the viewer's while it has focus. */
  get spaceIsModifier(): boolean {
    return this.viewer.bindings.bindings.some((b) => Array.isArray(b.modifiers) && b.modifiers.includes('space'));
  }

  /** The Space key, held for a drag that pans on any target. */
  setSpace(held: boolean): void {
    if (held && !this._space) this._spaceUsed = false;
    this._space = held;
  }

  /**
   * Space released: a tap (no drag began while it was held) is the Space key's command (the cursor's toggle); a hold
   * that panned is not. True when a command ran.
   */
  spaceUp(e: KeyState = {}): boolean {
    const tap = this._space && !this._spaceUsed;
    this.setSpace(false);
    return tap ? this.key({ ...e, key: ' ' }) : false;
  }

  /**
   * A key: the first command whose shortcut it is and that can run now (Esc cancels a gesture's preview first, then
   * clears the selection), on the section under the pointer (else the first). True when one ran. Keys go to the
   * commands' shortcuts; the keymap of `@molstar/ui` replaces this when the commands join a manager.
   */
  key(e: KeyState & { key: string }): boolean {
    if (this.viewer.retired) return false;
    const section =
      (this.resting ? this.sectionAt(this._y) : undefined) ?? this.viewer.view.resolved.value?.sections[0]?.address;
    const ctx = { section };
    // a gesture in progress first (a preview's command), then the rest in the spec's order
    for (const previews of [true, false])
      for (const c of this.viewer.commands) {
        const k = c.shortcut;
        if (!k || k.key !== e.key || (changeKindOf(c.id) === 'preview') !== previews) continue;
        if (!!k.alt !== !!e.altKey || !!k.shift !== !!e.shiftKey) continue;
        if (!!k.mod !== (this.mac ? !!e.metaKey : !!e.ctrlKey)) continue;
        const can = c.when?.(ctx as never);
        if (can !== undefined && can !== true) continue;
        this.viewer.runCommand(c.id, ctx);
        return true;
      }
    return false;
  }
}
