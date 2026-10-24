/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// DOM code both default UIs share (seqspec:interaction#the-dom-element): the popover both menus and the issue list open
// with, the keyboard-only focus indicator on a surface, and the placement of the chrome that follows the pointer, the
// keyboard cursor or the viewport (the tooltip, a track header's menu). Pointer-rate: these write styles only, never
// structure. Imports no React (the React entry imports these).

import { AxisGeometry } from '../axis.js';
import { type HeaderLine, headerLine, headerTrack, trackLine } from '../chrome.js';
import {
  HoverKind,
  type LegendRamp,
  type LegendSize,
  type SeqState,
  type SeqViewerContext,
  type Unsubscribable,
} from '../index.js';
import { cursorRect } from '../keyboard.js';
import { SizeUnit } from '../themes.js';
/**
 * A state's changes, not its replay: a state calls back with its value on subscribe (as a BehaviorSubject), which is not a
 * change. For a subscription whose first run is made some other way (a direct call, a `watch`).
 */
export function later<T>(s: { subscribe(fn: (v: T) => void): Unsubscribable }, fn: (v: T) => void): Unsubscribable {
  let first = true;
  return s.subscribe((v) => {
    if (first) first = false;
    else fn(v);
  });
}
/**
 * `fn` once now, then on each later change of any of the states: one function over several states runs once at mount,
 * not once per state's replay.
 */
export function watch(states: readonly Pick<SeqState<unknown>, 'subscribe'>[], fn: () => void): Unsubscribable[] {
  // Subscribed first: a state that changes while `fn` runs is not missed.
  const subs = states.map((s) => later(s, fn));
  fn();
  return subs;
}

// The last input's modality, shared by every viewer on the page: Tab sets it, a press clears it.
let keyboard = false;
let listeners = 0;
const onKey = (e: KeyboardEvent) => {
  if (e.key === 'Tab') keyboard = true;
};
const onPress = () => {
  keyboard = false;
};
/** `data-focus-visible` on the host while its surface holds focus that came from the keyboard. */
export function focusIndicator(host: HTMLElement, surface: HTMLElement): () => void {
  if (listeners++ === 0) {
    document.addEventListener('keydown', onKey, true);
    document.addEventListener('pointerdown', onPress, true);
  }
  const on = () => {
    if (keyboard) host.dataset.focusVisible = '';
  };
  const off = () => delete host.dataset.focusVisible;
  surface.addEventListener('focus', on);
  surface.addEventListener('blur', off);
  surface.addEventListener('pointerdown', off);
  return () => {
    surface.removeEventListener('focus', on);
    surface.removeEventListener('blur', off);
    surface.removeEventListener('pointerdown', off);
    if (--listeners > 0) return;
    document.removeEventListener('keydown', onKey, true);
    document.removeEventListener('pointerdown', onPress, true);
  };
}

/**
 * How long the pointer rests on a hint's trigger before the hint shows, in ms (seqspec:interaction#accessibility;
 * Stories' `Tooltip` too).
 */
export const HINT_DELAY_MS = 500;
/** A hint that hid less than this long ago, in ms, lets the next trigger's show at once (moving along a row). */
export const HINT_GROUP_MS = 400;
/** The gap between a hint and its trigger, and the least it keeps from the viewport's edges, in px. */
const HINT_GAP = 4;
const HINT_EDGE = 4;
/** What a hint should not cover: a control the pointer may be on its way to. */
const CONTROL =
  'button, a[href], input, select, textarea, summary, [contenteditable=""], [contenteditable="true"], [role="button"], [role="menuitem"], [role="menuitemradio"], [role="menuitemcheckbox"], [role="tab"], [role="switch"], [role="checkbox"], [role="radio"], [role="option"], [role="slider"]';
// When the last shown hint hid, shared by every hint on the page (the group window).
let hintHid = Number.NEGATIVE_INFINITY;

/**
 * A hint is a hint, never an obstacle: `bubble` (under `root`, `pointer-events: none` in the skin) shows
 * `HINT_DELAY_MS` after the pointer comes to `root`, at once within `HINT_GROUP_MS` of another hint hiding, and at once
 * on keyboard focus (`:focus-visible`); it hides at once when the pointer leaves, on a press anywhere, a scroll, a
 * wheel or a key, and when focus leaves `root`. A press or a key keeps it hidden until the pointer comes again or focus
 * moves. Never while a control under `root` is `aria-expanded` (its menu is open) or for the pointer over a menu
 * under `root`. Shown as `data-open` on `root`, on the side that covers no control (`placeHint`). Not pointer-rate:
 * the handlers run on enter, leave and focus changes.
 */
export function hintBehavior(root: HTMLElement, bubble: HTMLElement): () => void {
  const win = root.ownerDocument.defaultView;
  let timer = 0;
  let hovered = false;
  let held = false;
  let shown = false;
  const menu = (t: EventTarget | null) => t instanceof Element && !!t.closest('[role="menu"]') && root.contains(t);
  const show = () => {
    clearTimeout(timer);
    if (shown || held || !root.isConnected || root.querySelector('[aria-expanded="true"]')) return;
    shown = true;
    root.setAttribute('data-open', '');
    placeHint(root, bubble);
    win?.addEventListener('pointerdown', dismiss, true);
    win?.addEventListener('keydown', dismiss, true);
    win?.addEventListener('wheel', hide, { capture: true, passive: true });
    win?.addEventListener('scroll', hide, true);
  };
  function hide() {
    clearTimeout(timer);
    if (!shown) return;
    shown = false;
    hintHid = performance.now();
    root.removeAttribute('data-open');
    win?.removeEventListener('pointerdown', dismiss, true);
    win?.removeEventListener('keydown', dismiss, true);
    win?.removeEventListener('wheel', hide, true);
    win?.removeEventListener('scroll', hide, true);
  }
  function dismiss() {
    held = true;
    hide();
  }
  const over = (e: PointerEvent) => {
    if (e.pointerType === 'touch') return;
    if (menu(e.target)) return out();
    if (hovered) return;
    hovered = true;
    held = false;
    if (performance.now() - hintHid < HINT_GROUP_MS) show();
    else timer = win?.setTimeout(show, HINT_DELAY_MS) ?? 0;
  };
  // Leaving hides it, unless keyboard focus inside opened it and is still there.
  const out = () => {
    hovered = false;
    if (!root.querySelector(':focus-visible')) hide();
  };
  const leave = (e: PointerEvent) => {
    if (!root.contains(e.relatedTarget as Node | null) || menu(e.relatedTarget)) out();
  };
  const focusIn = (e: FocusEvent) => {
    held = false;
    if (e.target instanceof Element && e.target.matches(':focus-visible') && !menu(e.target)) show();
  };
  const focusOut = (e: FocusEvent) => {
    if (!root.contains(e.relatedTarget as Node | null) && !hovered) hide();
  };
  root.addEventListener('pointerover', over);
  root.addEventListener('pointerout', leave);
  root.addEventListener('focusin', focusIn);
  root.addEventListener('focusout', focusOut);
  return () => {
    hide();
    root.removeEventListener('pointerover', over);
    root.removeEventListener('pointerout', leave);
    root.removeEventListener('focusin', focusIn);
    root.removeEventListener('focusout', focusOut);
  };
}

/**
 * A shown hint's side: above its trigger, else below, the first that stays in the viewport and covers no control
 * outside `root` (probed with `elementFromPoint`, which passes through the hint), else the first that stays
 * in, else above; centered on the trigger and shifted to stay inside the viewport's width. Written as `data-side` on
 * `root` and the shift as `--seqv-hint-dx`; once per show, not pointer-rate.
 */
export function placeHint(root: HTMLElement, bubble: HTMLElement) {
  const doc = root.ownerDocument;
  const hit = root.getRootNode() as Document | ShadowRoot;
  const r = root.getBoundingClientRect();
  const vw = doc.documentElement.clientWidth;
  const vh = doc.documentElement.clientHeight;
  const w = bubble.offsetWidth;
  const h = bubble.offsetHeight;
  const center = r.left + r.width / 2 - w / 2;
  const x = Math.max(HINT_EDGE, Math.min(center, vw - HINT_EDGE - w));
  const fits = (y: number) => y >= HINT_EDGE && y + h <= vh - HINT_EDGE;
  const covers = (y: number) => {
    for (let py = y + 1; py < y + h; py += Math.max(1, (h - 2) / 2))
      for (let px = x + 1; px < x + w; px += Math.min(12, Math.max(1, w - 2))) {
        const el = hit.elementFromPoint(px, py);
        if (el && !root.contains(el) && el.closest(CONTROL)) return true;
      }
    return false;
  };
  const above = r.top - HINT_GAP - h;
  const below = r.bottom + HINT_GAP;
  let side: 'top' | 'bottom' = 'top';
  if (fits(above) && !covers(above)) side = 'top';
  else if (fits(below) && !covers(below)) side = 'bottom';
  else if (!fits(above) && fits(below)) side = 'bottom';
  if (root.getAttribute('data-side') !== side) root.setAttribute('data-side', side);
  root.style.setProperty('--seqv-hint-dx', `${Math.round(x - center)}px`);
}

export interface Popover {
  readonly isOpen: boolean;
  open(): void;
  close(): void;
  dispose(): void;
}
/**
 * What opens and closes a list under `root` from its `trigger` (seqspec:interaction#slots-and-chrome): `list` is put
 * under `root` while open; the trigger click, Escape, Tab, a press outside (a capture-phase `pointerdown` on the
 * document, so Safari closes too) and focus moving out of `root` are asked of `o.ask`, which the owner answers with
 * `open` or `close`
 * (a controlled menu answers through its props). `roving`: opening focuses the first item (`[role^=menuitem]`) and the
 * arrows, Home and End move between them. Closing with focus inside the list returns it to the trigger, except for a
 * press outside, whose target takes it.
 */
export function popover(
  root: HTMLElement,
  trigger: HTMLElement,
  list: HTMLElement,
  o: { roving: boolean; top?: boolean; ask(open: boolean): void },
): Popover {
  let isOpen = false;
  let back = true;
  const doc = root.ownerDocument;
  const win = doc.defaultView;
  // `top`: the list is shown in the top layer (a `manual` popover, still under `root` in the tree), so no ancestor's
  // overflow clips it, placed at the trigger by `placeMenu` while open (a scroll or a resize places it again).
  const place = () => placeMenu(trigger, list);
  const items = () => Array.from(list.querySelectorAll<HTMLElement>('[role^=menuitem]'));
  const focused = () => (root.getRootNode() as Document | ShadowRoot).activeElement;
  // `back` is the asker's to say for the close that follows; the close puts it back (an owner may answer a task later).
  const ask = (open: boolean, focus = true) => {
    back = focus;
    o.ask(open);
  };
  const press = (e: PointerEvent) => {
    if (!e.composedPath().includes(root)) ask(false, false);
  };
  const key = (e: KeyboardEvent) => {
    if (isOpen && e.key === 'Escape') {
      e.preventDefault();
      e.stopPropagation();
      return ask(false);
    }
    if (isOpen && e.key === 'Tab') return ask(false);
    const down = e.key === 'ArrowDown';
    if (!o.roving || !(down || e.key === 'ArrowUp' || e.key === 'Home' || e.key === 'End')) return;
    e.preventDefault();
    if (!isOpen) return down ? ask(true) : undefined;
    const all = items();
    const at = all.indexOf(focused() as HTMLElement);
    const next =
      e.key === 'Home'
        ? 0
        : e.key === 'End'
          ? all.length - 1
          : at < 0
            ? down
              ? 0
              : all.length - 1
            : at + (down ? 1 : -1);
    all[(next + all.length) % all.length]?.focus();
  };
  const out = (e: FocusEvent) => {
    const to = e.relatedTarget as Node | null;
    if (isOpen && to && !root.contains(to)) ask(false, false);
  };
  const click = () => ask(!isOpen);
  trigger.setAttribute('aria-expanded', 'false');
  trigger.addEventListener('click', click);
  root.addEventListener('keydown', key);
  root.addEventListener('focusout', out);
  const pop: Popover = {
    get isOpen() {
      return isOpen;
    },
    open() {
      if (isOpen) return;
      isOpen = true;
      back = true;
      root.append(list);
      if (o.top && typeof list.showPopover === 'function') {
        list.popover = 'manual';
        list.showPopover();
        place();
        win?.addEventListener('scroll', place, true);
        win?.addEventListener('resize', place);
      }
      trigger.setAttribute('aria-expanded', 'true');
      trigger.setAttribute('aria-controls', list.id);
      doc.addEventListener('pointerdown', press, true);
      if (o.roving) items()[0]?.focus();
    },
    close() {
      if (!isOpen) return;
      isOpen = false;
      const inside = list.contains(focused());
      doc.removeEventListener('pointerdown', press, true);
      win?.removeEventListener('scroll', place, true);
      win?.removeEventListener('resize', place);
      list.remove();
      trigger.setAttribute('aria-expanded', 'false');
      trigger.removeAttribute('aria-controls');
      const refocus = inside && back;
      back = true;
      if (refocus) trigger.focus();
    },
    dispose() {
      back = false;
      pop.close();
      trigger.removeEventListener('click', click);
      root.removeEventListener('keydown', key);
      root.removeEventListener('focusout', out);
    },
  };
  return pop;
}

/** The gap between a menu and its trigger, and the least a menu keeps from the viewport's edges, in px. */
const MENU_GAP = 2;
const MENU_EDGE = 4;
/**
 * A top-layer menu at its trigger (`popover`'s `top`), in viewport px as the tooltip is placed in the viewer: below the
 * trigger from its left edge; flipped above where it would leave the viewport's bottom and there is room, and to end at
 * the trigger's right edge where it would leave the right; kept inside the viewport either way.
 */
export function placeMenu(trigger: HTMLElement, list: HTMLElement) {
  const r = trigger.getBoundingClientRect();
  const view = trigger.ownerDocument.documentElement;
  const [vw, vh, w, h] = [view.clientWidth, view.clientHeight, list.offsetWidth, list.offsetHeight];
  let x = r.left;
  if (x + w > vw - MENU_EDGE) x = r.right - w;
  x = Math.max(MENU_EDGE, Math.min(x, vw - MENU_EDGE - w));
  let y = r.bottom + MENU_GAP;
  if (y + h > vh - MENU_EDGE && r.top - MENU_GAP - h >= MENU_EDGE) y = r.top - MENU_GAP - h;
  setStyle(list, 'left', `${Math.round(x)}px`);
  setStyle(list, 'top', `${Math.round(y)}px`);
}

/** A style written only when it changes (an unchanged write is still a mutation). */
export function setStyle(
  el: HTMLElement,
  prop: 'transform' | 'visibility' | 'left' | 'top' | 'height' | 'width',
  v: string,
) {
  if (el.style[prop] !== v) el.style[prop] = v;
}

/** What `placeTooltip` last wrote, so it writes only a change; one per tooltip. */
export interface TipSpot {
  x: number;
  y: number;
  shown: boolean;
}
export const tipSpot = (): TipSpot => ({ x: Number.NaN, y: Number.NaN, shown: false });
const RECT = new Float64Array(6);
/** The hovered header's line, for the tooltip of a track header's legend (`placeTooltip`); one, reused. */
const LINE = headerLine();
/** "12px", made once for a whole number in [0, 8192) (a place on a viewer is a few thousand values at most). */
const PX: string[] = [];
function px(n: number): string {
  const i = n | 0;
  if (i !== n || i < 0 || i >= 8192) return `${n}px`;
  PX[i] ??= `${i}px`;
  return PX[i];
}

/**
 * The tooltip at its owner's anchor: beside the pointer for the pointer's hover, at the keyboard cursor's cell or
 * header (`cursorRect`) for the keyboard's; hidden otherwise. It goes right of and below the anchor, and flips to the
 * left and above where it would leave the viewer (its offset parent). Measured from `offsetLeft` and `offsetTop` up to
 * that parent and its `clientWidth` and `clientHeight`, read again at each call, so a part that grew above the sections
 * moves the tooltip with the surface; no rects, no arrays; `at` holds what was written, and a call that changes nothing
 * writes nothing.
 */
export function placeTooltip(
  ctx: SeqViewerContext,
  surfaces: ReadonlyMap<string, HTMLElement>,
  t: HTMLElement,
  at: TipSpot,
) {
  const owner = ctx.hover.live.current.owner;
  const p = ctx.pointer.current;
  const parent = t.offsetParent as HTMLElement | null;
  const kb = owner === 'keyboard';
  let surface: HTMLElement | undefined;
  let ax = 0;
  let ay = 0;
  let aw = 0;
  let ah = 0;
  let line = false;
  if (owner === 'pointer' && p.inside) {
    const ref = ctx.layout.sectionRef(p.section);
    surface = ref === null ? undefined : surfaces.get(ref);
    ax = p.x;
    ay = p.y;
    // A track header's legend: below (or above) the header's whole line, never over its menu button on that line.
    const h = ctx.hover.live.current;
    const view = ref !== null && h.kind === HoverKind.header && h.track >= 0 ? ctx.graph.assembled.value?.view : null;
    if (ref !== null && view) {
      const track = view.set.strings.get(view.resolved.tracks.ref[h.track]);
      if (trackLine(ctx, ref, track, LINE)) {
        ay = LINE.top;
        ah = LINE.height;
        line = true;
      }
    }
  } else if (kb) {
    const c = ctx.ui.keyCursor.value;
    if (c && cursorRect(ctx, c, RECT)) {
      surface = surfaces.get(c.section);
      ax = RECT[0];
      ay = RECT[1];
      aw = RECT[2];
      ah = RECT[3];
    }
  }
  let ox = 0;
  let oy = 0;
  let e: Element | null = surface ?? null;
  while (e instanceof HTMLElement && e !== parent) {
    ox += e.offsetLeft;
    oy += e.offsetTop;
    e = e.offsetParent;
  }
  const shown = !!surface && !!parent && e === parent;
  if (shown !== at.shown) {
    at.shown = shown;
    setStyle(t, 'visibility', shown ? '' : 'hidden');
  }
  if (!shown || !parent) return;
  const tw = t.offsetWidth;
  const th = t.offsetHeight;
  let x = ox + ax + aw + (kb ? 6 : 14);
  if (x + tw > parent.clientWidth) x = ox + ax - (kb ? 6 : 12) - tw;
  const gap = kb || line ? 4 : 18;
  let y = oy + ay + ah + gap;
  x = Math.max(0, Math.round(x));
  // Where it would leave the viewer below: a header's legend goes above its line; the pointer's or the cursor's label
  // stays beside the anchor, as low as the viewer holds it, so it covers no ruler above the tracks, where it is beside
  // it (clear of it across); a viewer too narrow for that puts it above the anchor, off what it labels. A viewer too
  // short for it takes it from the top.
  if (y + th > parent.clientHeight) {
    const beside = x + tw <= ox + ax || x >= ox + ax + aw;
    y = !line && beside ? parent.clientHeight - th : oy + ay - (kb || line ? 4 : 12) - th;
  }
  y = Math.max(0, Math.round(y));
  if (x === at.x && y === at.y) return;
  at.x = x;
  at.y = y;
  // Two custom properties of cached strings (skin.css moves the tooltip by `transform` from them): nothing is allocated.
  t.style.setProperty('--seqv-x', px(x));
  t.style.setProperty('--seqv-y', px(y));
}

/**
 * A track header's menu at the header's right edge, on the track's line (`trackLine`: a rows track's hovered row, else
 * its entry; `line` is the caller's, kept across calls so the menu stays put while the pointer is over it); hidden
 * while the line is out of view. Called on every hover change, not only when the track changes.
 */
export function placeTrackActions(
  ctx: SeqViewerContext,
  section: string,
  track: string,
  t: HTMLElement,
  line: HeaderLine,
) {
  const l = ctx.layout.get(section)?.current;
  if (!l) return;
  const shown = trackLine(ctx, section, track, line);
  setStyle(t, 'visibility', shown ? '' : 'hidden');
  // On every hover change and viewport frame while shown: whole numbers come from the cache, nothing is allocated.
  setStyle(t, 'left', px(l.headerWidth - ctx.layout.props(section).controlWidth));
  if (!shown) return;
  setStyle(t, 'top', px(line.top));
  setStyle(t, 'height', px(line.height));
}

/** What `trackBox` asks of its UI: the box now shown (null: none), and the change of track to show. */
export interface TrackBoxOwner {
  box(): HTMLElement | null;
  /** The track whose box is to be shown, null: none. Called when it changes; the owner builds or removes the box. */
  show(track: string | null): void;
}
export interface TrackBox {
  /** Decides the track again (a hover or cursor change, the layout's); `now`: a removal that was waiting may go. */
  update(now?: boolean): void;
  /** Puts the box on the track's line (the viewport moved, or the box is new). */
  place(): void;
  /** The pointer entered (the track) or left (null) the box. */
  hold(track: string | null): void;
  /** A removal that waits a task, for focus or the pointer to settle (the box's `focusout`). */
  later(): void;
  /** The layout row of the line the box is on (a split track's value row, ≥ 0; else −1 or −2). */
  row(): number;
  dispose(): void;
}
/**
 * Which track's box `TrackActions` is shown for (seqspec:interaction#the-default-ui), for both UIs: the track under the
 * pointer's header, else the keyboard cursor's; kept while the pointer is over the box, while a control in it is
 * `aria-expanded` (so a host's Menu counts) and, with nothing hovered, while focus is inside it and the track still has
 * a line on screen. Leaving the surface for the box clears the hover before the box's `pointerenter` holds it, so a
 * removal waits a task (React batches the two; DOM parts write at once). Pointer-rate: the same track is only placed
 * again.
 */
export function trackBox(ctx: SeqViewerContext, section: string, o: TrackBoxOwner): TrackBox {
  let held: string | null = null;
  let track: string | null = null;
  let timer = 0;
  const line = headerLine();
  const place = () => {
    const box = o.box();
    if (box && track) placeTrackActions(ctx, section, track, box, line);
  };
  const settle = () => {
    timer = 0;
    update(true);
  };
  const update = (now = false) => {
    const box = o.box();
    const open = !!box?.querySelector('[aria-expanded=true]');
    // Focus inside keeps the box only while its track still has a line on screen (a track hidden by its own menu does not).
    const inside = !!box?.matches(':focus-within') && !!track && trackLine(ctx, section, track, line);
    const next = open ? track : (held ?? headerTrack(ctx, section) ?? (inside ? track : null));
    if (next === track) return place();
    if (!next && box && !now && !box.matches(':focus-within')) {
      timer ||= window.setTimeout(settle, 0);
      return;
    }
    track = next;
    line.entry = -1;
    o.show(next);
  };
  return {
    update,
    place,
    hold(t) {
      held = t;
      update();
    },
    later() {
      timer ||= window.setTimeout(settle, 0);
    },
    row: () => (line.entry >= 0 ? (ctx.layout.get(section)?.current.entryRow[line.entry] ?? -1) : -1),
    dispose: () => window.clearTimeout(timer),
  };
}

// Quarter pixels within a 1,024 px page, and the pages: two tables a scroll fills once per value it reaches (4,096 strings
// at most, and one per page), so a scroll anywhere, a tall alignment's far rows too, allocates nothing once warm.
const SCROLL_PAGE = 4096;
const SCROLL_LOW: string[] = [];
const SCROLL_PAGES = new Map<number, string>();
/**
 * A scrolling box's offset, in quarter pixels, as two custom properties skin.css adds in its `transform`:
 * `--seqv-scroll-page` (whole 1,024 px pages) and `--seqv-scroll-y` (the rest), each written only when it changes, from
 * cached strings (a `translateY()` string per frame cost 70 bytes, like the tooltip's). `at` is the caller's last place.
 */
export function placeScroll(el: HTMLElement, y: number, at: { quarter: number; page: number }) {
  const k = Math.round(y * 4);
  if (k === at.quarter) return;
  at.quarter = k;
  const page = Math.floor(k / SCROLL_PAGE);
  if (page !== at.page) {
    at.page = page;
    let v = SCROLL_PAGES.get(page);
    if (v === undefined) {
      v = `${(-page * SCROLL_PAGE) / 4}px`;
      SCROLL_PAGES.set(page, v);
    }
    el.style.setProperty('--seqv-scroll-page', v);
  }
  const low = k - page * SCROLL_PAGE;
  SCROLL_LOW[low] ??= `${-low / 4}px`;
  el.style.setProperty('--seqv-scroll-y', SCROLL_LOW[low]);
}

/** How long `Busy` stays after the last region in flight lands: a stream's next region usually starts within it. */
export const BUSY_HOLD_MS = 250;
/**
 * Whether `Busy` shows, for both UIs, told on change: at once while tasks run or regions in view are in flight; off
 * when neither, but `BUSY_HOLD_MS` after the last region landed, so a stream fetching region after region (a genome
 * panned) does not flicker it.
 */
export function busyShown(ctx: SeqViewerContext, tell: (busy: boolean) => void): { dispose(): void } {
  let shown: boolean | null = null;
  let held = false;
  let timer = 0;
  const release = () => {
    timer = 0;
    held = false;
    update();
  };
  const update = () => {
    if (ctx.streams.stats.current.inFlight > 0) {
      held = true;
      window.clearTimeout(timer);
      timer = 0;
    } else if (held) timer ||= window.setTimeout(release, BUSY_HOLD_MS);
    const busy = held || ctx.tasks.busy.value;
    if (busy === shown) return;
    shown = busy;
    tell(busy);
  };
  const subs = watch([ctx.tasks.busy, ctx.streams.stats.version], update);
  return {
    dispose() {
      for (const s of subs) s.unsubscribe();
      window.clearTimeout(timer);
    },
  };
}

/** The `Search` part's "Keep" icon, a pin (16 units): its outline and its needle. */
export const PIN = {
  outline: 'M9.8 1.8 14.2 6.2 12.5 6.9 10.2 9.2 10.4 11.4 9.4 12.4 3.6 6.6 4.6 5.6 6.8 5.8 9.1 3.5Z',
  needle: 'M6.4 9.6 2.2 13.8',
} as const;

/**
 * A control that hides while it has focus (Keep, Clear kept: `data-idle`) hands focus to `to` first, so focus is never
 * lost to the page (seqspec:interaction#accessibility).
 */
export function handOff(control: HTMLElement | null, idle: boolean, to: HTMLElement | null): void {
  if (idle && control && (control.getRootNode() as Document | ShadowRoot).activeElement === control) to?.focus();
}

/** The `Search` part's "Clear kept" icon (16 units): a cross. */
export const CROSS = 'M4.5 4.5 11.5 11.5M11.5 4.5 4.5 11.5';

/** The legend title's palette icon (16 units): its outline, and its three paint dots (cx, cy). */
export const PALETTE = {
  outline:
    'M8 1.5a6.5 6.5 0 1 0 0 13c.9 0 1.4-.6 1.4-1.3 0-.4-.2-.7-.4-1s-.4-.6-.4-1c0-.7.6-1.3 1.3-1.3h1.6a3 3 0 0 0 3-3C14.5 3.9 11.6 1.5 8 1.5Z',
  dots: [
    [4.8, 7.6],
    [6.4, 4.6],
    [9.8, 4.6],
  ],
} as const;

/** A legend ramp as a CSS gradient (two stops at one `t`, a discrete palette's edge, a hard step). */
export const legendGradient = (r: LegendRamp) =>
  `linear-gradient(to right, ${r.stops.map((s) => `${s.color} ${Math.round(s.t * 1000) / 10}%`).join(', ')})`;
/** The slot a band's share is drawn in, CSS px (a sample bar `value` of it tall). */
const SIZE_SLOT = 14;
/**
 * A legend sample size drawn as its channel draws it, as CSS for its swatch: a stroke `width` a line of that thickness;
 * a band's share (`thickness`, `height`) a bar that share of a 14 px slot; a `radius` a dot twice as wide (else, a
 * renderer's own size, a dot `value` px across). In the sample's color (a merged ramp's), else the ink.
 */
export function legendSizeStyle(s: LegendSize): Record<'width' | 'height' | 'borderRadius' | 'background', string> {
  const clamp = (x: number, lo: number, hi: number) => `${Math.min(hi, Math.max(lo, Math.round(x)))}px`;
  const background = s.color ?? 'var(--seq-ink)';
  if (s.unit === SizeUnit.band)
    return { width: '8px', height: clamp(s.value * SIZE_SLOT, 1, SIZE_SLOT), borderRadius: '0', background };
  if (s.channel === 'width') return { width: '18px', height: clamp(s.value, 1, 12), borderRadius: '1px', background };
  const d = clamp(s.channel === 'radius' ? 2 * s.value : s.value, 3, 24);
  return { width: d, height: d, borderRadius: '50%', background };
}
/** The caption before sample `i`: its field, where a run of one scale's samples starts (null: none). */
export const sizeCaption = (sizes: readonly LegendSize[], i: number): string | null =>
  sizes[i].field !== null &&
  (i === 0 || sizes[i - 1].field !== sizes[i].field || sizes[i - 1].channel !== sizes[i].channel)
    ? sizes[i].field
    : null;

/** A custom property written only when it changes. */
function setVar(el: HTMLElement, name: string, v: string) {
  if (el.style.getPropertyValue(name) !== v) el.style.setProperty(name, v);
}
/** The narrowest label column a legend group's name sits in: the label inset, its padding and a few letters. */
const LEGEND_LABEL_MIN = 40;
/**
 * The width a group's title needs in the label column: its name's whole text, the cell's padding, and a picker's icon
 * and padding (the button's width beyond the name). 0: no name, or none laid out (a closed legend).
 */
function titleWidth(el: HTMLElement): number {
  const name = el.querySelector<HTMLElement>('.seqv-legend-name');
  const cell = name?.closest<HTMLElement>('.seqv-legend-track');
  if (!name || !cell || name.scrollWidth === 0) return 0;
  const cs = getComputedStyle(cell);
  let w = name.scrollWidth + Number.parseFloat(cs.paddingLeft) + Number.parseFloat(cs.paddingRight);
  const b = name.parentElement;
  if (b && b !== cell) w += b.offsetWidth - name.offsetWidth;
  return Math.ceil(w);
}
/**
 * The legend on the canvas's grid (seqspec:interaction#the-legend, `Legend`): on `el` (a group, or the legend for its
 * toggle), its section's data left and width, row height and label font as `--seqv-*` properties, read from the
 * section's layout as the
 * canvas reads them, so a group's name sits in the label column and its items start at the data area. A wrapped
 * section has no label column (its gutter holds line numbers), nor has one narrower than `LEGEND_LABEL_MIN`:
 * `data-seqv-stacked` puts the name on its own line above the items. A name the label column cut ("Transcripts (10)"
 * past the tracks' own labels) gets its own line at the label inset (`data-seqv-long`), the items under it at the data
 * area's left, so a title is cut only past the legend's width (`data-seqv-cut`, and `onCut`: the part gives the title
 * a hint with the full name).
 * Follows the layout (a resize, the label column's width) and the title's size
 * (a legend opened, a new name); writes only on change.
 */
export function legendGrid(
  ctx: SeqViewerContext,
  section: string,
  el: HTMLElement,
  onCut?: (cut: boolean) => void,
): () => void {
  let follow: Unsubscribable | undefined;
  const write = () => {
    const l = ctx.layout.get(section)?.current;
    if (!l) return;
    const p = ctx.layout.props(section);
    setVar(el, '--seqv-data-left', `${l.dataLeft}px`);
    setVar(el, '--seqv-data-width', `${Math.max(0, l.dataRight - l.dataLeft)}px`);
    setVar(el, '--seqv-row', `${p.rowHeight}px`);
    setVar(el, '--seqv-label-px', `${p.labelFontPx}px`);
    const stacked = l.geometry === AxisGeometry.wrapped || l.dataLeft < LEGEND_LABEL_MIN;
    if (el.hasAttribute('data-seqv-stacked') !== stacked) el.toggleAttribute('data-seqv-stacked', stacked);
    if (!el.classList.contains('seqv-legend-group')) return;
    const long = !stacked && titleWidth(el) > l.dataLeft;
    if (el.hasAttribute('data-seqv-long') !== long) el.toggleAttribute('data-seqv-long', long);
    // Cut even so (a name wider than the whole legend): its hint shows the full name.
    const name = el.querySelector<HTMLElement>('.seqv-legend-name');
    const cut = !!name && name.scrollWidth > name.clientWidth;
    if (el.hasAttribute('data-seqv-cut') === cut) return;
    el.toggleAttribute('data-seqv-cut', cut);
    onCut?.(cut);
  };
  // A title's size changes with no layout change (the legend opened, its name rendered anew): measured again.
  const sized = typeof ResizeObserver === 'function' ? new ResizeObserver(() => write()) : null;
  sized?.observe(el);
  // The section's layout can come after the legend (the view assembles first): watched once the section is there.
  const sections = ctx.viewports.sections.subscribe(() => {
    follow ??= ctx.layout.get(section)?.version.subscribe(write);
  });
  return () => {
    sections.unsubscribe();
    follow?.unsubscribe();
    sized?.disconnect();
  };
}
