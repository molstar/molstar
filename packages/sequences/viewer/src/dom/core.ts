/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The vanilla UI's contracts and helpers (seqspec:interaction#the-dom-element): a part is a `DomPart`, a function that
// mounts into a host and returns its `dispose`; slots replace parts with other `DomPart`s, and a host mounts any part
// into a host element of its own (with a mounted UI's `ui`, or none: the default primitives). No virtual DOM: a part keeps
// the nodes it made and writes only those its state owns, an attribute or a text only when it changed. The primitives a
// host maps to its own (`Hint`, `ScrollArea`, `Menu`) are parts too. A node a part swaps keeps keyboard focus by its
// `data-seqv-key` (`spot`, seqspec:interaction#accessibility).

import type { MenuItem } from "../chrome.ts";
import type { LegendMode, SearchFn, SeqState, SeqViewerContext, Unsubscribable } from "../index.ts";
import { sameShallow } from "../same.ts";
import { hintBehavior, later, popover, watch } from "./shared.ts";

export { later, watch };

export interface DomHandle<P> {
  update?(props: P): void;
  /** Unsubscribes every state and removes every node the part added; the context outlives it. */
  dispose(): void;
}
export type DomPart<P = object> = (host: HTMLElement, context: SeqViewerContext, props: P) => DomHandle<P>;

/** What a surface factory returns (the canvas engine's `mountSection` by default; tests pass a stub). */
export interface MountedSurface {
  /**
   * The surface's element, its accessibility its own: one tab stop, `role="group"`, `aria-roledescription` and the
   * section's label, as `SectionSurface` sets them (the UIs do not; a host's factory sets them too).
   */
  readonly element: HTMLElement;
  dispose(): void;
}
export type SurfaceFactory = (context: SeqViewerContext, section: string, host: HTMLElement) => MountedSurface;

export interface PartProps {
  /** The UI it is part of (`mountSeqViewer`'s): its slots and surfaces; absent, the default primitives. */
  readonly ui?: DomUI;
}
export interface SectionProps extends PartProps {
  readonly section: string;
}
export interface TrackProps extends SectionProps {
  readonly track: string;
  /** The header line's layout row, read when the menu is built (a split row's index; −1: the track's own line). */
  readonly row?: () => number;
}
export interface RowProps extends SectionProps {
  readonly displayRow: number;
  readonly dataRow: number;
}
/** `children` fills the element the hint is for (React's children); the control keeps its own accessible name. */
export interface HintProps {
  readonly label: string;
  readonly children: (host: HTMLElement) => void;
}
export interface ScrollAreaProps {
  readonly className?: string;
  readonly children: (viewport: HTMLElement) => void;
}
/** Controlled, in the shape of Base UI's: the part that owns the menu holds `open` and passes it again through `update`. */
export interface MenuProps {
  /** The button that opens it; the menu adds its own handlers and ARIA state. */
  readonly trigger: HTMLButtonElement;
  readonly items: readonly MenuItem[];
  readonly open: boolean;
  /** Asked for on the trigger, Escape, Tab, an item chosen and a press outside; the owner answers with `update`. */
  onOpenChange(open: boolean): void;
}

/** seqspec:interaction#slots-and-chrome: every part, and the primitives they are drawn with. */
export interface DomSlots {
  ViewTabs: DomPart<PartProps>;
  /** A section's label, its filter chip, its `HiddenChips` and its `SectionSwitches`. */
  SectionHeader: DomPart<SectionProps>;
  /** A section's switches, each an action. */
  SectionSwitches: DomPart<SectionProps>;
  TrackActions: DomPart<TrackProps>;
  HiddenChips: DomPart<SectionProps>;
  RowControl: DomPart<RowProps>;
  Tooltip: DomPart<PartProps>;
  /** Item actions on hover (seqspec:interaction#the-hover-card); inert without `UI.ItemActions`. */
  HoverCard: DomPart<PartProps>;
  Legend: DomPart<PartProps>;
  SnapshotPlayer: DomPart<PartProps>;
  Description: DomPart<PartProps>;
  /** Only when the UI spec gives `search`. */
  Search: DomPart<PartProps>;
  Export: DomPart<PartProps>;
  Issues: DomPart<PartProps>;
  Busy: DomPart<PartProps>;
  EmptyState: DomPart<PartProps>;
  ErrorState: DomPart<PartProps>;
  Hint: DomPart<HintProps>;
  ScrollArea: DomPart<ScrollAreaProps>;
  Menu: DomPart<MenuProps>;
}
export type DomSlotName = keyof DomSlots;
export type Primitive = "Hint" | "ScrollArea" | "Menu";

export interface DomUISpec {
  /**
   * A part, "none" (left out wherever it is placed; a primitive keeps its default), or "default" (its default, also
   * placed under `clean` chrome).
   */
  readonly slots?: { readonly [K in DomSlotName]?: DomSlots[K] | "none" | "default" };
  /** On the root, where a host scopes its token mapping. */
  readonly className?: string;
  /** The parts `mountSeqViewer` places (`ChromeCompositions`); unset: the context's `UI.Chrome`. Never written to it. */
  readonly chrome?: "full" | "clean";
  /**
   * This UI's legend (seqspec:interaction#the-legend): collapsed, expanded or hidden; unset: the context's `UI.Legend`.
   */
  readonly legend?: LegendMode;
  readonly surface?: SurfaceFactory;
  /** An application's search; absent: no `Search` part. */
  readonly search?: SearchFn;
}
/** What every part receives as `props.ui`. */
export interface DomUI {
  readonly slots: { readonly [K in Exclude<DomSlotName, Primitive>]: DomSlots[K] | null } & Pick<DomSlots, Primitive>;
  /** The parts the default layout places (`placedParts`): the chrome's composition and the host's slots. */
  readonly placed: ReadonlySet<string>;
  readonly className: string;
  /** The spec's `legend`, if given. */
  readonly legend?: LegendMode;
  readonly surface: SurfaceFactory;
  readonly search?: SearchFn;
  /** Mounted surface elements by section ref (the tooltip places itself from them). */
  readonly surfaces: Map<string, HTMLElement>;
}

export const cx = (...names: (string | false | null | undefined)[]) => names.filter(Boolean).join(" ");

export type Attr = string | number | boolean | null | undefined;
/** As React writes them: `aria-*` booleans as "true"/"false"; other booleans present or absent; null/undefined absent. */
export function setAttr(e: Element, name: string, v: Attr) {
  const aria = name.startsWith("aria-");
  const s = v == null || (v === false && !aria) ? null : typeof v === "boolean" ? (aria ? String(v) : "") : String(v);
  if (s === null) e.removeAttribute(name);
  else if (e.getAttribute(name) !== s) e.setAttribute(name, s);
}
export function setText(e: Node, s: string) {
  if (e.textContent !== s) e.textContent = s;
}

export type Child = Node | string | false | null | undefined;
export function h<K extends keyof HTMLElementTagNameMap>(
  tag: K,
  attrs?: Record<string, Attr> | null,
  ...children: Child[]
): HTMLElementTagNameMap[K] {
  const e = document.createElement(tag);
  for (const k in attrs) setAttr(e, k, attrs[k]);
  for (const c of children) if (c) e.append(c);
  return e;
}
/** A button with its click handler. */
export function button(attrs: Record<string, Attr>, onClick: (e: MouseEvent) => void, ...children: Child[]) {
  const b = h("button", { type: "button", ...attrs }, ...children);
  b.addEventListener("click", onClick);
  return b;
}

// Keyboard focus is never lost. Every node a part may rebuild carries a `data-seqv-key`, unique in its part, and the
// nodes of a row (the chips) a `data-seqv-group`. What changes a node reads the focused key before and focuses the same
// key after, synchronously, in the same task.
const activeIn = (n: Node) => (n.getRootNode() as Document | ShadowRoot).activeElement as HTMLElement | null;
const keyed = (k: string) => `[data-seqv-key="${CSS.escape(k)}"]`;
interface Held {
  readonly key: string;
  readonly group: string | null;
  /** The keys of the group in the node being replaced, in order. */
  readonly order: readonly string[];
}
/** The focused key inside `scope`: "" for focus on an element without one; null: focus is elsewhere. */
function hold(scope: Element): Held | null {
  const a = activeIn(scope);
  if (!a || !scope.contains(a)) return null;
  const at = a.closest<HTMLElement>("[data-seqv-key]");
  const group = at && scope.contains(at) ? (at.dataset.seqvGroup ?? null) : null;
  const order = group
    ? Array.from(
        scope.querySelectorAll<HTMLElement>(`[data-seqv-group="${CSS.escape(group)}"]`),
        (e) => e.dataset.seqvKey ?? "",
      )
    : [];
  return { key: at && scope.contains(at) ? (at.dataset.seqvKey ?? "") : "", group, order };
}
/** The section's surface, else the viewer's root: where focus lands when nothing of its own remains. */
const landing = (near: Node) => {
  const e = near as Element;
  return (
    e.closest?.(".seqv-section")?.querySelector<HTMLElement>(".seqv-surface [tabindex]") ??
    e.closest?.<HTMLElement>(".seqv-root") ??
    null
  );
};
/** Focus after a change: the same key if it exists in `scope`, else the next of its group, else the previous, else the landing. */
function land(scope: Element | null, was: Held, near: Node) {
  const find = (k: string): HTMLElement | null =>
    scope?.matches(keyed(k)) ? (scope as HTMLElement) : (scope?.querySelector<HTMLElement>(keyed(k)) ?? null);
  let t = was.key ? find(was.key) : null;
  if (!t && was.group) {
    const i = was.order.indexOf(was.key);
    const nearest = [...was.order.slice(i + 1), ...was.order.slice(0, Math.max(i, 0)).reverse()];
    t = nearest.map(find).find((e) => e !== null) ?? null;
  }
  (t ?? landing(near))?.focus();
}
/** Runs `change`, which rewrites the children of `scope`, then refocuses the key that held focus there. */
export function keepFocus(scope: Element, change: () => void) {
  const was = hold(scope);
  change();
  if (was) land(scope, was, scope);
}
/** Removes `node` (through `remove`); focus that was inside it lands on `near`'s section surface or the viewer's root. */
export function lose(node: Element, near: Node, remove: () => void = () => node.remove()) {
  const had = node.contains(activeIn(near));
  remove();
  if (had) landing(near)?.focus();
}

// The nodes `landFocus` took focus from in this task: a framework removing several nodes at once runs each one's cleanup
// before it removes any, so a sibling among them is no landing (it would pass focus back, then go with it).
const leaving = new Set<Element>();
/**
 * `el` is going away (a framework unmounts it, where `lose` removes the node itself): focus inside it lands on the next
 * node of its `data-seqv-group` after it, else the previous, else the surface of `near`'s section (default `el`'s), else
 * the viewer's root. Call it before the node is removed.
 */
export function landFocus(el: HTMLElement, near: Node = el) {
  if (!el.contains(activeIn(el))) return;
  if (leaving.size === 0) queueMicrotask(() => leaving.clear());
  leaving.add(el);
  const g = el.dataset.seqvGroup;
  const sibling = g
    ? ([el.nextElementSibling, el.previousElementSibling] as (HTMLElement | null)[]).find(
        (n) => n?.dataset.seqvGroup === g && !leaving.has(n),
      )
    : null;
  (sibling ?? landing(near))?.focus();
}

/** A part's place in its host, held by an empty comment: the part's root comes and goes there, in order. */
export interface Spot {
  /** Swaps the node. Focus inside the one replaced stays on the same `data-seqv-key` of the new one, else lands. */
  set(node: ChildNode | null): void;
  dispose(): void;
}
export function spot(host: Node): Spot {
  const mark = host.appendChild(document.createComment(""));
  let shown: ChildNode | null = null;
  return {
    set(n) {
      if (n === shown) return;
      const was = shown instanceof Element ? hold(shown) : null;
      shown?.remove();
      if (n) mark.before(n);
      shown = n;
      if (was) land(n instanceof Element ? n : null, was, host);
    },
    dispose() {
      const node = shown;
      if (node instanceof Element) lose(node, mark.parentNode ?? host, () => node.remove());
      else shown?.remove();
      shown = null;
      mark.remove();
    },
  };
}

export type Subscribable = Pick<SeqState<unknown>, "subscribe">;
/**
 * `apply` with a value derived from states, run when it changed by `compare` (shallow by default; `sameData` for a
 * list of records, or a caller's own): React's `useDerived`, so a part writes only when what it shows changed.
 * `compute` runs once at mount, then once per change.
 */
export function derived<T>(
  states: readonly Subscribable[],
  compute: () => T,
  apply: (v: T) => void,
  compare: (a: T, b: T) => boolean = sameShallow,
): Unsubscribable[] {
  let last: { value: T } | undefined;
  return watch(states, () => {
    const value = compute();
    if (last && compare(last.value, value)) return;
    last = { value };
    apply(value);
  });
}
/** Unsubscribes and disposes, in order. */
export function off(subs: readonly (Unsubscribable | { dispose(): void } | undefined | null)[]) {
  for (const s of subs)
    if (s && "unsubscribe" in s) s.unsubscribe();
    else s?.dispose();
}

/** An accessible hint on hover and on keyboard focus, never in the pointer's way (`hintBehavior`). */
export const Hint: DomPart<HintProps> = (host, _ctx, props) => {
  const root = h("span", { class: "seqv-hint" });
  props.children(root);
  const bubble = h("span", { class: "seqv-hint-bubble", "aria-hidden": "true" }, props.label);
  root.append(bubble);
  host.append(root);
  const behavior = hintBehavior(root, bubble);
  return {
    update: (p) => setText(bubble, p.label),
    dispose: () => {
      behavior();
      root.remove();
    },
  };
};

/** A plain region: the default UI scrolls nothing in the DOM (a host maps this to its own scroll area). */
export const ScrollArea: DomPart<ScrollAreaProps> = (host, _ctx, props) => {
  const root = h("div", { class: cx("seqv-scroll", props.className) });
  props.children(root);
  host.append(root);
  return { dispose: () => root.remove() };
};

let menus = 0;
/**
 * The trigger, then a list of buttons while `open`: opening focuses the first item, the arrows move between items,
 * Escape, Tab, choosing an item and a press outside ask the owner to close it, and closing with focus inside returns it
 * to the trigger (`popover`).
 */
export const Menu: DomPart<MenuProps> = (host, _ctx, props) => {
  const { trigger: t } = props;
  let p = props;
  const list = h("span", { id: `seqv-menu-${++menus}`, role: "menu", class: "seqv-menu-list" });
  const root = h("span", { class: "seqv-menu" }, t);
  const pop = popover(root, t, list, { roving: true, top: true, ask: (open) => p.onOpenChange(open) });
  const item = (i: MenuItem) =>
    button(
      {
        role: i.checked === undefined ? "menuitem" : "menuitemradio",
        "aria-checked": i.checked,
        class: cx("seqv-menu-item", i.danger && "seqv-danger"),
        "data-seqv-key": `item:${i.key}`,
        tabindex: -1,
      },
      () => {
        p.onOpenChange(false);
        i.run();
      },
      i.label,
    );
  const fill = () => list.replaceChildren(...p.items.map(item));
  const sync = () => {
    if (p.open !== pop.isOpen) {
      if (!p.open) return pop.close();
      fill();
      return pop.open();
    }
    if (p.open) keepFocus(list, fill);
  };
  setAttr(t, "aria-haspopup", "menu");
  host.append(root);
  sync();
  return {
    update(next) {
      p = next;
      sync();
    },
    dispose() {
      pop.dispose();
      root.remove();
    },
  };
};

/** The primitives a part draws with: its UI's slots, else the defaults (a part a host placed without a UI). */
export const primitives = (ui: DomUI | undefined): Pick<DomSlots, Primitive> => ui?.slots ?? { Hint, ScrollArea, Menu };

/** A `Menu` slot whose open state its owner (a header's switch) keeps itself, in this closure. */
export function ownedMenu(
  M: DomSlots["Menu"],
  host: HTMLElement,
  ctx: SeqViewerContext,
  trigger: HTMLButtonElement,
  items: readonly MenuItem[],
) {
  let open = false;
  const props = (): MenuProps => ({ trigger, items, open, onOpenChange });
  const onOpenChange = (o: boolean) => {
    open = o;
    menu.update?.(props());
  };
  const menu = M(host, ctx, props());
  return menu;
}
