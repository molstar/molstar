/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Observable } from 'rxjs';
import type { IconName } from '../components/icon/icon.js';
import type { Child } from '../dom/jsx-types.js';
import type { ViewRenderer } from '../dom/view.js';
import { type ObservableRecordState, type ObservableState, observableState } from '../reactive/observable-state.js';
import { ReactiveModel } from '../reactive/reactive-model.js';
import type { ReadonlyState } from '../reactive/readonly-state.js';

export type DockId = 'left' | 'right' | 'top' | 'bottom';
export const DockIds: readonly DockId[] = ['left', 'right', 'top', 'bottom'];

export type ShellMode = 'desktop' | 'compact' | 'embedded';

/** What a tab, overlay or canvas item is rendered with. */
export interface ShellViewContext {
  shell: ShellModel;
  /** The dock a tab is shown in (tabs only). */
  dock?: DockId;
  /** The section to bring into view (overlays opened at a section). */
  focus?: string;
  /** What a section is about (the inspected scene node for `inspector:<kind>` sections). */
  subject?: unknown;
}

/** A panel tab, contributed by the default UI or an extension. */
export interface ShellTab {
  id: string;
  /** Its name: the rail's tooltip, and the panel's header unless `heading` says otherwise. */
  title: string;
  /** The panel's header when it differs from the title: text, or a view (Home's is the Mol* wordmark). */
  heading?: string | (() => Child);
  icon: IconName;
  /** Where it goes by default ('right'). */
  dock?: DockId;
  order?: number;
  /**
   * Where its button sits on a floating dock's rail: its head (Home, as the Mol* mark, a section of its own), the tab
   * strip (default), or its foot by the command menu (Settings).
   */
  place?: 'head' | 'strip' | 'foot';
  /**
   * Which rail buttons go first under More when the rail is too short for all of them, even at its smallest size: the
   * lower the number, the longer the button stays. Defaults to its position in the rail (top to bottom), so the foot goes
   * first; the command menu says -1 (it is the main way in).
   */
  priority?: number;
  /** In the compact layout's bottom bar (Home, Structure); the rest go under More there. */
  compactBar?: boolean;
  /**
   * The look of its content on a floating top or bottom dock: `card` (default) is the glass card, `none` has no card at
   * all (no fill, shadow, filter or hairline) and the canvas shows through, for a view that draws on the canvas itself
   * (the sequence). Ignored in a side dock and in the docked (not floating) layout.
   */
  frame?: 'card' | 'none';
  render: ViewRenderer<ShellViewContext>;
}

/** A card over the canvas, opened by a toolbar button or a command. */
export interface ShellOverlay {
  id: string;
  title: string;
  render: ViewRenderer<ShellViewContext>;
}

/**
 * A page over the whole app, for what needs the room a panel or a card does not have (the agent's Inspector): opened by
 * a command or a link, closed with Esc or its ✕, the app unchanged underneath.
 */
export interface ShellPage {
  id: string;
  title: string;
  icon?: IconName;
  /** What goes in the header after the title (the page's own tabs); its content takes the rest. */
  header?: ViewRenderer<ShellViewContext>;
  render: ViewRenderer<ShellViewContext>;
}

export type DockState = {
  /** Tab ids in this dock, in order. */
  tabs: readonly string[];
  active: string | undefined;
  /** Width (left/right) or height (top/bottom) in pixels. */
  size: number;
  /** Collapsed docks show only their tab strip. */
  collapsed: boolean;
  /**
   * Floats over the canvas instead of taking space beside it (desktop mode only; compact mode keeps the bottom sheet).
   * The right and top docks float by default: the canvas runs under them and they read as part of the visualization.
   */
  floating: boolean;
  /**
   * A floating side dock pinned open: the canvas gives it room instead of running under it (some workflows keep a
   * panel open the whole time and want the whole molecule beside it). Only while floating; off by default.
   */
  pinned: boolean;
};

export interface ShellLayout {
  docks: Record<DockId, DockState>;
  mode: 'auto' | ShellMode;
}

export interface ShellOptions {
  mode?: ShellLayout['mode'];
  docks?: Partial<Record<DockId, Partial<DockState>>>;
  /** Below this width (px) an `auto` shell is compact. */
  compactWidth?: number;
  /**
   * Which toolbar buttons show, by item id; ids not listed show. From the UI spec (`layout.buttons`), so an app turns
   * off what it does not want (`{ xr: false, animation: false }`) without unregistering anything.
   */
  buttons?: Readonly<Record<string, boolean>>;
}

/**
 * Px a band tab's backdrop (the canvas' own color behind the sequence) fades out over below its content, so the
 * structure passes under it faded, with no hard cut. The band's view keeps this much clear at its bottom.
 */
export const BandFadePx = 28;
/** Px of the band's last lines over which its backdrop already starts to fade (the dissolve begins within the band). */
export const BandDissolvePx = 36;

const DefaultSizes: Record<DockId, number> = { left: 300, right: 340, top: 88, bottom: 220 };
// Side docks are resized by their splitter; a top or bottom dock is sized by its view (the sequence follows its content).
const MinSizes: Record<DockId, number> = { left: 180, right: 180, top: 40, bottom: 40 };

/**
 * A section contributed to a slot of a panel or overlay: the Export panel's `export` slot (image, session, animation,
 * geometry, models), an inspector's node-kind slot. Extensions add sections instead of buttons, so the chrome does not
 * grow with every extension; commands open the panel at a section.
 */
export interface ShellSection {
  id: string;
  /** The slot it goes to (`export`). */
  slot: string;
  title: string;
  icon?: IconName;
  /** Words it is also found by (the command menu, the section filter). */
  keywords?: readonly string[];
  order?: number;
  /** Shown only when this holds (assembly symmetry only for a symmetric assembly). */
  when?(ctx: ShellViewContext): boolean;
  render: ViewRenderer<ShellViewContext>;
}

/** A button on the canvas toolbar. */
export interface CanvasToolbarItem {
  id: string;
  icon: IconName;
  title: string;
  /**
   * The rail section it belongs to: a built-in one (`canvas`, `export`) or a new one; items without one come first (the
   * command menu). Extensions add their buttons to a section this way (MP4 export, geometry export → `export`). Sections
   * appear in the order of their first item; `order` places items within a section.
   */
  group?: string;
  order?: number;
  /** For rail items: the same as `ShellTab.priority`; which go under More first when the rail is too short. */
  priority?: number;
  /**
   * Where it shows: the rail (default: panels, the command menu, app actions like save and export) or the canvas stack at
   * the top right of the canvas (what acts on the view: camera, lighting, XR, expand), under the mode buttons.
   */
  place?: 'rail' | 'canvas';
  /** A small caption under the icon: the command menu's shortcut, so people learn how to bring it up. */
  caption?: string;
  /** Shown only while this holds (the Sequence button while the sequence card is hidden). */
  visible?: Observable<boolean>;
  /** Shown pressed (a toggle, or the open overlay's button). */
  pressed?: Observable<boolean>;
  /** A face that changes (the theme switch: the current theme's icon and name), in place of `icon` and `title`. */
  face?: Observable<{ icon: IconName; title: string }>;
  /** Opens (or closes) this overlay, anchored to the button. */
  overlay?: string;
  run?(anchor: Element): void;
  /**
   * Variants of the button's action, shown beside it on hover or focus (the camera: reset zoom, orient axes, reset
   * axes); a click on the button still runs its own action, or, without one, opens the menu (illumination and XR under
   * one button). Canvas stack items only.
   */
  menu?: readonly CanvasMenuItem[];
}

/** One entry of a toolbar button's menu. */
export interface CanvasMenuItem {
  id: string;
  label: string;
  icon?: IconName;
  /** A short explanation, shown under the label. */
  hint?: string;
  /** A toggle's (or a single choice's) state: shown with a check mark. */
  pressed?: Observable<boolean>;
  /**
   * One of a single choice (a group of items of which one is current): `pressed` marks the current one as the picking
   * levels' list marks its own, in the accent (icon and label), with no trailing check.
   */
  choice?: boolean;
  /** Items of one group sit together; a hairline separates a group from the one before it. */
  group?: string;
  /** Why it cannot be used here (shown greyed, the reason in place of the hint): fullscreen in an iframe that forbids it. */
  disabled?: string;
  run(): void;
}

/** The open page and what to show in it (`focus`: a tab of the page, an item). */
export interface OpenPage {
  page: ShellPage;
  focus?: string;
}

/** The open overlay, the element it is anchored to (for placement), and the section to bring into view. */
export interface OpenOverlay {
  overlay: ShellOverlay;
  anchor?: Element;
  focus?: string;
}

/** One dock's state and its resolved tabs. */
export class DockModel {
  readonly #state: ObservableRecordState<DockState>;
  readonly #tabs: ObservableState<readonly ShellTab[]>;
  readonly #active: ObservableState<ShellTab | undefined>;

  constructor(
    readonly id: DockId,
    initial: DockState,
  ) {
    this.#state = observableState<DockState>(initial);
    this.#tabs = observableState<readonly ShellTab[]>([]);
    this.#active = observableState<ShellTab | undefined>(undefined);
  }

  get state(): ReadonlyState<DockState> {
    return this.#state.asReadonly;
  }

  /** The tabs in this dock, in order. */
  get tabs(): ReadonlyState<readonly ShellTab[]> {
    return this.#tabs.asReadonly;
  }

  /** The active tab, unless the dock is collapsed. */
  get active(): ReadonlyState<ShellTab | undefined> {
    return this.#active.asReadonly;
  }

  /** @internal */
  update(partial: Partial<DockState>, registered: ReadonlyMap<string, ShellTab>): void {
    this.#state.update(partial);
    const state = this.#state.value;
    const tabs = state.tabs.map((id) => registered.get(id)).filter((t): t is ShellTab => !!t);
    this.#tabs.set(tabs);
    this.#active.set(state.collapsed ? undefined : tabs.find((t) => t.id === state.active));
  }
}

/**
 * The frame around the canvas: docks of tabs, the overlay over the canvas, and the mode (desktop, compact, embedded).
 * Owned by the UI context; persisted layout and contributions come later.
 */
export class ShellModel extends ReactiveModel {
  readonly #registered = new Map<string, ShellTab>();
  readonly #overlays = new Map<string, ShellOverlay>();
  readonly #pages = new Map<string, ShellPage>();
  readonly #page = observableState<OpenPage | undefined>(undefined);
  readonly #docks: Record<DockId, DockModel>;
  readonly #mode = observableState<ShellLayout['mode']>('auto');
  readonly #width = observableState(0);
  readonly #effectiveMode = observableState<ShellMode>('desktop');
  readonly #overlay = observableState<OpenOverlay | undefined>(undefined);
  /** The last request to show a tab at a section (each request a new object, so asking again scrolls again). */
  readonly #tabFocus = observableState<{ tab: string; section: string } | undefined>(undefined);
  /** Every registered item, and the shown ones (`buttons` applied). */
  readonly #toolbarAll = observableState<readonly CanvasToolbarItem[]>([]);
  readonly #toolbar = observableState<readonly CanvasToolbarItem[]>([]);
  #buttons: Readonly<Record<string, boolean>>;
  readonly #sections = observableState<readonly ShellSection[]>([]);
  readonly #slots = new Map<string, ReadonlyState<readonly ShellSection[]>>();
  readonly #compactWidth: number;

  constructor(options: ShellOptions = {}) {
    super();
    this.#compactWidth = options.compactWidth ?? 720;
    this.#buttons = options.buttons ?? {};
    this.#docks = {} as Record<DockId, DockModel>;
    for (const id of DockIds) {
      this.#docks[id] = new DockModel(id, {
        tabs: [],
        active: undefined,
        size: DefaultSizes[id],
        collapsed: false,
        floating: id === 'right' || id === 'top',
        pinned: false,
        ...options.docks?.[id],
      });
    }
    this.#mode.set(options.mode ?? 'auto');
    this.#resolveMode();
  }

  dock(id: DockId): DockModel {
    return this.#docks[id];
  }

  /** The mode in effect: `auto` resolved by the shell's width. */
  get mode(): ReadonlyState<ShellMode> {
    return this.#effectiveMode.asReadonly;
  }

  /**
   * Which sections the user opened or closed (by slot and section id), remembered while the shell lives: a panel
   * shows them as left, across tab switches and re-renders.
   */
  readonly #sectionOpen = new Map<string, boolean>();

  /** Whether the user left a section open or closed; undefined until they touched it. */
  sectionOpen(slot: string, section: string): boolean | undefined {
    return this.#sectionOpen.get(`${slot}\n${section}`);
  }

  setSectionOpen(slot: string, section: string, open: boolean): void {
    this.#sectionOpen.set(`${slot}\n${section}`, open);
  }

  /** A tab shown at a section (`show(tab, section)`): the tab's content gets it as its `focus`. */
  get tabFocus(): ReadonlyState<{ tab: string; section: string } | undefined> {
    return this.#tabFocus.asReadonly;
  }

  get overlay(): ReadonlyState<OpenOverlay | undefined> {
    return this.#overlay.asReadonly;
  }

  /** The page open over the app, if any. */
  get page(): ReadonlyState<OpenPage | undefined> {
    return this.#page.asReadonly;
  }

  /** The canvas toolbar's buttons, in order. */
  get toolbar(): ReadonlyState<readonly CanvasToolbarItem[]> {
    return this.#toolbar.asReadonly;
  }

  /** The sections in a slot, in order. */
  sections(slot: string): ReadonlyState<readonly ShellSection[]> {
    let state = this.#slots.get(slot);
    if (!state) {
      const own = observableState<readonly ShellSection[]>([]);
      const update = (all: readonly ShellSection[]) =>
        own.set(all.filter((s) => s.slot === slot).sort((a, b) => (a.order ?? 0) - (b.order ?? 0)));
      update(this.#sections.value);
      this.subscribePermanently(this.#sections.asReadonly, update);
      state = own.asReadonly;
      this.#slots.set(slot, state);
    }
    return state;
  }

  registerSections(sections: readonly ShellSection[]): () => void {
    this.#sections.set([...this.#sections.value, ...sections]);
    return () => this.#sections.set(this.#sections.value.filter((s) => !sections.includes(s)));
  }

  registerToolbar(items: readonly CanvasToolbarItem[]): () => void {
    const sorted = (list: readonly CanvasToolbarItem[]) => [...list].sort((a, b) => (a.order ?? 0) - (b.order ?? 0));
    this.#toolbarAll.set(sorted([...this.#toolbarAll.value, ...items]));
    this.#updateToolbar();
    return () => {
      this.#toolbarAll.set(this.#toolbarAll.value.filter((i) => !items.includes(i)));
      this.#updateToolbar();
    };
  }

  /** Shows or hides toolbar buttons by id (the UI spec's `layout.buttons`, or a settings toggle). */
  setButtons(buttons: Readonly<Record<string, boolean>>): void {
    this.#buttons = { ...this.#buttons, ...buttons };
    this.#updateToolbar();
  }

  #updateToolbar() {
    this.#toolbar.set(this.#toolbarAll.value.filter((i) => this.#buttons[i.id] !== false));
  }

  /** The layout as it can be persisted. */
  get layout(): ShellLayout {
    const docks = {} as Record<DockId, DockState>;
    for (const id of DockIds) docks[id] = this.#docks[id].state.value;
    return { docks, mode: this.#mode.value };
  }

  /** Adds tabs (to their dock, after the tabs there); returns a function that removes them. */
  registerTabs(tabs: readonly ShellTab[]): () => void {
    const added = [...tabs].sort((a, b) => (a.order ?? 0) - (b.order ?? 0));
    for (const tab of added) {
      this.#registered.set(tab.id, tab);
      const dock = this.#dockOf(tab.id) ?? this.#docks[tab.dock ?? 'right'];
      const state = dock.state.value;
      // A new tab goes by its order among those already there (an extension's tab lands where its order says, not
      // at the end); one already in the dock (a restored layout) keeps its place.
      let ids = state.tabs;
      if (!ids.includes(tab.id)) {
        const order = tab.order ?? 0;
        const at = ids.findIndex((id) => (this.#registered.get(id)?.order ?? 0) > order);
        ids = at < 0 ? [...ids, tab.id] : [...ids.slice(0, at), tab.id, ...ids.slice(at)];
      }
      dock.update({ tabs: ids, active: state.active ?? tab.id }, this.#registered);
    }
    return () => {
      for (const tab of added) {
        if (this.#registered.get(tab.id) !== tab) continue;
        this.#registered.delete(tab.id);
        const dock = this.#dockOf(tab.id);
        if (!dock) continue;
        const state = dock.state.value;
        const ids = state.tabs.filter((id) => id !== tab.id);
        dock.update({ tabs: ids, active: state.active === tab.id ? ids[0] : state.active }, this.#registered);
      }
    };
  }

  registerOverlays(overlays: readonly ShellOverlay[]): () => void {
    for (const o of overlays) this.#overlays.set(o.id, o);
    return () => {
      for (const o of overlays) {
        if (this.#overlays.get(o.id) !== o) continue;
        this.#overlays.delete(o.id);
        if (this.#overlay.value?.overlay === o) this.closeOverlay();
      }
    };
  }

  registerPages(pages: readonly ShellPage[]): () => void {
    for (const p of pages) this.#pages.set(p.id, p);
    return () => {
      for (const p of pages) {
        if (this.#pages.get(p.id) !== p) continue;
        this.#pages.delete(p.id);
        if (this.#page.value?.page === p) this.closePage();
      }
    };
  }

  /** Opens a page over the app (an open card closes); `focus` says what to show in it. */
  openPage(id: string, focus?: string): void {
    const page = this.#pages.get(id);
    if (!page) return;
    this.closeOverlay();
    this.#page.set({ page, ...(focus ? { focus } : {}) });
  }

  closePage(): void {
    this.#page.set(undefined);
  }

  /** Shows a tab: activates it and expands its dock; with `section`, at that section (the Selection tab at Measure). */
  show(tabId: string, section?: string): void {
    const dock = this.#dockOf(tabId);
    if (!dock) return;
    // At a section, or plainly: a plain show forgets an earlier section request for this tab (it was done).
    if (section) this.#tabFocus.set({ tab: tabId, section });
    else if (this.#tabFocus.value?.tab === tabId) this.#tabFocus.set(undefined);
    dock.update({ active: tabId, collapsed: false }, this.#registered);
  }

  /** A tab strip click: activates the tab, or collapses the dock when it is already the active one. */
  select(tabId: string): void {
    const dock = this.#dockOf(tabId);
    if (!dock) return;
    if (this.#tabFocus.value?.tab === tabId) this.#tabFocus.set(undefined);
    const state = dock.state.value;
    if (state.active === tabId && !state.collapsed) dock.update({ collapsed: true }, this.#registered);
    else dock.update({ active: tabId, collapsed: false }, this.#registered);
  }

  toggle(dock: DockId): void {
    const d = this.#docks[dock];
    d.update({ collapsed: !d.state.value.collapsed }, this.#registered);
  }

  /** Moves a tab to another dock (at the end, or at `index`), making it active there. */
  move(tabId: string, to: DockId, index?: number): void {
    const from = this.#dockOf(tabId);
    if (!from || !this.#registered.has(tabId)) return;
    const fromState = from.state.value;
    const rest = fromState.tabs.filter((id) => id !== tabId);
    from.update({ tabs: rest, active: fromState.active === tabId ? rest[0] : fromState.active }, this.#registered);
    const target = this.#docks[to];
    const ids = [...target.state.value.tabs.filter((id) => id !== tabId)];
    ids.splice(index ?? ids.length, 0, tabId);
    target.update({ tabs: ids, active: tabId, collapsed: false }, this.#registered);
  }

  setSize(dock: DockId, size: number): void {
    this.#docks[dock].update({ size: Math.max(MinSizes[dock], Math.round(size)) }, this.#registered);
  }

  /** Floats a dock over the canvas, or docks it beside the canvas. */
  setFloating(dock: DockId, floating: boolean): void {
    this.#docks[dock].update({ floating }, this.#registered);
  }

  /** Pins a floating dock beside the canvas (the canvas shrinks to make room), or lets it float over it again. */
  setPinned(dock: DockId, pinned: boolean): void {
    this.#docks[dock].update({ pinned }, this.#registered);
  }

  setMode(mode: ShellLayout['mode']): void {
    this.#mode.set(mode);
    this.#resolveMode();
  }

  /** The shell's width, reported by its view; resolves `auto`. */
  setWidth(width: number): void {
    this.#width.set(width);
    this.#resolveMode();
  }

  /** Opens an overlay; `focus` names a section to bring into view (a command opening the Export panel at MP4). */
  openOverlay(id: string, anchor?: Element, focus?: string): void {
    const overlay = this.#overlays.get(id);
    if (overlay) this.#overlay.set({ overlay, ...(anchor ? { anchor } : {}), ...(focus ? { focus } : {}) });
  }

  /** Opens the overlay, or closes it when it is the open one. */
  toggleOverlay(id: string, anchor?: Element): void {
    if (this.#overlay.value?.overlay.id === id) this.closeOverlay();
    else this.openOverlay(id, anchor);
  }

  closeOverlay(): void {
    this.#overlay.set(undefined);
  }

  #dockOf(tabId: string): DockModel | undefined {
    for (const id of DockIds) if (this.#docks[id].state.value.tabs.includes(tabId)) return this.#docks[id];
    return undefined;
  }

  #resolveMode(): void {
    const mode = this.#mode.value;
    const width = this.#width.value;
    this.#effectiveMode.set(mode !== 'auto' ? mode : width > 0 && width < this.#compactWidth ? 'compact' : 'desktop');
  }
}
