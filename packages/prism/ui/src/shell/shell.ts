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
import { NavigationManager } from '../managers/navigation.js';
import {
  type DockId,
  DockIds,
  type DockLayout,
  type ShellMode,
  type ViewLayout,
  ViewStateManager,
} from '../managers/view-state.js';
import { batch, type ObservableState, observableState } from '../reactive/observable-state.js';
import { ReactiveModel } from '../reactive/reactive-model.js';
import type { ReadonlyState } from '../reactive/readonly-state.js';

export { type DockId, DockIds, type ShellMode } from '../managers/view-state.js';

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
  /** Stays open on a click outside it, as a pinned card (Esc or its ✕ close it): what is read beside the scene. */
  persistent?: boolean;
  /** Its view takes the card's whole height and scrolls its own parts (the card's body does not scroll). */
  fill?: boolean;
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
   * panel open the whole time and want the whole molecule beside it). Only while floating, on the desktop; the right
   * dock's default (compact and embedded never pin: the layout keeps it for the desktop).
   */
  pinned: boolean;
};

export interface ShellOptions {
  /**
   * The layout's state (api.md §11): docks, sizes, mode, rail, buttons. Given, the shell reads and writes it and
   * the options below that it holds are ignored; absent, the shell makes one from them.
   */
  viewState?: ViewStateManager;
  /** What is open: each dock's tab, the card, the page. Given or made, as `viewState`. */
  navigation?: NavigationManager;
  mode?: ViewLayout['mode'];
  docks?: Partial<Record<DockId, Partial<DockState>>>;
  /** Below this width (px) an `auto` shell is compact. */
  compactWidth?: number;
  /**
   * Which toolbar buttons show, by item id; ids not listed show. From the UI spec (`layout.buttons`), so an app turns
   * off what it does not want (`{ xr: false, animation: false }`) without unregistering anything.
   */
  buttons?: Readonly<Record<string, boolean>>;
  /**
   * A compact rail (not the compact `mode`, the narrow layout): the whole rail column is hidden behind a button in the
   * canvas's corner, which shows and hides it, and the canvas takes its width. Off by default; starts hidden
   * (`railOpen`).
   */
  compactRail?: boolean;
  /**
   * The canvas toolbar (`ShellModel.toolbarCollapsed`, shell.md §2): `collapsed` is one button whose menu lists the
   * tools, `expanded` the full bar; `auto` (the default) is the button in the compact mode (below `compactWidth`, an
   * embedded layout too) and the bar otherwise, following the mode as the shell resizes. The host's choice: the person
   * has no control for it and it is not part of the kept layout.
   */
  toolbar?: 'auto' | 'expanded' | 'collapsed';
}

/**
 * Px a band tab's backdrop (the canvas' own color behind the sequence) fades out over below its content, so the
 * structure passes under it faded, with no hard cut. The band's view keeps this much clear at its bottom.
 */
export const BandFadePx = 28;

/**
 * Ms the canvas stack (its toolbar) takes to fade out when the canvas empties (`--msp-motion-fade`): what is in it is
 * kept that long, so the fade has something to fade (`CanvasStack`'s `quiet`).
 */
export const CanvasStackFadeMs = 360;
/** Px of the band's last lines over which its backdrop already starts to fade (the dissolve begins within the band). */
export const BandDissolvePx = 36;

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
  /** Its header's words when they follow what it offers (Export's "Image & animation"); `title` otherwise. */
  heading?: Observable<string>;
  icon?: IconName;
  /** Words it is also found by (the command menu, the section filter). */
  keywords?: readonly string[];
  order?: number;
  /** Shown only when this holds (assembly symmetry only for a symmetric assembly). */
  when?(ctx: ShellViewContext): boolean;
  /**
   * Told when `when` may answer otherwise for reasons a list cannot see (a PAE fetched for a structure: no change of the
   * scene); a list that shows its sections by `when` asks it again then. Not needed for what the list already follows.
   */
  whenChanged?: ReadonlyState<unknown>;
  /** Open until the person folds it, in a list whose sections start folded (`SectionList`'s `collapsed`). */
  open?: boolean;
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
  /** Its click leaves an open card open (the command menu's button: the menu acts on what is open). */
  keepsCard?: boolean;
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
  /**
   * A submenu: its rows (a list that may change while the menu is open), shown beside the menu on hover, a click or →
   * (the canvas stack's menus), or under the row where a menu has no room beside it (the rail's). Its `run` is not
   * called.
   */
  items?: Observable<readonly CanvasMenuItem[]>;
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

/** One dock's state and its resolved tabs: its layout (`viewState`) and what is open in it (`navigation`). */
export class DockModel {
  readonly #state: ObservableState<DockState>;
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
  sync(state: DockState, registered: ReadonlyMap<string, ShellTab>, hidden: readonly string[]): void {
    this.#state.set(state);
    // Hidden tabs (`ShellModel.hideTab`) keep their place in the layout and are left out here.
    const tabs = state.tabs
      .filter((id) => !hidden.includes(id))
      .map((id) => registered.get(id))
      .filter((t): t is ShellTab => !!t);
    this.#tabs.set(tabs);
    this.#active.set(state.collapsed ? undefined : tabs.find((t) => t.id === state.active));
  }
}

const LayoutKeys = ['tabs', 'size', 'floating', 'pinned'] as const;

/**
 * The frame around the canvas: the parts registered into it (tabs, overlays, pages, sections, toolbar items) and what
 * a person does with them (show, select, move, resize, open a card). Its state is two managers' (api.md §11):
 * `viewState` (the layout) and `navigation` (what is open), so an app keeps, restores and drives them without a view.
 */
export class ShellModel extends ReactiveModel {
  readonly #registered = new Map<string, ShellTab>();
  readonly #overlays = new Map<string, ShellOverlay>();
  readonly #pages = new Map<string, ShellPage>();
  readonly #viewState: ViewStateManager;
  readonly #navigation: NavigationManager;
  readonly #page = observableState<OpenPage | undefined>(undefined);
  readonly #docks: Record<DockId, DockModel>;
  readonly #railShown = observableState(true);
  readonly #width = observableState(0);
  readonly #effectiveMode = observableState<ShellMode>('desktop');
  // Compared by identity: each request is a new object, so asking again for the same section is not swallowed.
  readonly #overlay = observableState<OpenOverlay | undefined>(undefined, { compare: Object.is });
  /** The element the open card is anchored to: a view's, not state. */
  #anchor: { id: string; element: Element } | undefined;
  readonly #tabFocus = observableState<{ tab: string; section: string } | undefined>(undefined, {
    compare: Object.is,
  });
  readonly #compactRail = observableState(false);
  readonly #railOpen = observableState(false);
  readonly #hiddenTabs = observableState<readonly ShellTab[]>([]);
  /** Every registered item, and the shown ones (`buttons` applied). */
  readonly #toolbarAll = observableState<readonly CanvasToolbarItem[]>([]);
  readonly #toolbar = observableState<readonly CanvasToolbarItem[]>([]);
  readonly #sections = observableState<readonly ShellSection[]>([]);
  readonly #slots = new Map<string, ReadonlyState<readonly ShellSection[]>>();
  readonly #compactWidth: number;
  readonly #toolbarOption: 'auto' | 'expanded' | 'collapsed';
  readonly #toolbarCollapsed = observableState(false);

  constructor(options: ShellOptions = {}) {
    super();
    this.#compactWidth = options.compactWidth ?? 720;
    this.#toolbarOption = options.toolbar ?? 'auto';
    this.#viewState = options.viewState ?? new ViewStateManager(layoutOptions(options));
    this.#navigation = options.navigation ?? new NavigationManager(navigationOptions(options));
    this.#docks = {} as Record<DockId, DockModel>;
    for (const id of DockIds) this.#docks[id] = new DockModel(id, this.#dockState(id));
    this.subscribePermanently(this.#viewState.layout, () => this.#sync());
    this.subscribePermanently(this.#navigation.state, () => this.#sync());
  }

  /** The layout's state (docks, sizes, mode, rail, buttons). */
  get viewState(): ViewStateManager {
    return this.#viewState;
  }

  /** What is open (each dock's tab, the card, the page). */
  get navigation(): NavigationManager {
    return this.#navigation;
  }

  /** The rail is a column hidden behind the canvas corner's button (`ShellOptions.compactRail`). */
  get compactRail(): ReadonlyState<boolean> {
    return this.#compactRail.asReadonly;
  }

  /** In a compact rail: whether the column is shown (it is hidden until the button shows it). */
  get railOpen(): ReadonlyState<boolean> {
    return this.#railOpen.asReadonly;
  }

  /** Whether the rail column is shown: always, unless the rail is compact and not open. */
  get railShown(): ReadonlyState<boolean> {
    return this.#railShown.asReadonly;
  }

  /** The canvas toolbar is one button with a menu (`ShellOptions.toolbar`, shell.md §2); follows the mode under `auto`. */
  get toolbarCollapsed(): ReadonlyState<boolean> {
    return this.#toolbarCollapsed.asReadonly;
  }

  setCompactRail(compact: boolean): void {
    this.#viewState.setCompactRail(compact);
    this.setRailOpen(false);
  }

  /**
   * Shows or hides a compact rail's column. Hiding it closes what opened from it: the overlay card and the right panel,
   * floating or pinned (a pinned panel without its rail would have no way to switch tabs; the pin is kept for when it
   * is opened again).
   */
  setRailOpen(open: boolean): void {
    // What opened from it closes first, while the rail is still shown: hidden UI is paused (styling.md §5), so a panel
    // closed after the rail hid would only see it when the rail came back, and flash open and shut then.
    if (this.#compactRail.value && !open) {
      this.closeOverlay();
      this.#navigation.setDock('right', { collapsed: true });
    }
    this.#navigation.setRailOpen(open);
  }

  toggleRail(): void {
    if (this.#compactRail.value) this.setRailOpen(!this.#navigation.state.value.railOpen);
  }

  /**
   * The canvas fills the viewer (`navigation.expanded`, `view.expand`): the view leaves out the rail, the docks and the
   * sequence band; the canvas stack stays, with the way back (and Escape).
   */
  get expanded(): ReadonlyState<boolean> {
    return this.#navigation.expanded;
  }

  setExpanded(expanded: boolean): void {
    if (expanded) this.closeOverlay();
    this.#navigation.setExpanded(expanded);
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
    this.#viewState.setButtons(buttons);
  }

  /** Adds tabs (to their dock, after the tabs there); returns a function that removes them. */
  registerTabs(tabs: readonly ShellTab[]): () => void {
    const added = [...tabs].sort((a, b) => (a.order ?? 0) - (b.order ?? 0));
    for (const tab of added) {
      this.#registered.set(tab.id, tab);
      this.#place(tab);
    }
    this.#sync();
    return () => {
      for (const tab of added) {
        if (this.#registered.get(tab.id) !== tab) continue;
        this.#registered.delete(tab.id);
        const dock = this.#dockOf(tab.id);
        if (!dock) continue;
        const state = this.#dockState(dock);
        const ids = state.tabs.filter((id) => id !== tab.id);
        this.#update(dock, { tabs: ids, active: state.active === tab.id ? ids[0] : state.active });
      }
      this.#sync();
    };
  }

  /**
   * Puts a registered tab in its dock: where the layout has it, else its own dock (`dock`, default right) by its order
   * among those already there (an extension's tab lands where its order says, not at the end).
   */
  #place(tab: ShellTab): void {
    const dock = this.#dockOf(tab.id) ?? tab.dock ?? 'right';
    const state = this.#dockState(dock);
    let ids = state.tabs;
    if (!ids.includes(tab.id)) {
      const order = tab.order ?? 0;
      const at = ids.findIndex((id) => (this.#registered.get(id)?.order ?? 0) > order);
      ids = at < 0 ? [...ids, tab.id] : [...ids.slice(0, at), tab.id, ...ids.slice(at)];
    }
    // A hidden tab does not become the dock's open one.
    const hidden = this.#viewState.layout.value.hiddenTabs.includes(tab.id);
    this.#update(dock, { tabs: ids, active: state.active ?? (hidden ? undefined : tab.id) });
  }

  /** The registered tabs taken out of view (`hideTab`), in the order they were hidden. */
  get hiddenTabs(): ReadonlyState<readonly ShellTab[]> {
    return this.#hiddenTabs.asReadonly;
  }

  /**
   * Takes a tab out of view without unregistering it (shell.md §3): its button and panel go, it keeps its dock and
   * place, and the layout keeps it hidden (`viewState`, so a kept layout keeps it so). Its open panel closes.
   */
  hideTab(tabId: string): void {
    const dock = this.#dockOf(tabId);
    batch(() => {
      this.#viewState.setTabHidden(tabId, true);
      if (!dock) return;
      const state = this.#dockState(dock);
      if (state.active !== tabId) return;
      const hidden = this.#viewState.layout.value.hiddenTabs;
      const next = state.tabs.find((id) => id !== tabId && !hidden.includes(id) && this.#registered.has(id));
      this.#update(dock, { active: next, collapsed: true });
    });
  }

  /** Brings a hidden tab back where it was (closed: `show` opens it). */
  showTab(tabId: string): void {
    this.#viewState.setTabHidden(tabId, false);
  }

  /** Brings every hidden tab back. */
  showAllTabs(): void {
    this.#viewState.showAllTabs();
  }

  /**
   * The layout back to its defaults (shell.md §3): `viewState.reset()` (the app's seed over `defaultViewLayout()`),
   * then every registered tab the seed does not place goes back to its own dock by its order; hidden tabs show again.
   * Each dock keeps its open tab where that tab is still in it, and whether it is collapsed.
   */
  resetLayout(): void {
    batch(() => {
      this.#viewState.reset();
      const tabs = [...this.#registered.values()].sort((a, b) => (a.order ?? 0) - (b.order ?? 0));
      for (const tab of tabs) this.#place(tab);
      for (const id of DockIds) {
        const state = this.#dockState(id);
        if (state.active !== undefined && !state.tabs.includes(state.active)) {
          this.#navigation.setDock(id, { active: state.tabs.find((t) => this.#registered.has(t)) });
        }
      }
    });
  }

  registerOverlays(overlays: readonly ShellOverlay[]): () => void {
    for (const o of overlays) this.#overlays.set(o.id, o);
    this.#sync();
    return () => {
      for (const o of overlays) {
        if (this.#overlays.get(o.id) !== o) continue;
        this.#overlays.delete(o.id);
        if (this.#navigation.state.value.card?.id === o.id) this.closeOverlay();
      }
      this.#sync();
    };
  }

  registerPages(pages: readonly ShellPage[]): () => void {
    for (const p of pages) this.#pages.set(p.id, p);
    this.#sync();
    return () => {
      for (const p of pages) {
        if (this.#pages.get(p.id) !== p) continue;
        this.#pages.delete(p.id);
        if (this.#navigation.state.value.page?.id === p.id) this.closePage();
      }
      this.#sync();
    };
  }

  /** Opens a page over the app (an open card closes); `focus` says what to show in it. */
  openPage(id: string, focus?: string): void {
    if (!this.#pages.has(id)) return;
    this.closeOverlay();
    this.#navigation.setPage({ id, ...(focus ? { focus } : {}) });
  }

  closePage(): void {
    this.#navigation.setPage(undefined);
  }

  /** Shows a tab: activates it and expands its dock; with `section`, at that section (the Selection tab at Measure). */
  show(tabId: string, section?: string): void {
    const dock = this.#dockOf(tabId);
    if (!dock) return;
    // Asked for by name (a command, a link): a hidden tab comes back.
    this.showTab(tabId);
    // At a section, or plainly: a plain show forgets an earlier section request for this tab (it was done).
    if (section) this.#navigation.setFocus({ tab: tabId, section });
    else if (this.#navigation.state.value.focus?.tab === tabId) this.#navigation.setFocus(undefined);
    this.#update(dock, { active: tabId, collapsed: false });
  }

  /** A tab strip click: activates the tab, or collapses the dock when it is already the active one. */
  select(tabId: string): void {
    const dock = this.#dockOf(tabId);
    if (!dock) return;
    if (this.#navigation.state.value.focus?.tab === tabId) this.#navigation.setFocus(undefined);
    const state = this.#dockState(dock);
    if (state.active === tabId && !state.collapsed) this.#update(dock, { collapsed: true });
    else this.#update(dock, { active: tabId, collapsed: false });
  }

  toggle(dock: DockId): void {
    this.#update(dock, { collapsed: !this.#dockState(dock).collapsed });
  }

  /** Moves a tab to another dock (at the end, or at `index`), making it active there. */
  move(tabId: string, to: DockId, index?: number): void {
    const from = this.#dockOf(tabId);
    if (!from || !this.#registered.has(tabId)) return;
    const fromState = this.#dockState(from);
    const rest = fromState.tabs.filter((id) => id !== tabId);
    this.#update(from, { tabs: rest, active: fromState.active === tabId ? rest[0] : fromState.active });
    const ids = [...this.#dockState(to).tabs.filter((id) => id !== tabId)];
    ids.splice(index ?? ids.length, 0, tabId);
    this.#update(to, { tabs: ids, active: tabId, collapsed: false });
  }

  setSize(dock: DockId, size: number): void {
    this.#viewState.setDock(dock, { size });
  }

  /** Floats a dock over the canvas, or docks it beside the canvas. */
  setFloating(dock: DockId, floating: boolean): void {
    this.#viewState.setDock(dock, { floating });
  }

  /** Pins a floating dock beside the canvas (the canvas shrinks to make room), or lets it float over it again. */
  setPinned(dock: DockId, pinned: boolean): void {
    this.#viewState.setDock(dock, { pinned });
  }

  /** The layout's mode; `embedded` starts with the rail hidden (shell.md §3). */
  setMode(mode: ViewLayout['mode']): void {
    this.#viewState.setMode(mode);
    if (mode === 'embedded') this.setRailOpen(false);
  }

  /** The shell's width, reported by its view; resolves `auto` (and a narrow `embedded`). */
  setWidth(width: number): void {
    this.#width.set(width);
    if (this.#resolveMode()) this.#sync();
  }

  /** Opens an overlay; `focus` names a section to bring into view (a command opening the Export panel at MP4). */
  openOverlay(id: string, anchor?: Element, focus?: string): void {
    if (!this.#overlays.has(id)) return;
    this.#anchor = anchor ? { id, element: anchor } : this.#anchor?.id === id ? this.#anchor : undefined;
    const before = this.#overlay.value;
    this.#navigation.setCard({ id, ...(focus ? { section: focus } : {}) });
    // Asked again at the same section (the sync found nothing new): the card brings it into view again.
    if (this.#overlay.value === before) this.#syncOverlay(true);
  }

  /** Opens the overlay, or closes it when it is the open one. */
  toggleOverlay(id: string, anchor?: Element): void {
    if (this.#navigation.state.value.card?.id === id) this.closeOverlay();
    else this.openOverlay(id, anchor);
  }

  closeOverlay(): void {
    this.#anchor = undefined;
    this.#navigation.setCard(undefined);
  }

  #dockOf(tabId: string): DockId | undefined {
    const docks = this.#viewState.layout.value.docks;
    for (const id of DockIds) if (docks[id].tabs.includes(tabId)) return id;
    return undefined;
  }

  #dockState(id: DockId): DockState {
    const layout = this.#viewState.layout.value.docks[id];
    const nav = this.#navigation.state.value.docks[id];
    const state = { ...layout, active: nav.active, collapsed: nav.collapsed };
    // Embedded (shell.md §3): the right panel floats and nothing is pinned, whatever the layout says (it is kept for
    // another mode); a viewer inside a page has no width to give a panel.
    if (this.#effectiveMode.value === 'embedded') {
      return { ...state, floating: id === 'right' ? true : state.floating, pinned: false };
    }
    // A compact rail lies over the canvas and takes no width (shell.md §2): nothing pinned beside it either (the pin is
    // kept for a full rail). Pinned, the canvas ended at the panel while the rail's button stayed in the canvas's corner,
    // over the toolbar, and the rail came and went beside a fixed panel. Only while the dock floats: docked, there is
    // no compact rail (the view's `compactRail` is floating's), and the pin keeps what it does.
    if (this.#viewState.layout.value.compactRail && state.floating && state.pinned) return { ...state, pinned: false };
    // Compact (a phone): the panel is a bottom sheet, never pinned (the person's pin is kept for the desktop): Escape
    // and a click outside close it; a card opens in the sheet in its place, closed the same ways (not a modal dialog:
    // Tab leaves it, as it leaves the panel).
    if (this.#effectiveMode.value === 'compact' && state.pinned) return { ...state, pinned: false };
    return state;
  }

  /** Writes a dock's change to the manager that holds each part. */
  #update(id: DockId, partial: Partial<DockState>): void {
    const layout: Partial<{ -readonly [K in keyof DockLayout]: DockLayout[K] }> = {};
    for (const key of LayoutKeys) if (key in partial) (layout as Record<string, unknown>)[key] = partial[key];
    if (Object.keys(layout).length) this.#viewState.setDock(id, layout);
    const nav: { active?: string; collapsed?: boolean } = {};
    if ('active' in partial) nav.active = partial.active;
    if ('collapsed' in partial) nav.collapsed = partial.collapsed;
    if (Object.keys(nav).length) this.#navigation.setDock(id, nav);
  }

  /** The views' states, from the managers and the registries. */
  #sync(): void {
    if (!this.#docks) return;
    this.#resolveMode();
    const layout = this.#viewState.layout.value;
    for (const id of DockIds) this.#docks[id].sync(this.#dockState(id), this.#registered, layout.hiddenTabs);
    const hidden = layout.hiddenTabs.map((id) => this.#registered.get(id)).filter((t): t is ShellTab => !!t);
    const before = this.#hiddenTabs.value;
    if (hidden.length !== before.length || hidden.some((t, i) => t !== before[i])) this.#hiddenTabs.set(hidden);
    const nav = this.#navigation.state.value;
    // Embedded, the rail is always the compact one (shell.md §3).
    const compactRail = layout.compactRail || this.#effectiveMode.value === 'embedded';
    this.#compactRail.set(compactRail);
    this.#railOpen.set(nav.railOpen);
    this.#railShown.set(!compactRail || nav.railOpen);
    this.#tabFocus.set(nav.focus);
    const page = nav.page ? this.#pages.get(nav.page.id) : undefined;
    const open = this.#page.value;
    if (!page) this.#page.set(undefined);
    else if (open?.page !== page || open.focus !== nav.page?.focus) {
      this.#page.set({ page, ...(nav.page?.focus ? { focus: nav.page.focus } : {}) });
    }
    this.#syncOverlay(false);
    this.#updateToolbar();
  }

  #syncOverlay(again: boolean): void {
    const card = this.#navigation.state.value.card;
    const overlay = card ? this.#overlays.get(card.id) : undefined;
    if (!overlay || !card) {
      this.#overlay.set(undefined);
      return;
    }
    const current = this.#overlay.value;
    const anchor = this.#anchor?.id === card.id ? this.#anchor.element : undefined;
    if (!again && current?.overlay === overlay && current.anchor === anchor && current.focus === card.section) return;
    this.#overlay.set({ overlay, ...(anchor ? { anchor } : {}), ...(card.section ? { focus: card.section } : {}) });
  }

  #updateToolbar() {
    const buttons = this.#viewState.layout.value.buttons;
    this.#toolbar.set(this.#toolbarAll.value.filter((i) => buttons[i.id] !== false));
  }

  /** The mode in effect; whether it changed. `auto` and `embedded` are compact below `compactWidth` (shell.md §3). */
  #resolveMode(): boolean {
    const mode = this.#viewState.layout.value.mode;
    const width = this.#width.value;
    const narrow = width > 0 && width < this.#compactWidth;
    const next: ShellMode =
      mode === 'auto' ? (narrow ? 'compact' : 'desktop') : mode === 'embedded' && narrow ? 'compact' : mode;
    const option = this.#toolbarOption;
    this.#toolbarCollapsed.set(option === 'auto' ? next === 'compact' : option === 'collapsed');
    if (this.#effectiveMode.value === next) return false;
    this.#effectiveMode.set(next);
    return true;
  }
}

function layoutOptions(options: ShellOptions): ConstructorParameters<typeof ViewStateManager>[0] {
  const docks: Partial<Record<DockId, Partial<DockLayout>>> = {};
  for (const id of DockIds) {
    const d = options.docks?.[id];
    if (!d) continue;
    const layout: Record<string, unknown> = {};
    for (const key of LayoutKeys) if (key in d) layout[key] = d[key];
    docks[id] = layout as Partial<DockLayout>;
  }
  return {
    docks,
    ...(options.mode ? { mode: options.mode } : {}),
    compactRail: !!options.compactRail,
    buttons: options.buttons ?? {},
  };
}

function navigationOptions(options: ShellOptions): ConstructorParameters<typeof NavigationManager>[0] {
  const docks: Partial<Record<DockId, { active?: string; collapsed?: boolean }>> = {};
  for (const id of DockIds) {
    const d = options.docks?.[id];
    if (!d) continue;
    docks[id] = {
      ...('active' in d ? { active: d.active } : {}),
      ...(d.collapsed !== undefined ? { collapsed: d.collapsed } : {}),
    };
  }
  return { docks };
}
