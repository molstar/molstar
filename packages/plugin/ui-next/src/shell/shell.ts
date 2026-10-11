/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { IconName } from '../components/icon/icon.js';
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
  title: string;
  icon: IconName;
  /** Where it goes by default ('right'). */
  dock?: DockId;
  order?: number;
  render: ViewRenderer<ShellViewContext>;
}

/** A card over the canvas, opened by a toolbar button or a command. */
export interface ShellOverlay {
  id: string;
  title: string;
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
}

const DefaultSizes: Record<DockId, number> = { left: 300, right: 340, top: 88, bottom: 220 };
const MinSize = 180;

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
  /** Words it is also found by (the palette, the section filter). */
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
   * palette). Extensions add their buttons to a section this way (MP4 export, geometry export → `export`). Sections
   * appear in the order of their first item; `order` places items within a section.
   */
  group?: string;
  order?: number;
  /** Shown pressed (a toggle, or the open overlay's button). */
  pressed?: ReadonlyState<boolean>;
  /** Opens (or closes) this overlay, anchored to the button. */
  overlay?: string;
  run?(anchor: Element): void;
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
  readonly #docks: Record<DockId, DockModel>;
  readonly #mode = observableState<ShellLayout['mode']>('auto');
  readonly #width = observableState(0);
  readonly #effectiveMode = observableState<ShellMode>('desktop');
  readonly #overlay = observableState<OpenOverlay | undefined>(undefined);
  readonly #toolbar = observableState<readonly CanvasToolbarItem[]>([]);
  readonly #sections = observableState<readonly ShellSection[]>([]);
  readonly #slots = new Map<string, ReadonlyState<readonly ShellSection[]>>();
  readonly #compactWidth: number;

  constructor(options: ShellOptions = {}) {
    super();
    this.#compactWidth = options.compactWidth ?? 720;
    this.#docks = {} as Record<DockId, DockModel>;
    for (const id of DockIds) {
      this.#docks[id] = new DockModel(id, {
        tabs: [],
        active: undefined,
        size: DefaultSizes[id],
        collapsed: false,
        floating: id === 'right' || id === 'top',
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

  get overlay(): ReadonlyState<OpenOverlay | undefined> {
    return this.#overlay.asReadonly;
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
    this.#toolbar.set(sorted([...this.#toolbar.value, ...items]));
    return () => this.#toolbar.set(this.#toolbar.value.filter((i) => !items.includes(i)));
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
      const ids = state.tabs.includes(tab.id) ? state.tabs : [...state.tabs, tab.id];
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

  /** Shows a tab: activates it and expands its dock. */
  show(tabId: string): void {
    const dock = this.#dockOf(tabId);
    if (dock) dock.update({ active: tabId, collapsed: false }, this.#registered);
  }

  /** A tab strip click: activates the tab, or collapses the dock when it is already the active one. */
  select(tabId: string): void {
    const dock = this.#dockOf(tabId);
    if (!dock) return;
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
    this.#docks[dock].update({ size: Math.max(MinSize, Math.round(size)) }, this.#registered);
  }

  /** Floats a dock over the canvas, or docks it beside the canvas. */
  setFloating(dock: DockId, floating: boolean): void {
    this.#docks[dock].update({ floating }, this.#registered);
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
