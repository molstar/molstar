/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { batch, derive, listen, type ReadonlyState, State } from '@molstar/core/reactive';
import type { Child } from '../dom/index.js';
import type { ViewRenderer } from '../dom/view.js';
import type { IconName } from '../components/icon/icon.js';
import {
  allGroups,
  clampRect,
  defaultLayout,
  findGroup,
  findView,
  type LayoutTarget,
  mergeGroup,
  moveInColumn,
  moveView,
  placeWindow,
  type Rect,
  RightGroup,
  railGroup,
  raiseWindow,
  removeView,
  resizeNode,
  selectView,
  setCollapsed,
  setGroupState,
  setHidden,
  setRail,
  setWindowRect,
  type ShellLayout,
  stacked,
  validateLayout,
  type ViewRef,
} from './layout.js';

/** Where a view is shown (shell §3): it may look different in each. */
export type ViewPlacement = 'docked' | 'window' | 'card' | 'page' | 'modal' | 'sheet';

/** What a view is rendered with. */
export interface ShellViewContext {
  readonly shell: ShellModel;
  /** The view's ref in the layout (its id, or `id#key` for an instance). */
  readonly ref: ViewRef;
  /** The instance's key (`instances: 'many'`). */
  readonly key?: string;
  /** What it is about (given to `show`): the structure a plot is of, an item an inspector shows. */
  readonly subject: ReadonlyState<unknown>;
  /** The section to bring into view (`show(id, { section })`), anew at each such `show`. */
  readonly focus: ReadonlyState<string | undefined>;
  readonly placement: ViewPlacement;
  /** Closes it where it is (a card, a window, an instance; a docked view is hidden). */
  close(): void;
}

/** Where a view goes when it is first shown (shell §2): a docked group, a window, or one of the shell's layers. */
export type ViewHome = string | 'float' | 'card' | 'page' | 'modal';

/** A view contributed by the app or an extension (shell §2: what a tab, a card, a dialog or a page is). */
export interface ShellViewDef {
  readonly id: string;
  /** Its name: the rail's tooltip, a tab's label, and its header unless `heading` says otherwise. */
  readonly title: string;
  /** Its header when it differs from its name: shorter text, or a view (a wordmark). */
  readonly heading?: string | (() => Child);
  readonly icon: IconName;
  readonly render: ViewRenderer<ShellViewContext>;
  /** What goes in its header after the title (a page's own tabs, a view's actions). */
  readonly actions?: ViewRenderer<ShellViewContext>;
  /** Its place among the views of its group (across contributors). */
  readonly order?: number;
  /** On a rail: its head (Home, a section of its own), the strip (default), or its foot (Settings, Help). */
  readonly place?: 'head' | 'strip' | 'foot';
  /** While false the view leaves its strip (it keeps its place and comes back there); a view that exists only sometimes. */
  readonly visible?: ReadonlyState<boolean>;
  /** Which rail buttons go under More first when the rail is too short: the higher, the sooner. */
  readonly priority?: number;
  /** In the compact layout's bar (the rest go under More there). */
  readonly compactBar?: boolean;
  /** Its size as a window (default 360 × 320) and the least it takes. */
  readonly size?: { readonly width: number; readonly height: number };
  readonly minSize?: { readonly width: number; readonly height: number };
  /** Where it goes when first shown (default the right group). */
  readonly home?: ViewHome;
  /** `many`: the app opens instances by a key, each its own place and state (two plots side by side). */
  readonly instances?: 'single' | 'many';
  /** False: its geometry is not part of a saved layout (it opens at its default place and size each time). */
  readonly remember?: boolean;
  /** A card that stays on a press outside it (Escape or its ✕ close it). */
  readonly persistent?: boolean;
  /** Its content takes the body's whole height and scrolls its own parts (the body does not scroll). */
  readonly fill?: boolean;
}

/** A section contributed to a slot of a view (an Export view's `export` slot): extensions add sections, not views. */
export interface ShellSection {
  readonly id: string;
  readonly slot: string;
  readonly title: string;
  /** Its header's words when they follow what it offers. */
  readonly heading?: ReadonlyState<string>;
  readonly icon?: IconName;
  /** Words it is also found by. */
  readonly keywords?: readonly string[];
  readonly order?: number;
  /** Shown only while this holds. */
  readonly when?: (ctx: ShellViewContext) => boolean;
  /** Open until folded, in a list whose sections start folded. */
  readonly open?: boolean;
  readonly render: ViewRenderer<ShellViewContext>;
}

/** A button on the rail that is no view: an app-wide tool (the theme, the command menu). */
export interface RailItem {
  readonly id: string;
  readonly icon: IconName;
  readonly title: string;
  /** Items of one group sit together, a hairline between groups; items without one go at the foot. */
  readonly group?: string;
  readonly order?: number;
  readonly priority?: number;
  /** A small caption under its icon (a shortcut). */
  readonly caption?: string;
  readonly pressed?: ReadonlyState<boolean>;
  readonly visible?: ReadonlyState<boolean>;
  run(anchor: Element): void;
}

/** A view open in one of the shell's layers (a card, a page, a modal). */
export interface OpenLayer {
  readonly ref: ViewRef;
  readonly def: ShellViewDef;
  /** Where it was opened from (focus goes back there). */
  readonly anchor?: Element;
}

export type ShellMode = 'desktop' | 'compact' | 'embedded';

export interface ShellOptions {
  /** `auto` (default): compact below `compactWidth` of the shell's own width, desktop above. */
  readonly mode?: 'auto' | ShellMode;
  readonly compactWidth?: number;
  /** The layout to start from (a stored one: `validateLayout` is applied); default: the right group. */
  readonly layout?: unknown;
  /** The right group's width (default 340). */
  readonly size?: number;
}

/** Options of `show`. */
export interface ShowOptions {
  /** The instance (`instances: 'many'`). */
  readonly key?: string;
  readonly subject?: unknown;
  /** The section to bring into view. */
  readonly section?: string;
  /** What opened it: a card or a page gives focus back there; a new window opens beside it (`anchorRect`). */
  readonly anchor?: Element;
  /** The anchor's rectangle relative to the main area (where a new window opens). */
  readonly anchorRect?: Rect;
  /** Where it goes, for this once (a window instead of its home). */
  readonly target?: LayoutTarget;
}

const DefaultWindow = { width: 360, height: 320 } as const;

/** A view's ref: its id, or `id#key`. */
export function viewRef(id: string, key?: string): ViewRef {
  return key === undefined ? id : `${id}#${key}`;
}

/** The view id of a ref. */
export function viewIdOf(ref: ViewRef): string {
  const i = ref.indexOf('#');
  return i < 0 ? ref : ref.slice(0, i);
}

const isLayerHome = (home: ViewHome | undefined) => home === 'card' || home === 'page' || home === 'modal';

/**
 * The shell's state (shell §2), DOM-free: what is registered (views, sections, rail items), the layout (a `State` of
 * a `ShellLayout`, the value an app saves and restores), what is open in the layers over it (a card, a page, a modal),
 * the mode. Operations are named for what a person does (show, select, move, float, dock, close, hide); each is one
 * change of the layout.
 */
export class ShellModel {
  readonly #views = new State<ReadonlyMap<string, ShellViewDef>>(new Map());
  readonly #sections = new State<readonly ShellSection[]>([]);
  readonly #railItems = new State<readonly RailItem[]>([]);
  /** Views whose `visible` is false now: out of their strips, their places kept. */
  readonly #absent = new State<ReadonlySet<string>>(new Set());
  readonly #layout: State<ShellLayout>;
  readonly #defaultLayout: ShellLayout;
  readonly #card = new State<OpenLayer | undefined>(undefined);
  readonly #page = new State<OpenLayer | undefined>(undefined);
  readonly #modal = new State<OpenLayer | undefined>(undefined);
  readonly #width = new State(0);
  readonly #area = new State({ width: 0, height: 0 });
  readonly #modeOption: 'auto' | ShellMode;
  readonly #compactWidth: number;
  /** The right group's width when it is made (`ShellOptions.size`). */
  readonly #size: number;
  readonly #mode: ReadonlyState<ShellMode>;
  readonly #subjects = new Map<ViewRef, State<unknown>>();
  readonly #focus = new Map<ViewRef, State<string | undefined>>();
  readonly #sectionOpen = new Map<string, boolean>();
  /** Each view's last window rectangle this session (a window closed and opened again comes back there). */
  readonly #lastRects = new Map<ViewRef, Rect>();

  constructor(options: ShellOptions = {}) {
    this.#modeOption = options.mode ?? 'auto';
    this.#compactWidth = options.compactWidth ?? 720;
    this.#size = options.size ?? 340;
    this.#defaultLayout = defaultLayout([], { size: options.size });
    this.#layout = new State(validateLayout(options.layout) ?? this.#defaultLayout);
    this.#mode = derive(this.#width, (w) =>
      this.#modeOption !== 'auto' ? this.#modeOption : w > 0 && w < this.#compactWidth ? 'compact' : 'desktop',
    );
    if (this.#modeOption === 'embedded') this.#layout.set(setRail(this.#layout.value, false));
  }

  // ---------------------------------------------------------------------------------------------------------------
  // What is registered

  get views(): ReadonlyState<ReadonlyMap<string, ShellViewDef>> {
    return this.#views.asReadonly;
  }

  view(id: string): ShellViewDef | undefined {
    return this.#views.value.get(viewIdOf(id));
  }

  /** Views whose `visible` flag is false: their strips leave them out. */
  get absent(): ReadonlyState<ReadonlySet<string>> {
    return this.#absent.asReadonly;
  }

  /** Whether a view shows in its strip: registered, not hidden by the person, its flag not false. */
  shows(ref: ViewRef): boolean {
    return !!this.view(ref) && !this.#absent.value.has(viewIdOf(ref)) && !this.#layout.value.hidden.includes(ref);
  }

  #setAbsent(id: string, absent: boolean): void {
    const now = this.#absent.value;
    if (now.has(id) === absent) return;
    const next = new Set(now);
    if (absent) next.add(id);
    else next.delete(id);
    this.#absent.set(next);
  }

  get railItems(): ReadonlyState<readonly RailItem[]> {
    return this.#railItems.asReadonly;
  }

  /**
   * Adds views; returns their undo (an extension unloads). A view with a docked home is placed there now, by its order
   * (unless the layout has it already: a stored layout keeps where the person put it); a view with another home waits
   * for `show`.
   */
  registerViews(defs: readonly ShellViewDef[]): () => void {
    const stops: (() => void)[] = [];
    batch(() => {
      for (const d of defs) {
        if (!d.visible) continue;
        const flag = d.visible;
        const write = () => this.#setAbsent(d.id, !flag.value);
        write();
        stops.push(listen(flag, write));
      }
      const next = new Map(this.#views.value);
      for (const d of defs) next.set(d.id, d);
      this.#views.set(next);
      this.#layout.set(this.#placeMissing(this.#layout.value, defs));
    });
    return () => {
      for (const stop of stops) stop();
      const next = new Map(this.#views.value);
      for (const d of defs) if (next.get(d.id) === d) next.delete(d.id);
      for (const d of defs) if (d.visible) this.#setAbsent(d.id, false);
      this.#views.set(next);
      // Its place stays in the layout (a view that registers again lands there); what showed it in a layer goes.
      for (const layer of [this.#card, this.#page, this.#modal]) {
        if (layer.value && defs.includes(layer.value.def)) layer.set(undefined);
      }
    };
  }

  /** Every single view with a docked home that the layout does not have, placed at its home by its order. */
  #placeMissing(layout: ShellLayout, defs: readonly ShellViewDef[]): ShellLayout {
    let l = layout;
    for (const d of defs) {
      if ((d.instances ?? 'single') !== 'single' || isLayerHome(d.home) || d.home === 'float') continue;
      if (findView(l, d.id)) continue;
      l = this.#placeAtHome(l, d.id, d, false);
    }
    return l;
  }

  /** `ref` into its def's docked home group (made beside the main area when it is not there), at its order. */
  #placeAtHome(layout: ShellLayout, ref: ViewRef, def: ShellViewDef, activate: boolean): ShellLayout {
    const home = typeof def.home === 'string' && !isLayerHome(def.home) && def.home !== 'float' ? def.home : RightGroup;
    const place = findGroup(layout, home);
    let l: ShellLayout;
    if (place) {
      const views = place.group.views;
      const order = def.order ?? Number.POSITIVE_INFINITY;
      let index = views.length;
      for (let i = 0; i < views.length; i++) {
        const o = this.view(views[i])?.order ?? Number.POSITIVE_INFINITY;
        if (o > order) {
          index = i;
          break;
        }
      }
      l = moveView(layout, ref, { kind: 'group', id: home, index });
      // Placing is not opening: the group keeps showing what it showed, folded or not.
      if (!activate) l = setGroupState(l, home, place.group.active ?? ref, place.group.collapsed === true);
    } else {
      // The default layout's empty right group is not kept (normalized away): it is made here, at the option's width.
      l = moveView(layout, ref, { kind: 'edge', of: 'main', side: 'right', id: home, size: this.#size });
      if (!activate) l = setCollapsed(l, home, true);
    }
    return l;
  }

  get sectionsAll(): ReadonlyState<readonly ShellSection[]> {
    return this.#sections.asReadonly;
  }

  /** The sections of a slot, by order. */
  sections(slot: string): ReadonlyState<readonly ShellSection[]> {
    return derive(this.#sections, (all) =>
      all.filter((s) => s.slot === slot).sort((a, b) => (a.order ?? 0) - (b.order ?? 0)),
    );
  }

  registerSections(sections: readonly ShellSection[]): () => void {
    this.#sections.set([...this.#sections.value, ...sections]);
    return () => this.#sections.set(this.#sections.value.filter((s) => !sections.includes(s)));
  }

  registerRailItems(items: readonly RailItem[]): () => void {
    this.#railItems.set([...this.#railItems.value, ...items]);
    return () => this.#railItems.set(this.#railItems.value.filter((i) => !items.includes(i)));
  }

  /** Whether a section is open, as the person left it this session (`undefined`: not touched yet). */
  sectionOpen(view: string, section: string): boolean | undefined {
    return this.#sectionOpen.get(`${view}/${section}`);
  }

  setSectionOpen(view: string, section: string, open: boolean): void {
    this.#sectionOpen.set(`${view}/${section}`, open);
  }

  // ---------------------------------------------------------------------------------------------------------------
  // The layout

  get layout(): ReadonlyState<ShellLayout> {
    return this.#layout.asReadonly;
  }

  /**
   * The layout as a value to keep (shell §2, persistence): views whose geometry is not to be remembered
   * (`remember: false`) and instances whose key the app maps back on restore are left as they are; the app stores it.
   */
  save(): ShellLayout {
    let l = this.#layout.value;
    for (const ref of allGroups(l).flatMap((g) => g.views)) {
      if (this.view(ref)?.remember === false) l = removeView(l, ref);
    }
    return l;
  }

  /**
   * A stored layout back (validated: another version or a wrong shape is ignored, and the current one stays); the views
   * it does not have go to their homes. Returns whether it was taken.
   */
  restore(value: unknown): boolean {
    const l = validateLayout(value);
    if (!l) return false;
    this.#layout.set(this.#placeMissing(this.#clampWindows(l), [...this.#views.value.values()]));
    return true;
  }

  /** Back to the default layout: every view at its home, nothing hidden, windows closed, the rail shown. */
  resetLayout(): void {
    batch(() => {
      this.#lastRects.clear();
      const l = this.#placeMissing(this.#defaultLayout, [...this.#views.value.values()]);
      // Embedded keeps its rail hidden (the host page drives it).
      this.#layout.set(this.#modeOption === 'embedded' ? setRail(l, false) : l);
    });
  }

  // ---------------------------------------------------------------------------------------------------------------
  // What a person does

  /**
   * Shows a view (shell §2): a card, a page or a modal opens in its layer; a placed view becomes its group's active one
   * (the group opens, a window comes to the top, a hidden view comes back); a view not placed yet goes to its home, or
   * to `target`. With `section`, that section is brought into view.
   */
  show(id: string, options: ShowOptions = {}): void {
    const def = this.view(id);
    if (!def) return;
    const ref = viewRef(def.id, def.instances === 'many' ? (options.key ?? 'default') : undefined);
    batch(() => {
      if ('subject' in options) this.#subjectOf(ref).set(options.subject);
      if (options.section !== undefined) this.#focusOf(ref).set(options.section);
      if (!options.target && isLayerHome(def.home)) {
        const open: OpenLayer = { ref, def, anchor: options.anchor };
        if (def.home === 'card') {
          this.#page.set(undefined);
          this.#card.set(open);
        } else if (def.home === 'page') {
          this.#card.set(undefined);
          this.#page.set(open);
        } else this.#modal.set(open);
        return;
      }
      // Opening a docked view leaves a card (it shows in the card's place).
      if (this.#card.value && !options.target) this.#card.set(undefined);
      let l = this.#layout.value;
      const place = findView(l, ref);
      if (place && !options.target) {
        l = selectView(l, ref);
        if (place.window) l = raiseWindow(l, place.window.group.id);
      } else if (options.target) {
        l = moveView(l, ref, options.target);
      } else if (def.home === 'float') {
        l = moveView(l, ref, { kind: 'float', rect: this.#newWindowRect(ref, def, options.anchorRect) });
      } else {
        l = this.#placeAtHome(l, ref, def, true);
      }
      this.#layout.set(l.rail || !railGroup(l)?.views.includes(ref) ? l : setRail(l, true));
    });
  }

  /**
   * A press on a view's button in its strip: shows it, or, when it is already the open one of an open group, folds the
   * group (a rail's panel closes to the rail).
   */
  select(ref: ViewRef): void {
    const place = findView(this.#layout.value, ref);
    if (!place) return this.show(ref);
    if (place.group.active === ref && !place.group.collapsed && !place.window) {
      this.#layout.set(setCollapsed(this.#layout.value, place.group.id, true));
      return;
    }
    this.show(viewIdOf(ref), { key: keyOf(ref) });
  }

  /** Folds a group to its strip, or opens it. */
  collapse(groupId: string, collapsed = true): void {
    this.#layout.set(setCollapsed(this.#layout.value, groupId, collapsed));
  }

  /** Closes a view where it is: a card, page or modal leaves; a window's view and an instance leave the layout; a docked single view is hidden. */
  close(ref: ViewRef): void {
    batch(() => {
      for (const layer of [this.#card, this.#page, this.#modal]) {
        if (layer.value?.ref === ref) {
          layer.set(undefined);
          return;
        }
      }
      const l = this.#layout.value;
      const place = findView(l, ref);
      if (!place) return;
      const def = this.view(ref);
      if (place.window) this.#lastRects.set(ref, place.window.rect);
      if (!place.window && (def?.instances ?? 'single') === 'single' && def?.home !== 'float') {
        this.#layout.set(setHidden(l, ref, true));
        return;
      }
      if (place.window && def && (def.instances ?? 'single') === 'single' && def.home !== 'float') {
        // A docked view popped out goes back to its home, hidden (it was closed, not docked).
        this.#layout.set(setHidden(this.#placeAtHome(removeView(l, ref), ref, def, false), ref, true));
        return;
      }
      this.#layout.set(removeView(l, ref));
      this.#subjects.delete(ref);
      this.#focus.delete(ref);
    });
  }

  /** Moves a view (shell §2: pop out, dock, split, into another group). */
  move(ref: ViewRef, target: LayoutTarget): void {
    this.#layout.set(moveView(this.#layout.value, ref, target));
  }

  /** Pops a view out into a window of its own (beside `anchorRect`, at its last place, or cascaded). */
  float(ref: ViewRef, anchorRect?: Rect): void {
    const def = this.view(ref);
    if (!def) return;
    this.move(ref, { kind: 'float', rect: this.#newWindowRect(ref, def, anchorRect) });
  }

  /** A view back to its home group (from a window, from another group). */
  dock(ref: ViewRef): void {
    const def = this.view(ref);
    if (!def) return;
    const place = findView(this.#layout.value, ref);
    if (place?.window) this.#lastRects.set(ref, place.window.rect);
    this.#layout.set(this.#placeAtHome(removeView(this.#layout.value, ref), ref, def, true));
  }

  /** A window's views into another group (a window docked whole). */
  merge(fromGroup: string, intoGroup: string): void {
    this.#layout.set(mergeGroup(this.#layout.value, fromGroup, intoGroup));
  }

  /**
   * Closes a docked group that is not the rail's (its ✕): its views go back into the rail's group, which keeps showing
   * what it showed, folded or not. Why: closing a split must not open anything, nor leave views out of reach.
   */
  closeGroup(groupId: string): void {
    const l = this.#layout.value;
    const rail = railGroup(l);
    const from = findGroup(l, groupId);
    if (!rail || !from || rail.id === groupId || from.window) return;
    let next = l;
    for (const v of from.group.views) next = moveView(next, v, { kind: 'group', id: rail.id });
    this.#layout.set(setGroupState(next, rail.id, rail.active, rail.collapsed === true));
  }

  /** A docked view into a new group under its own (Split below). */
  splitBelow(ref: ViewRef): void {
    const place = findView(this.#layout.value, ref);
    if (!place || place.window || place.group.views.length < 2) return;
    this.move(ref, { kind: 'edge', of: place.group.id, side: 'bottom', size: 240 });
  }

  /** A docked view into the group above or below its own in its column. */
  moveInColumn(ref: ViewRef, dir: -1 | 1): void {
    this.#layout.set(moveInColumn(this.#layout.value, ref, dir));
  }

  resize(nodeId: string, size: number): void {
    this.#layout.set(resizeNode(this.#layout.value, nodeId, size));
  }

  setWindowRect(groupId: string, rect: Rect): void {
    const area = this.#area.value;
    // At least its views' least size (`minSize`).
    for (const v of findGroup(this.#layout.value, groupId)?.group.views ?? []) {
      const min = this.view(v)?.minSize;
      if (min) rect = { ...rect, width: Math.max(rect.width, min.width), height: Math.max(rect.height, min.height) };
    }
    this.#layout.set(setWindowRect(this.#layout.value, groupId, area.width ? clampRect(rect, area) : rect));
  }

  raise(groupId: string): void {
    this.#layout.set(raiseWindow(this.#layout.value, groupId));
  }

  /** Out of view, its place kept (a view a person never uses). */
  hide(ref: ViewRef): void {
    this.#layout.set(setHidden(this.#layout.value, ref, true));
  }

  unhide(ref: ViewRef): void {
    this.#layout.set(setHidden(this.#layout.value, ref, false));
  }

  get hiddenViews(): ReadonlyState<readonly ViewRef[]> {
    return derive(this.#layout, (l) => l.hidden);
  }

  /** Shows or hides the rail (and the docked column with it). */
  toggleRail(shown = !this.#layout.value.rail): void {
    this.#layout.set(setRail(this.#layout.value, shown));
  }

  get railShown(): ReadonlyState<boolean> {
    return derive(this.#layout, (l) => l.rail);
  }

  // ---------------------------------------------------------------------------------------------------------------
  // Layers

  get card(): ReadonlyState<OpenLayer | undefined> {
    return this.#card.asReadonly;
  }

  get page(): ReadonlyState<OpenLayer | undefined> {
    return this.#page.asReadonly;
  }

  get modal(): ReadonlyState<OpenLayer | undefined> {
    return this.#modal.asReadonly;
  }

  closeCard(): void {
    this.#card.set(undefined);
  }

  closePage(): void {
    this.#page.set(undefined);
  }

  closeModal(): void {
    this.#modal.set(undefined);
  }

  // ---------------------------------------------------------------------------------------------------------------
  // A view's own state

  /** What an instance is about (`show`'s `subject`). */
  subject(ref: ViewRef): ReadonlyState<unknown> {
    return this.#subjectOf(ref).asReadonly;
  }

  /** The section a `show` asked for (a new value each time it asks). */
  focus(ref: ViewRef): ReadonlyState<string | undefined> {
    return this.#focusOf(ref).asReadonly;
  }

  #subjectOf(ref: ViewRef): State<unknown> {
    let s = this.#subjects.get(ref);
    if (!s) this.#subjects.set(ref, (s = new State<unknown>(undefined)));
    return s;
  }

  #focusOf(ref: ViewRef): State<string | undefined> {
    let s = this.#focus.get(ref);
    if (!s) this.#focus.set(ref, (s = new State<string | undefined>(undefined, () => false)));
    return s;
  }

  // ---------------------------------------------------------------------------------------------------------------
  // Size and mode

  /** `auto`'s mode follows this (the shell's own width, px). */
  setWidth(width: number): void {
    if (width !== this.#width.value) this.#width.set(width);
  }

  /** The main area's size: windows open and stay inside it. */
  setArea(width: number, height: number): void {
    const a = this.#area.value;
    if (a.width === width && a.height === height) return;
    batch(() => {
      this.#area.set({ width, height });
      this.#layout.set(this.#clampWindows(this.#layout.value));
    });
  }

  get area(): ReadonlyState<{ width: number; height: number }> {
    return this.#area.asReadonly;
  }

  get mode(): ReadonlyState<ShellMode> {
    return this.#mode;
  }

  #clampWindows(l: ShellLayout): ShellLayout {
    const area = this.#area.value;
    if (!area.width || !area.height) return l;
    let next = l;
    for (const w of l.windows) {
      const r = clampRect(w.rect, area);
      if (r.x !== w.rect.x || r.y !== w.rect.y || r.width !== w.rect.width || r.height !== w.rect.height) {
        next = setWindowRect(next, w.group.id, r);
      }
    }
    return next;
  }

  #newWindowRect(ref: ViewRef, def: ShellViewDef, anchorRect?: Rect): Rect {
    const area = this.#area.value.width ? this.#area.value : { width: 1200, height: 800 };
    const last = def.remember === false ? undefined : this.#lastRects.get(ref);
    if (last) return clampRect(last, area);
    const embedded = this.#mode.value === 'embedded';
    const size = def.size ?? DefaultWindow;
    const scaled = embedded ? { width: Math.round(size.width * 0.8), height: Math.round(size.height * 0.8) } : size;
    return placeWindow(
      scaled,
      area,
      stacked(this.#layout.value).map((w) => w.rect),
      anchorRect,
    );
  }
}

function keyOf(ref: ViewRef): string | undefined {
  const i = ref.indexOf('#');
  return i < 0 ? undefined : ref.slice(i + 1);
}
