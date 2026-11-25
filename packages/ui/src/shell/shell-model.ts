/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { ContextMenuItem } from '../components/context-menu/context-menu.js';
import { batch, derive, listen, type ReadonlyState, State } from '@molstar/core/reactive';
import type { Child } from '../dom/index.js';
import type { ViewRenderer } from '../dom/view.js';
import type { IconName } from '../components/icon/icon.js';
import {
  allGroups,
  arrangeGrid,
  clampRect,
  DefaultColumnSize,
  defaultLayout,
  DefaultSplitSize,
  DefaultTopSize,
  DefaultWindowSize,
  findGroup,
  findView,
  type GroupNode,
  gridPeers,
  type LayoutIndex,
  type LayoutNode,
  inRegion,
  inRightColumn,
  keepHidden,
  layoutIndex,
  mainNode,
  MinGroupSize,
  MinMainSize,
  MinPaneSize,
  RailWidth,
  type SplitNode,
  type LayoutTarget,
  mergeGroup,
  moveInColumn,
  moveView,
  placeWindow,
  type Rect,
  RightGroup,
  railGroup,
  raiseWindow,
  regionHas,
  removeView,
  resizeNode,
  selectView,
  setCollapsed,
  setGroupState,
  setHidden,
  setRail,
  setSizes,
  setWindowRect,
  type ShellLayout,
  stacked,
  validateLayout,
  type ViewRef,
} from './layout.js';

/**
 * Where a view shows (shell §3; it may look different in each): `stage`, a pane of the main area's region; `top`,
 * overlaid on the content's top (shell §2.1).
 */
export type ViewPlacement = 'docked' | 'window' | 'card' | 'page' | 'modal' | 'sheet' | 'stage' | 'top';

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
  /** Closes it where it is (`ShellModel.close`: a layer's view leaves it; a docked single view is hidden). */
  close(): void;
  /** Its pane's port (shell §2.1): what the layout tells the pane, and the little the pane tells the layout. */
  readonly pane: PanePort;
}

/**
 * Where a view goes when it is first shown (shell §2): a docked group, a window, one of the shell's layers, or the main
 * area's region (`stage`: beside the app's content; `top`: on the content's top, overlaid or split by the shell's `top`
 * option; shell §2.1).
 */
export type ViewHome = string | 'float' | 'card' | 'page' | 'modal' | 'stage' | 'top';

/** A side of a node in the main area's region. */
export type StageSide = 'left' | 'right' | 'top' | 'bottom';

/**
 * A pane's status as the layout knows it (shell §2.1), for the pane's own content to react to (a stage shows its
 * toolbar while active or hovered): the layout gives the status, the panel decides what it shows.
 */
export interface PaneStatus {
  /** The pane a person or an agent works in (`activePane`): one at a time (the app's content when nothing was chosen). */
  readonly active: boolean;
  /** The pointer is over it. */
  readonly hovered: boolean;
  /** The focus is in it. */
  readonly focused: boolean;
  /** It is the one maximized pane. */
  readonly maximized: boolean;
  /** Several panes show (not one alone, nor one maximized): one on its own needs no mark of being the active one. */
  readonly several: boolean;
}

/**
 * What of a pane's bottom edge the shell's own chrome covers (px from each corner; shell §2.1): `ShellView`'s `over`
 * sits over the main area's bottom corners, so the pane in a corner keeps its own bottom bar clear of it.
 */
export interface PaneInsets {
  readonly bottomLeft: number;
  readonly bottomRight: number;
}

const NoInsets: PaneInsets = { bottomLeft: 0, bottomRight: 0 };

/**
 * A pane's port (shell §2.1): the layout's side of a pane. The layout tells the pane where it stands (`status`, `insets`)
 * and does to it what a person does (maximize, float, close); the pane tells the layout only its title. What a pane
 * holds (a stage, a plugin) is the app's: the layout never stores it (`PaneMap` keeps it, app side).
 */
export interface PanePort {
  /** The pane's id: `content` (the app's main view) or a view's ref (one kind of id in every pane API, shell §2.1). */
  readonly id: string;
  readonly status: ReadonlyState<PaneStatus>;
  readonly insets: ReadonlyState<PaneInsets>;
  /** Its title (an instance's own: "Stage D"): tabs, a window's header, a frame's name. Returns the undo. */
  setTitle(title: string): () => void;
  /**
   * Shows it alone in the main area, or every pane again (shell §2.1): a pane of the region alone, a window's out of its
   * window over the main area; nothing for a view that is no pane (docked beside the main area). The shell's chrome
   * stays (that is `toggleExpanded`, the shell's).
   */
  toggleMaximized(): void;
  float(): void;
  close(): void;
}

/**
 * Where a `top` view goes (shell §2.1): overlaid on the content's top (the app places it: a stage's band slot), above
 * the content in a split of its own, or either (both offered).
 */
export type TopMode = 'overlay' | 'split' | 'both';

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
  /** Its size as a window (default `DefaultWindowSize`, 360 × 320). */
  readonly size?: { readonly width: number; readonly height: number };
  /** The least size of a window holding it (every window rectangle the shell sets, `MinWindow` at least). */
  readonly minSize?: { readonly width: number; readonly height: number };
  /** Where it goes when first shown (default the right group). */
  readonly home?: ViewHome;
  /**
   * Its place in the main area's region, for `home: 'stage'` (`side` of the app's content, default right; `size` px,
   * default `DefaultGroupSize`) and `home: 'top'` in a split (`size`, default `DefaultTopSize`; shell §2.1).
   */
  readonly region?: { readonly side?: StageSide; readonly size?: number };
  /** `many`: the app opens instances by a key, each its own place and state (two plots side by side). */
  readonly instances?: 'single' | 'many';
  /** False: its geometry is not part of a saved layout (it opens at its default place and size each time). */
  readonly remember?: boolean;
  /** A card that stays on a press outside it (Escape or its ✕ close it). */
  readonly persistent?: boolean;
  /** Its content takes the body's whole height and scrolls its own parts (the body does not scroll). */
  readonly fill?: boolean;
  /**
   * A canvas, not a panel (shell §2.1): alone in its pane (in the main area's region, or maximized out of a window) it
   * has no header (a lone canvas has no chrome above it; its menu is a right-click on its tab when it has one).
   */
  readonly canvas?: boolean;
}

/** A section contributed to a slot of a view (an Export view's `export` slot): extensions add sections, not views. */
export interface ShellSection {
  readonly id: string;
  readonly slot: string;
  readonly title: string;
  /** Its header's words when they follow what it offers. */
  readonly heading?: ReadonlyState<string>;
  readonly order?: number;
  /** Shown only while this holds: asked again when the view's subject changes, and when `whenChanged` changes. */
  readonly when?: (ctx: ShellViewContext) => boolean;
  /** Changes when `when` may answer otherwise for reasons the list cannot see (a value fetched for the subject). */
  readonly whenChanged?: ReadonlyState<unknown>;
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
  /**
   * In the shell's bottom-right corner: the rail's last button, under its Hide button (the command menu; shell §3.2).
   * Its `group` is ignored; it stays when the rail runs short.
   */
  readonly corner?: boolean;
  readonly order?: number;
  readonly priority?: number;
  /** A small caption under its icon (a shortcut). */
  readonly caption?: string;
  readonly pressed?: ReadonlyState<boolean>;
  readonly visible?: ReadonlyState<boolean>;
  /**
   * What a press does; `anchor` is the button (to place what it opens). An item with a `menu` and no `run` opens the
   * menu on a press (not on hover: the pointer travels along the rail, and the hover is the tooltip's).
   */
  run?(anchor: Element): void;
  /** Its menu, read as it opens (a choice's rows mark the chosen one with `checked`): beside the button, on a press. */
  menu?(): readonly ContextMenuItem[];
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
  /** The right group's width when it is made (default `DefaultColumnSize`, 340). */
  readonly size?: number;
  /** Where a `top` view goes (shell §2.1): `split` (default; above the content), `overlay` (on its top), or `both`. */
  readonly top?: TopMode;
  /** The shell's words (shell §6): the app's replace the defaults (`DefaultShellLabels`, English) by name. */
  readonly labels?: Partial<ShellLabels>;
}

/**
 * Every word the shell shows or says (shell §6): its chrome, its menus, its commands. Strings, or functions for the
 * parametrized ones (a view's title in them). The app passes its own (`ShellOptions.labels`), merged once in the
 * model (`shell.labels`), so the view, the menus and the commands read the same words, headless too.
 */
export interface ShellLabels {
  /** The shell's root: its default name, and what it is (`aria-roledescription`). */
  readonly workspace: string;
  readonly workspaceRole: string;
  readonly mainArea: string;
  /** The app's content as a pane, when it gave itself no title. */
  readonly mainTitle: string;
  /** A splitter. */
  readonly resize: string;
  /** The rail (a `nav`). */
  readonly rail: string;
  /** The rail's lists of views (its head, its strip, its foot: tab lists). */
  readonly railViews: string;
  /** The rail's groups of tools (the app's rail items, by group; its foot's; its corner's). */
  readonly railTools: string;
  readonly more: string;
  /** The rail's Hide button, and the button that shows it again. */
  readonly hideViews: string;
  readonly showViews: string;
  /** A window none of whose views shows. */
  readonly window: string;
  /** The busy overlay's button. */
  readonly cancel: string;
  /** The compact layout's bar. */
  readonly compactBar: string;
  /** The file drop's overlay (`ShellView`'s `dropFormats` goes under it). */
  readonly dropFiles: string;
  /** A frame's ✕ and its ⋯ button. */
  readonly close: (title: string) => string;
  readonly layoutMenu: (title: string) => string;
  /** Why fullscreen cannot be asked for (`fullscreenAvailable`): headless, no API, not allowed by the page. */
  readonly fullscreenUnavailable: string;
  readonly fullscreenNoApi: string;
  readonly fullscreenNotAllowed: string;
  // The menus (`viewMenu`, `paneMenu`, `railMenu`)
  readonly dock: string;
  readonly moveWithKeyboard: string;
  readonly resizeWithKeyboard: string;
  readonly resetPosition: string;
  readonly float: string;
  readonly splitBelow: string;
  readonly splitRight: string;
  readonly moveUp: string;
  readonly moveDown: string;
  readonly dockRight: string;
  /** Into an open window, named by its view's title. */
  readonly moveToWindow: (title: string) => string;
  readonly hide: string;
  readonly showView: (title: string) => string;
  readonly hideRail: string;
  readonly resetLayout: string;
  readonly besideMain: string;
  readonly belowMain: string;
  readonly aboveMain: string;
  readonly overMainTop: string;
  readonly maximize: string;
  readonly showEveryPane: string;
  readonly arrangeGrid: string;
  readonly closePane: string;
  // The commands (`shellCommands`): titles, then what they say back
  /** Their category (default; `shellCommands`' own argument replaces it). */
  readonly commandCategory: string;
  readonly expandTitle: string;
  readonly expandedMessage: string;
  readonly notExpandedMessage: string;
  readonly fullscreenTitle: string;
  readonly toggleRailTitle: string;
  readonly viewsShownMessage: string;
  readonly viewsHiddenMessage: string;
  readonly maximizePaneTitle: string;
  readonly maximizedMessage: string;
  readonly allPanesMessage: string;
  readonly onePane: string;
  readonly gridTitle: string;
  readonly gridMessage: string;
  readonly closePaneTitle: string;
  readonly nothingToClose: string;
  readonly paneClosedMessage: string;
  readonly resetLayoutTitle: string;
  readonly layoutResetMessage: string;
}

/** The shell's words in English (shell §6). */
export const DefaultShellLabels: ShellLabels = {
  workspace: 'Workspace',
  workspaceRole: 'workspace',
  mainArea: 'Main area',
  mainTitle: 'Main',
  resize: 'Resize',
  rail: 'Views and tools',
  railViews: 'Views',
  railTools: 'Tools',
  more: 'More',
  hideViews: 'Hide the views',
  showViews: 'Show the views',
  window: 'Window',
  cancel: 'Cancel',
  compactBar: 'Views',
  dropFiles: 'Drop files to open them',
  close: (title) => `Close ${title}`,
  layoutMenu: (title) => `${title}: layout`,
  fullscreenUnavailable: 'Fullscreen is not available here',
  fullscreenNoApi: 'This browser has no fullscreen',
  fullscreenNotAllowed: 'Fullscreen is not allowed here (the page embedding this does not allow it)',
  dock: 'Dock',
  moveWithKeyboard: 'Move with the keyboard',
  resizeWithKeyboard: 'Resize with the keyboard',
  resetPosition: 'Reset position',
  float: 'Float',
  splitBelow: 'Split below',
  splitRight: 'Split right',
  moveUp: 'Move up',
  moveDown: 'Move down',
  dockRight: 'Dock right',
  moveToWindow: (title) => `Move to the ${title} window`,
  hide: 'Hide',
  showView: (title) => `Show ${title}`,
  hideRail: 'Hide the rail',
  resetLayout: 'Reset layout',
  besideMain: 'Beside the main view',
  belowMain: 'Below the main view',
  aboveMain: 'Above the main view',
  overMainTop: 'Over the main view’s top',
  maximize: 'Maximize',
  showEveryPane: 'Show every pane',
  arrangeGrid: 'Arrange as a grid',
  closePane: 'Close',
  commandCategory: 'View',
  expandTitle: 'Expand the main area',
  expandedMessage: 'Expanded: Escape to come back',
  notExpandedMessage: 'Back from expanded',
  fullscreenTitle: 'Fullscreen',
  toggleRailTitle: 'Show or hide the views',
  viewsShownMessage: 'Views shown',
  viewsHiddenMessage: 'Views hidden',
  maximizePaneTitle: 'Maximize the pane',
  maximizedMessage: 'Maximized: Escape to come back',
  allPanesMessage: 'All panes',
  onePane: 'There is one pane',
  gridTitle: 'Arrange the panes as a grid',
  gridMessage: 'Panes arranged as a grid',
  closePaneTitle: 'Close the pane',
  nothingToClose: 'Nothing to close',
  paneClosedMessage: 'Pane closed',
  resetLayoutTitle: 'Reset the layout',
  layoutResetMessage: 'Layout reset',
};

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

/** The main area's size while it is not measured (headless, before the first layout): where new windows are placed. */
export const UnmeasuredArea = { width: 1200, height: 800 } as const;

/** A view's ref: its id, or `id#key`. */
export function viewRef(id: string, key?: string): ViewRef {
  return key === undefined ? id : `${id}#${key}`;
}

/** The view id of a ref. */
export function viewIdOf(ref: ViewRef): string {
  const i = ref.indexOf('#');
  return i < 0 ? ref : ref.slice(0, i);
}

/** The instance key of a ref (`id#key`), if it has one. */
export function keyOf(ref: ViewRef): string | undefined {
  const i = ref.indexOf('#');
  return i < 0 ? undefined : ref.slice(i + 1);
}

/** A number for each model (`ShellModel.uid`: ids of its elements unique on a page with several shells). */
let nextShellId = 1;

const isLayerHome = (home: ViewHome | undefined) => home === 'card' || home === 'page' || home === 'modal';
/** A home that places a view when it is shown, not when it is registered: a window, the main area's region. */
const isShownHome = (home: ViewHome | undefined) => home === 'float' || home === 'stage' || home === 'top';

/**
 * The shell's state (shell §2), DOM-free: what is registered (views, sections, rail items), the layout (a `State` of
 * a `ShellLayout`, the value an app saves and restores), what is open in the layers over it (a card, a page, a modal),
 * the mode. Operations are named for what a person does (show, select, move, float, dock, close); each is one
 * change of the layout.
 */
export class ShellModel {
  readonly #views = new State<ReadonlyMap<string, ShellViewDef>>(new Map());
  readonly #sections = new State<readonly ShellSection[]>([]);
  readonly #railItems = new State<readonly RailItem[]>([]);
  /** Views whose `visible` is false now: out of their strips, their places kept. */
  readonly #absent = new State<ReadonlySet<string>>(new Set());
  /** The listener of each registered view's `visible` flag (one a view id: registering anew replaces it). */
  readonly #visibleStops = new Map<string, () => void>();
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
  readonly #slots = new Map<string, ReadonlyState<readonly ShellSection[]>>();
  /** Expanded (shell §3.9): the shell's chrome away, the main area takes the shell. */
  readonly #expanded = new State(false);
  /** The pane shown alone (its id: `content` or a view's ref, shell §2.1). */
  readonly #maximized = new State<string | undefined>(undefined);
  /** A grid asked for before the main area had a size (`arrangeGrid`). */
  #gridPending: { columns?: number } | undefined;
  readonly #activePane = new State('content');
  readonly #hoveredPane = new State<string | undefined>(undefined);
  readonly #focusedPane = new State<string | undefined>(undefined);
  readonly #paneStatus = new Map<string, ReadonlyState<PaneStatus>>();
  readonly #titles = new State<ReadonlyMap<string, string>>(new Map());
  readonly #paneTitles = new Map<string, ReadonlyState<string>>();
  readonly #ports = new Map<string, PanePort>();
  readonly #paneInsets = new Map<string, ReadonlyState<PaneInsets>>();
  /** The widths of the shell's chrome over the main area's bottom corners (`over`), the view's input. */
  readonly #overChrome = new State<PaneInsets>(NoInsets);
  readonly #top: TopMode;
  readonly #fullscreen = new State(false);
  readonly #fullscreenRequest = new State<boolean | undefined>(undefined);
  readonly #fullscreenAvailable: State<true | string>;
  readonly #hiddenViews: ReadonlyState<readonly ViewRef[]>;
  readonly #railShown: ReadonlyState<boolean>;
  /** Each view's last window rectangle this session (a window closed and opened again comes back there). */
  readonly #lastRects = new Map<ViewRef, Rect>();
  /** Several panes show (in the main area's region, none maximized): every pane's `several`. */
  readonly #several: ReadonlyState<boolean>;
  /** The panes in the main area's bottom corners (or the one shown alone): every pane's insets read it. */
  readonly #corners: ReadonlyState<{ readonly left?: string; readonly right?: string }>;
  /** The shell's words (shell §6): the app's (`ShellOptions.labels`) over `DefaultShellLabels`, merged once. */
  readonly labels: ShellLabels;
  /** Unique among the models on a page: the view's element ids start with it (two shells, no clash). */
  readonly uid = `msui-shell-${nextShellId++}`;

  constructor(options: ShellOptions = {}) {
    this.labels = Object.freeze({ ...DefaultShellLabels, ...options.labels });
    this.#fullscreenAvailable = new State<true | string>(this.labels.fullscreenUnavailable);
    this.#modeOption = options.mode ?? 'auto';
    this.#compactWidth = options.compactWidth ?? 720;
    this.#size = options.size ?? DefaultColumnSize;
    this.#top = options.top ?? 'split';
    this.#defaultLayout = defaultLayout();
    this.#layout = new State(validateLayout(options.layout) ?? this.#defaultLayout);
    this.#mode = derive(this.#width, (w) =>
      this.#modeOption !== 'auto' ? this.#modeOption : w > 0 && w < this.#compactWidth ? 'compact' : 'desktop',
    );
    if (this.#modeOption === 'embedded') this.#layout.set(setRail(this.#layout.value, false));
    // A new list only when it changes (every layout write makes a new array).
    this.#hiddenViews = derive(
      this.#layout,
      (l) => l.hidden,
      (a, b) => a.length === b.length && a.every((v, i) => v === b[i]),
    );
    this.#railShown = derive(this.#layout, (l) => l.rail);
    // Of the panes that show (a group whose views are all flagged off or hidden is none), from the region's tree: a
    // layout write that keeps it (a window's move, the right column's drag) asks nothing again (`normalize` keeps it).
    const region = derive(this.#layout, (l) => mainNode(l).root);
    const sources = [region, this.#maximized, this.#absent, this.#views, this.#hiddenViews] as const;
    this.#several = derive(sources, (root, m) => !m && this.#regionPanes(root, 2) >= 2);
    this.#corners = derive(
      sources,
      // Alone (maximized, in the region or out of its window): every corner is its.
      (root, maximized) => (maximized ? { left: maximized, right: maximized } : this.#cornerPanes(root)),
      (a, b) => a.left === b.left && a.right === b.right,
    );
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
   * for `show`. A view registered again under its id replaces the one before (its `visible` flag followed instead);
   * the undo of the one before then leaves it.
   */
  registerViews(defs: readonly ShellViewDef[]): () => void {
    batch(() => {
      const next = new Map(this.#views.value);
      for (const d of defs) {
        next.set(d.id, d);
        this.#followVisible(d.id, d.visible);
      }
      this.#views.set(next);
      this.#setLayout(this.#placeMissing(this.#layout.value, defs));
    });
    return () =>
      batch(() => {
        const next = new Map(this.#views.value);
        const removed = new Set<string>();
        for (const d of defs) {
          // Registered anew since: that one stays.
          if (next.get(d.id) !== d) continue;
          next.delete(d.id);
          removed.add(d.id);
          this.#followVisible(d.id, undefined);
        }
        if (!removed.size) return;
        this.#views.set(next);
        // Its place stays in the layout (a view that registers again lands there); what showed it in a layer goes.
        for (const layer of [this.#card, this.#page, this.#modal]) {
          if (layer.value && removed.has(viewIdOf(layer.value.ref))) layer.set(undefined);
        }
      });
  }

  /** Follows a view's `visible` flag (none: present), the listener of the view registered before under its id stopped. */
  #followVisible(id: string, flag: ReadonlyState<boolean> | undefined): void {
    this.#visibleStops.get(id)?.();
    this.#visibleStops.delete(id);
    if (!flag) return this.#setAbsent(id, false);
    const write = () => this.#setAbsent(id, !flag.value);
    write();
    this.#visibleStops.set(id, listen(flag, write));
  }

  /** Every single view with a docked home that the layout does not have, placed at its home by its order. */
  #placeMissing(layout: ShellLayout, defs: readonly ShellViewDef[]): ShellLayout {
    let l = layout;
    for (const d of defs) {
      if ((d.instances ?? 'single') !== 'single' || isLayerHome(d.home) || isShownHome(d.home)) continue;
      if (findView(l, d.id)) continue;
      l = this.#placeAtHome(l, d.id, d, false);
    }
    return l;
  }

  /**
   * `ref` into a docked group (`into`, default its def's home: a group, else the right group; a `stage` or `top` home is
   * the main area's region), at its order; a group that is not there is made beside the main area, on its right.
   */
  #placeAtHome(layout: ShellLayout, ref: ViewRef, def: ShellViewDef, activate: boolean, into?: string): ShellLayout {
    if (into === undefined && (def.home === 'stage' || def.home === 'top')) return this.#placeInStage(layout, ref, def);
    const home = into ?? (def.home && !isLayerHome(def.home) && !isShownHome(def.home) ? def.home : RightGroup);
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
      // Made here (the default layout has no right group), at the option's width.
      l = moveView(layout, ref, { kind: 'edge', of: 'main', side: 'right', id: home, size: this.#size });
      if (!activate) l = setCollapsed(l, home, true);
    }
    return l;
  }

  /** A view with a home in the main area's region, placed there (shell §2.1): beside the content, or on its top. */
  #placeInStage(layout: ShellLayout, ref: ViewRef, def: ShellViewDef): ShellLayout {
    const size = def.region?.size;
    if (def.home === 'top') {
      return this.#top === 'split'
        ? moveView(layout, ref, { kind: 'edge', of: 'content', side: 'top', size: size ?? DefaultTopSize })
        : moveView(layout, ref, { kind: 'top' });
    }
    return moveView(layout, ref, { kind: 'edge', of: 'content', side: def.region?.side ?? 'right', size });
  }

  /** The sections of a slot, by order. */
  sections(slot: string): ReadonlyState<readonly ShellSection[]> {
    // One derived state a slot, made once (every list of the slot shares it).
    let d = this.#slots.get(slot);
    if (!d) {
      d = derive(this.#sections, (all) =>
        all.filter((s) => s.slot === slot).sort((a, b) => (a.order ?? 0) - (b.order ?? 0)),
      );
      this.#slots.set(slot, d);
    }
    return d;
  }

  registerSections(sections: readonly ShellSection[]): () => void {
    this.#sections.set([...this.#sections.value, ...sections]);
    return () => this.#sections.set(this.#sections.value.filter((s) => !sections.includes(s)));
  }

  /** Adds the app's rail items; their `visible` flags are followed as views' are (`absent`, `item:<id>`). */
  registerRailItems(items: readonly RailItem[]): () => void {
    const stops: (() => void)[] = [];
    batch(() => {
      for (const i of items) {
        if (!i.visible) continue;
        const flag = i.visible;
        const write = () => this.#setAbsent(`item:${i.id}`, !flag.value);
        write();
        stops.push(listen(flag, write));
      }
      this.#railItems.set([...this.#railItems.value, ...items]);
    });
    return () => {
      for (const stop of stops) stop();
      batch(() => {
        for (const i of items) if (i.visible) this.#setAbsent(`item:${i.id}`, false);
        this.#railItems.set(this.#railItems.value.filter((i) => !items.includes(i)));
      });
    };
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
   * Every layout write (shell §2, §2.1): the windows fitted (their views' least sizes, inside the main area), and what
   * refers to the layout follows it in the same batch. The maximized and the active pane follow what their region
   * group shows (a tab switch); a maximize ends when its pane is no longer one (closed, docked out of the main area)
   * and when `reveal` (a view shown or moved) lands in the region outside it (it would stay hidden behind the
   * maximize); the active pane falls back to `content` when its ref is no longer a pane. What was kept for a ref that
   * left the layout and every layer goes (its port's states, its title, its subject and section focus; a port an app
   * holds keeps working, and a new one is made on demand).
   */
  #setLayout(layout: ShellLayout, reveal?: ViewRef): void {
    const prev = this.#layout.value;
    const next = this.#fitWindows(layout);
    if (next === prev) return;
    batch(() => {
      this.#layout.set(next);
      const ix = layoutIndex(next);
      const m = this.#maximized.value;
      if (m !== undefined) {
        // Out of its window into the region, or out of the region into a window: no longer the pane maximized.
        const moved = m !== 'content' && !!layoutIndex(prev).views.get(m)?.window !== !!ix.views.get(m)?.window;
        this.#maximized.set(moved ? undefined : this.#followPane(ix, m, true, reveal));
      }
      const active = this.#activePane.value;
      if (active !== 'content') this.#activePane.set(this.#followPane(ix, active, false) ?? 'content');
      const before = layoutIndex(prev).views;
      for (const ref of before.keys()) {
        if (ix.views.has(ref) || this.#layered(ref)) continue;
        this.#forgetPane(ref);
        this.#forgetState(ref);
      }
    });
  }

  /**
   * A pane id after a layout write (shell §2.1): itself, the view its region group shows now (a tab switch), or
   * undefined when it is no pane (a pane: a view in a window or in the region, its top included; the content). The
   * maximized one (`maximized`) is a pane that can show alone: not the top group's, not the content without a region;
   * a view shown in the region outside it (`reveal`) ends it.
   */
  #followPane(ix: LayoutIndex, pane: string, maximized: boolean, reveal?: ViewRef): string | undefined {
    let group: string | undefined;
    let next = pane;
    if (pane === 'content') {
      if (maximized && !ix.main.root) return undefined;
    } else {
      const place = ix.views.get(pane);
      if (!place) return undefined;
      // A window's maximized view is held by its ref (out of its window); its active pane follows its tabs.
      if (place.window) return maximized ? pane : (place.group.active ?? pane);
      if (!(maximized ? ix.regionNodes : ix.regionGroups).has(place.group.id)) return undefined;
      group = place.group.id;
      next = place.group.active ?? pane;
    }
    if (reveal === undefined) return next;
    const shown = ix.views.get(reveal);
    if (!shown || shown.window || !ix.regionGroups.has(shown.group.id)) return next;
    return shown.group.id === group || (pane === 'content' && ix.main.top?.id === shown.group.id) ? next : undefined;
  }

  /**
   * The windows as a layout write leaves them (shell §2): each at least its views' least size (`minSize`), inside the
   * main area when it is measured. The same layout when every window is.
   */
  #fitWindows(l: ShellLayout): ShellLayout {
    const area = this.#area.value;
    let next = l;
    for (const w of l.windows) {
      let { width, height } = w.rect;
      for (const v of w.group.views) {
        const min = this.view(v)?.minSize;
        if (min) {
          width = Math.max(width, min.width);
          height = Math.max(height, min.height);
        }
      }
      let r: Rect = width === w.rect.width && height === w.rect.height ? w.rect : { ...w.rect, width, height };
      if (area.width && area.height) r = clampRect(r, area);
      const same = r.x === w.rect.x && r.y === w.rect.y && r.width === w.rect.width && r.height === w.rect.height;
      if (!same) next = setWindowRect(next, w.group.id, r);
    }
    return next;
  }

  /** Whether a ref is open in a layer (a card, a page, a modal). */
  #layered(ref: ViewRef): boolean {
    return this.#card.value?.ref === ref || this.#page.value?.ref === ref || this.#modal.value?.ref === ref;
  }

  /** What was kept for a pane whose ref left (its port and its states are made again on demand). */
  #forgetPane(ref: ViewRef): void {
    this.#paneStatus.delete(ref);
    this.#ports.delete(ref);
    this.#paneInsets.delete(ref);
    this.#paneTitles.delete(ref);
    if (this.#titles.value.has(ref)) {
      const next = new Map(this.#titles.value);
      next.delete(ref);
      this.#titles.set(next);
    }
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
    // A grid asked for before the main area was measured was the layout's being replaced.
    this.#gridPending = undefined;
    this.#setLayout(this.#clampSizes(this.#placeMissing(l, [...this.#views.value.values()])));
    return true;
  }

  /** Back to the default layout: every view at its home, nothing hidden, windows closed, the rail shown. */
  resetLayout(): void {
    batch(() => {
      this.#lastRects.clear();
      this.#gridPending = undefined;
      this.#maximized.set(undefined);
      const l = this.#placeMissing(this.#defaultLayout, [...this.#views.value.values()]);
      // Embedded keeps its rail hidden (the host page drives it).
      this.#setLayout(this.#modeOption === 'embedded' ? setRail(l, false) : l);
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
      // Opening a docked view (or a window) leaves Expand: what it opens must show.
      if (this.#expanded.value) this.#expanded.set(false);
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
      // Shown in the right column while the rail (and the column with it) is hidden: the rail comes back.
      const now = findView(l, ref);
      this.#setLayout(l.rail || !now || now.window || !inRightColumn(l, now.group.id) ? l : setRail(l, true), ref);
    });
  }

  /**
   * A press on a view's button in its strip: shows it, or, when it is already the open one of an open group docked
   * beside the main area, folds the group (a rail's panel closes to the rail). A pane's tab (in the main area's region,
   * on its top) and a window's never fold: a pane has no strip to fold to.
   */
  select(ref: ViewRef): void {
    const l = this.#layout.value;
    const place = findView(l, ref);
    // Not placed: shown as what it is (an instance keeps its key).
    if (!place) return this.show(viewIdOf(ref), { key: keyOf(ref) });
    if (place.group.active === ref && !place.group.collapsed && !place.window && !inRegion(l, place.group.id)) {
      this.#setLayout(setCollapsed(this.#layout.value, place.group.id, true));
      return;
    }
    this.show(viewIdOf(ref), { key: keyOf(ref) });
  }

  /** Folds a group to its strip, or opens it. */
  collapse(groupId: string, collapsed = true): void {
    this.#setLayout(setCollapsed(this.#layout.value, groupId, collapsed));
  }

  /**
   * Closes a view where it is (shell §2, close): a card, page or modal leaves its layer; a single view docked beside the
   * main area is hidden (its place kept); one with a docked home closed from a window goes back home, hidden. Any other
   * leaves the layout: a view in the main area's region or on its top, an instance, and a view whose home places it at
   * `show` (`float`, `stage`, `top`); shown again, it goes to its home.
   */
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
      // Its subject and focus go with it (`#setLayout`), and so does a maximize of its pane.
      const docked = !!def && (def.instances ?? 'single') === 'single' && !isShownHome(def.home);
      if (!docked || (!place.window && inRegion(l, place.group.id))) return this.#setLayout(removeView(l, ref));
      if (!place.window) return this.#setLayout(setHidden(l, ref, true));
      // A docked view popped out goes back to its home, hidden (it was closed, not docked).
      this.#setLayout(setHidden(this.#placeAtHome(removeView(l, ref), ref, def, false), ref, true));
    });
  }

  /** What was kept of a view that left the layout and every layer: its subject and section focus (shell §2, close). */
  #forgetState(ref: ViewRef): void {
    this.#subjects.delete(ref);
    this.#focus.delete(ref);
  }

  /** Moves a view (shell §2: pop out, dock, split, into another group). */
  move(ref: ViewRef, target: LayoutTarget): void {
    this.#setLayout(moveView(this.#layout.value, ref, target), ref);
  }

  /** Pops a view out into a window of its own (at its last place this session, beside `anchorRect`, or cascaded). */
  float(ref: ViewRef, anchorRect?: Rect): void {
    const def = this.view(ref);
    if (!def) return;
    this.move(ref, { kind: 'float', rect: this.#newWindowRect(ref, def, anchorRect) });
  }

  /**
   * A view docked (from a window, the main area's region, another group): into `into` (a docked group's id, made beside
   * the main area on its right when it is not there; Dock right is `RightGroup`), default its home.
   */
  dock(ref: ViewRef, into?: string): void {
    const def = this.view(ref);
    if (!def) return;
    const place = findView(this.#layout.value, ref);
    if (place?.window) this.#lastRects.set(ref, place.window.rect);
    this.#setLayout(this.#placeAtHome(removeView(this.#layout.value, ref), ref, def, true, into), ref);
  }

  /** A window's views into another group (a window docked whole). */
  merge(fromGroup: string, intoGroup: string): void {
    this.#setLayout(mergeGroup(this.#layout.value, fromGroup, intoGroup));
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
    // What the person hid stays hidden (a move shows a view: not this one).
    next = keepHidden(next, l.hidden);
    this.#setLayout(setGroupState(next, rail.id, rail.active, rail.collapsed === true));
  }

  // ---------------------------------------------------------------------------------------------------------------
  // The main area's region (shell §2.1)

  /**
   * `ref` beside a node of the main area's region (`of`: the app's content, a region group's id, or a view's ref: its
   * group), on `side`, `size` px (default `DefaultGroupSize`, at least `MinPaneSize`). A view opened anew comes in; a
   * placed one moves (from the right column, a window, another pane). Ends Expand and a maximize.
   */
  split(ref: ViewRef, side: StageSide, of = 'content', size?: number): void {
    const l = this.#layout.value;
    const at = of === 'content' ? of : (findGroup(l, of)?.group.id ?? findView(l, of)?.group.id);
    if (!at || (at !== 'content' && !inRegion(l, at))) return;
    batch(() => {
      if (this.#expanded.value) this.#expanded.set(false);
      this.#maximized.set(undefined);
      this.#setLayout(moveView(l, ref, { kind: 'edge', of: at, side, size }));
    });
  }

  /**
   * The main area's panes as a grid (shell §2.1, `arrangeGrid`): equal cells over the main area, lined up (a line across
   * moves as one) when the panes fill it. `columns` defaults to the square root of the panes, rounded up.
   */
  arrangeGrid(columns?: number): void {
    batch(() => {
      this.#maximized.set(undefined);
      this.#setLayout(arrangeGrid(this.#layout.value, this.#area.value, columns, this.#groupShows));
      // Asked before the main area was measured (an app arranging as it starts): arranged again at its first size.
      this.#gridPending = this.#area.value.width > 0 ? undefined : { columns };
    });
  }

  /** `ref` overlaid on the top of the main area's content (shell §2.1): the app places it (a stage's band slot). */
  toTop(ref: ViewRef): void {
    this.#setLayout(moveView(this.#layout.value, ref, { kind: 'top' }), ref);
  }

  /** The shell's `top` option: where a `top` view goes. */
  get top(): TopMode {
    return this.#top;
  }

  /**
   * The pane shown alone in the main area (shell §2.1; its id, `content` or a view's ref): a pane of the region (the
   * rest away and paused; a tab switch in it follows) or a window's view (out of its window over the main area, back
   * after); Escape back; undefined: every pane.
   */
  get maximized(): ReadonlyState<string | undefined> {
    return this.#maximized.asReadonly;
  }

  /**
   * Shows a pane alone (`content` or a view's ref), or every pane (undefined). A view in a region group that shows
   * another becomes the shown one; nothing for what can not show alone (a view docked beside the main area, on the
   * content's top, the content without a region).
   */
  maximize(pane: string | undefined): void {
    if (pane === undefined) return this.#maximized.set(undefined);
    const l = this.#layout.value;
    if (pane === 'content') return this.#maximized.set(mainNode(l).root ? pane : undefined);
    const place = findView(l, pane);
    if (!place) return;
    if (place.window) return this.#maximized.set(pane);
    if (!layoutIndex(l).regionNodes.has(place.group.id)) return;
    batch(() => {
      if (place.group.active !== pane) this.#setLayout(selectView(l, pane));
      this.#maximized.set(pane);
    });
  }

  toggleMaximize(pane: string): void {
    this.maximize(this.#maximized.value === pane ? undefined : pane);
  }

  /**
   * The pane of the main area's region a person or an agent last worked in (pressed, focused, chosen): `content` or a
   * view's ref. What follows it (a shared toolbar, the status) is the app's (shell §2.1).
   */
  get activePane(): ReadonlyState<string> {
    return this.#activePane.asReadonly;
  }

  setActivePane(pane: string): void {
    if (pane !== this.#activePane.value) this.#activePane.set(pane);
  }

  /**
   * A pane's status (`content` or a view's ref; shell §2.1): active, hovered, focused, maximized, several panes. One
   * state a pane, told only when a field changes. Hover and focus are the view's inputs (`setPaneHovered`,
   * `setPaneFocused`); headless, they stay off.
   */
  paneStatus(pane: string): ReadonlyState<PaneStatus> {
    let st = this.#paneStatus.get(pane);
    if (!st) {
      // Whether it is a pane (in the region, in a window): a layout write that moves no pane (a drag's sizes, a
      // window's move) changes no status.
      const isPane = derive(this.#layout, (l) => !!this.#paneNodeIn(l, pane) || !!findView(l, pane)?.window);
      st = derive(
        [this.#activePane, this.#hoveredPane, this.#focusedPane, this.#maximized, isPane, this.#several],
        (active, hovered, focused, maximized, is, several): PaneStatus => ({
          active: is && active === pane,
          hovered: is && hovered === pane,
          focused: is && focused === pane,
          maximized: is && maximized === pane,
          several,
        }),
        (a, b) =>
          a.active === b.active &&
          a.hovered === b.hovered &&
          a.focused === b.focused &&
          a.maximized === b.maximized &&
          a.several === b.several,
      );
      this.#paneStatus.set(pane, st);
    }
    return st;
  }

  /**
   * A pane's port (`content` or a view's ref; shell §2.1): one a pane, made once. A view gets its own as `ctx.pane`, the
   * app's content as `main`'s `pane`.
   */
  panePort(pane: string): PanePort {
    let port = this.#ports.get(pane);
    if (!port) {
      port = {
        id: pane,
        status: this.paneStatus(pane),
        insets: this.paneInsets(pane),
        setTitle: (title) => this.setPaneTitle(pane, title),
        toggleMaximized: () => this.toggleMaximize(pane),
        float: () => {
          if (pane !== 'content') this.float(pane);
        },
        close: () => {
          if (pane !== 'content') this.close(pane);
        },
      };
      this.#ports.set(pane, port);
    }
    return port;
  }

  /** A pane's own title (the pane's port: `setTitle`). The last set wins; the undo takes only its own back. */
  setPaneTitle(pane: string, title: string): () => void {
    const map = new Map(this.#titles.value);
    map.set(pane, title);
    this.#titles.set(map);
    return () => {
      if (this.#titles.value.get(pane) !== title) return;
      const next = new Map(this.#titles.value);
      next.delete(pane);
      this.#titles.set(next);
    };
  }

  /** A pane's title: its own (`setPaneTitle`), else its view's (the app's content: "Main"). */
  paneTitle(pane: string): ReadonlyState<string> {
    let st = this.#paneTitles.get(pane);
    if (!st) {
      st = derive([this.#titles, this.#views], (m) => m.get(pane) ?? this.view(pane)?.title ?? this.labels.mainTitle);
      this.#paneTitles.set(pane, st);
    }
    return st;
  }

  /**
   * What of a pane's bottom edge the shell's chrome over the main area covers (shell §2.1): the pane at the main area's
   * bottom left gets the `over` status's width, the one at the bottom right the `over` corner's; a pane shown alone
   * (maximized) both. None in a window.
   */
  paneInsets(pane: string): ReadonlyState<PaneInsets> {
    let st = this.#paneInsets.get(pane);
    if (!st) {
      st = derive(
        [this.#corners, this.#overChrome],
        (corners, chrome): PaneInsets =>
          !chrome.bottomLeft && !chrome.bottomRight
            ? NoInsets
            : {
                bottomLeft: corners.left === pane ? chrome.bottomLeft : 0,
                bottomRight: corners.right === pane ? chrome.bottomRight : 0,
              },
        (a, b) => a.bottomLeft === b.bottomLeft && a.bottomRight === b.bottomRight,
      );
      this.#paneInsets.set(pane, st);
    }
    return st;
  }

  /** The widths of the shell's chrome over the main area's bottom corners (`ShellView`'s `over`): the view's input. */
  setOverChrome(chrome: PaneInsets): void {
    const c = this.#overChrome.value;
    if (c.bottomLeft !== chrome.bottomLeft || c.bottomRight !== chrome.bottomRight) this.#overChrome.set(chrome);
  }

  /** Whether a group has a view that shows (registered, not flagged off, not hidden): a pane that shows. */
  readonly #groupShows = (g: GroupNode): boolean => {
    for (const v of g.views) if (this.shows(v)) return true;
    return false;
  };

  /** How many panes of the main area's region show (the content and the groups with a view that shows), up to `max`. */
  #regionPanes(root: LayoutNode | undefined, max: number): number {
    if (!root) return 1;
    let count = 0;
    const walk = (n: LayoutNode) => {
      if (count >= max) return;
      if (n.kind === 'content' || (n.kind === 'group' && this.#groupShows(n))) count++;
      else if (n.kind === 'split') for (const c of n.children) walk(c);
    };
    walk(root);
    return count;
  }

  /**
   * The panes in the main area's bottom corners, of those that show: down the last child of each column, the first or
   * last of each row.
   */
  #cornerPanes(root: LayoutNode | undefined): { left?: string; right?: string } {
    if (!root) return { left: 'content', right: 'content' };
    const shows = (n: LayoutNode): boolean =>
      n.kind === 'content' || (n.kind === 'group' ? this.#groupShows(n) : n.kind === 'split' && n.children.some(shows));
    const walk = (n: LayoutNode, side: 'left' | 'right'): string | undefined => {
      if (n.kind === 'content') return 'content';
      if (n.kind === 'group') return n.active ?? n.views[0];
      if (n.kind !== 'split') return undefined;
      const fromEnd = n.dir === 'column' || side === 'right';
      let c: LayoutNode | undefined;
      if (!fromEnd) c = n.children.find(shows);
      else for (let i = n.children.length - 1; i >= 0 && !c; i--) if (shows(n.children[i])) c = n.children[i];
      return c && walk(c, side);
    };
    return { left: walk(root, 'left'), right: walk(root, 'right') };
  }

  /** The pane the pointer is over (`content`, a view's ref), or none: the view's input. */
  setPaneHovered(pane: string | undefined): void {
    if (pane !== this.#hoveredPane.value) this.#hoveredPane.set(pane);
  }

  /** The pane the focus is in, or none: the view's input. */
  setPaneFocused(pane: string | undefined): void {
    if (pane !== this.#focusedPane.value) this.#focusedPane.set(pane);
  }

  /** The region node a pane is in (`content`, or the group of a view's ref), if it is in the region. */
  paneNode(pane: string): string | undefined {
    return this.#paneNodeIn(this.#layout.value, pane);
  }

  #paneNodeIn(l: ShellLayout, pane: string): string | undefined {
    if (pane === 'content') return pane;
    const place = findView(l, pane);
    return place && !place.window && inRegion(l, place.group.id) ? place.group.id : undefined;
  }

  /** A docked view into a new group under its own (Split below). */
  splitBelow(ref: ViewRef): void {
    const place = findView(this.#layout.value, ref);
    if (!place || place.window || place.group.views.length < 2) return;
    this.move(ref, { kind: 'edge', of: place.group.id, side: 'bottom', size: DefaultSplitSize });
  }

  /** A docked view into the group above or below its own in its column. */
  moveInColumn(ref: ViewRef, dir: -1 | 1): void {
    this.#setLayout(moveInColumn(this.#layout.value, ref, dir));
  }

  /**
   * A docked node's size: at least `MinGroupSize`, at most what leaves the main area `MinMainSize` (the shell's width
   * known); a pane of the main area's region at least `MinPaneSize`, at most what leaves the content its least size
   * (`#clampRegion`, the main area's size known).
   */
  resize(nodeId: string, size: number): void {
    const l = resizeNode(this.#layout.value, nodeId, size);
    this.#setLayout(this.#clampSizes(regionHas(l, nodeId) ? this.#clampRegion(l, nodeId) : l));
  }

  /**
   * The region's sizes kept within the main area (shell §2.1): in each split, its sized children leave the child that
   * takes the rest its least size (its even share of the room, at most `MinMainSize`); over that, `prefer` (the node
   * being resized, with its grid peers) gives way first, then the rest shrink in proportion, each at least
   * `MinPaneSize`. Nothing while the main area's size is not known. Why: a main area made smaller (a window, a panel
   * opened beside it) must not squeeze the content to nothing.
   */
  #clampRegion(l: ShellLayout, prefer?: string): ShellLayout {
    const area = this.#area.value;
    const root = mainNode(l).root;
    if (!area.width || !area.height || root?.kind !== 'split') return l;
    const preferred = prefer ? new Set([prefer, ...gridPeers(l, prefer)]) : undefined;
    const sizes = new Map<string, number>();
    const walk = (n: LayoutNode, width: number, height: number) => {
      if (n.kind !== 'split') return;
      const row = n.dir === 'row';
      const room = row ? width : height;
      const now = new Map<string, number>();
      let total = 0;
      for (const c of n.children) {
        if ((c.kind === 'group' || c.kind === 'split') && c.size !== undefined) {
          now.set(c.id, c.size);
          total += c.size;
        }
      }
      const flex = n.children.length - now.size;
      const least = Math.min(MinMainSize, Math.max(MinPaneSize, Math.round(room / n.children.length)));
      let over = total - (room - flex * least);
      if (over > 0) {
        for (const [id, s] of now) {
          if (!preferred?.has(id) || over <= 0) continue;
          const to = Math.max(MinPaneSize, s - over);
          over -= s - to;
          now.set(id, to);
        }
        if (over > 0) {
          let rest = 0;
          for (const [id, s] of now) if (!preferred?.has(id)) rest += s;
          const k = rest > 0 ? Math.max(0, rest - over) / rest : 0;
          for (const [id, s] of now) if (!preferred?.has(id)) now.set(id, Math.max(MinPaneSize, Math.floor(s * k)));
        }
        for (const c of n.children) {
          const s = now.get(c.id);
          if (s !== undefined && s !== (c as GroupNode | SplitNode).size) sizes.set(c.id, s);
        }
      }
      let used = 0;
      for (const s of now.values()) used += s;
      for (const c of n.children) {
        const s = now.get(c.id) ?? Math.max(0, room - used);
        walk(c, row ? s : width, row ? height : s);
      }
    };
    walk(root, area.width, area.height);
    return sizes.size ? setSizes(l, sizes) : l;
  }

  /**
   * The docked sizes of the root row kept within the shell's width: their sum, the rail and the main area's least size
   * (`MinMainSize`) fit; the largest gives way first. Nothing while the width is not known.
   */
  #clampSizes(l: ShellLayout): ShellLayout {
    const width = this.#width.value;
    if (!width || l.root.kind !== 'split' || l.root.dir !== 'row') return l;
    const room = width - (l.rail ? RailWidth : 0) - MinMainSize;
    let next = l;
    for (let guard = 0; guard < 8; guard++) {
      const root = next.root as SplitNode;
      const sized = root.children.filter(
        (c): c is GroupNode | SplitNode => (c.kind === 'group' || c.kind === 'split') && c.size !== undefined,
      );
      const total = sized.reduce((t, c) => t + (c.size ?? 0), 0);
      if (total <= room || !sized.length) break;
      const largest = sized.reduce((a, b) => ((b.size ?? 0) > (a.size ?? 0) ? b : a));
      const over = total - room;
      const to = Math.max(MinGroupSize, (largest.size ?? 0) - over);
      if (to === largest.size) break;
      next = resizeNode(next, largest.id, to);
    }
    return next;
  }

  /** A window's rectangle: at least its views' least size, inside the main area (`#fitWindows`). */
  setWindowRect(groupId: string, rect: Rect): void {
    this.#setLayout(setWindowRect(this.#layout.value, groupId, rect));
  }

  raise(groupId: string): void {
    this.#setLayout(raiseWindow(this.#layout.value, groupId));
  }

  /** The views out of view by the person's choice (`close` of a docked single view; `show` brings one back). */
  get hiddenViews(): ReadonlyState<readonly ViewRef[]> {
    return this.#hiddenViews;
  }

  /** Shows or hides the rail (and the docked column with it). */
  toggleRail(shown = !this.#layout.value.rail): void {
    this.#setLayout(setRail(this.#layout.value, shown));
  }

  get railShown(): ReadonlyState<boolean> {
    return this.#railShown;
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

  // ---------------------------------------------------------------------------------------------------------------
  // Expand and Fullscreen (shell §3.9)

  /**
   * Expanded: the shell's own chrome steps away (the rail, the docked column, the windows: hidden, paused) and the main
   * area takes the shell, showing what it shows (every pane, or the maximized one); the layers (card, page, modal) still
   * open over it. The shell's presentation, not a pane's: a pane alone is `maximize` (shell §2.1); an app wanting both
   * from one button calls both. What is open, not the layout: not saved. Showing a docked view or a window leaves it.
   */
  get expanded(): ReadonlyState<boolean> {
    return this.#expanded.asReadonly;
  }

  /** Expands (`true`), steps back (`false`), or toggles (no argument). */
  toggleExpanded(on = !this.#expanded.value): void {
    this.#expanded.set(on);
  }

  /** Whether the shell is fullscreen now (the view writes back what the browser did, its own Esc included). */
  get fullscreen(): ReadonlyState<boolean> {
    return this.#fullscreen.asReadonly;
  }

  /** What was asked of fullscreen (the view does it, then writes the result); undefined: nothing pending. */
  get fullscreenRequest(): ReadonlyState<boolean | undefined> {
    return this.#fullscreenRequest.asReadonly;
  }

  /**
   * Whether fullscreen can be asked for here, and if not, why (an iframe that forbids it; headless: there is no screen).
   * The view says so as it mounts.
   */
  get fullscreenAvailable(): ReadonlyState<true | string> {
    return this.#fullscreenAvailable.asReadonly;
  }

  /** Asks for fullscreen (or out of it); nothing where it is not available. */
  requestFullscreen(on = !this.#fullscreen.value): void {
    if (this.#fullscreenAvailable.value !== true) return;
    this.#fullscreenRequest.set(on);
  }

  /** The view's: what the browser did (and the request is done). */
  setFullscreen(on: boolean): void {
    batch(() => {
      this.#fullscreen.set(on);
      this.#fullscreenRequest.set(undefined);
    });
  }

  /** The view's: whether fullscreen can be asked for, or why not. */
  setFullscreenAvailable(available: true | string): void {
    this.#fullscreenAvailable.set(available);
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
    if (width === this.#width.value) return;
    batch(() => {
      this.#width.set(width);
      // A narrower shell: the docked sizes give way so the main area keeps its least width.
      this.#setLayout(this.#clampSizes(this.#layout.value));
    });
  }

  /** The main area's size: windows open and stay inside it (`#fitWindows`), the region's panes give way. */
  setArea(width: number, height: number): void {
    const a = this.#area.value;
    if (a.width === width && a.height === height) return;
    batch(() => {
      this.#area.set({ width, height });
      let l = this.#layout.value;
      if (this.#gridPending && width > 0) {
        l = arrangeGrid(l, { width, height }, this.#gridPending.columns, this.#groupShows);
        this.#gridPending = undefined;
      }
      // The region's panes give way so the content keeps its least size (shell §2.1).
      this.#setLayout(this.#clampRegion(l));
    });
  }

  get area(): ReadonlyState<{ width: number; height: number }> {
    return this.#area.asReadonly;
  }

  get mode(): ReadonlyState<ShellMode> {
    return this.#mode;
  }

  /**
   * A new window's rectangle (shell §2, window geometry): its last this session, else beside `anchorRect`, else
   * cascaded; embedded at 0.8 of its size; at least its `minSize` (`#fitWindows` keeps it so).
   */
  #newWindowRect(ref: ViewRef, def: ShellViewDef, anchorRect?: Rect): Rect {
    const area = this.#area.value.width ? this.#area.value : UnmeasuredArea;
    const last = def.remember === false ? undefined : this.#lastRects.get(ref);
    if (last) return clampRect(last, area);
    const size = def.size ?? DefaultWindowSize;
    const k = this.#mode.value === 'embedded' ? 0.8 : 1;
    const min = def.minSize;
    const scaled = {
      width: Math.max(Math.round(size.width * k), min?.width ?? 0),
      height: Math.max(Math.round(size.height * k), min?.height ?? 0),
    };
    return placeWindow(
      scaled,
      area,
      stacked(this.#layout.value).map((w) => w.rect),
      anchorRect,
    );
  }
}
