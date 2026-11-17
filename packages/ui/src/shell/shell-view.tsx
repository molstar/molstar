/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { derive, listen, type ReadonlyState, State } from '@molstar/core/reactive';
import { Button } from '../components/button/button.js';
import { type ContextMenuItem, openContextMenu, openMenuBeside } from '../components/context-menu/context-menu.js';
import { Icon, type IconName } from '../components/icon/icon.js';
import { installTooltips } from '../components/tooltip/tooltip.js';
import { dismissable } from '../components/popover/dismiss.js';
import {
  activity,
  type Child,
  deepActiveElement,
  drag,
  effect,
  eventTarget,
  focusables,
  focusedByKeyboard,
  For,
  frameState,
  Motion,
  observeResize,
  on,
  onActive,
  onCleanup,
  rememberFocus,
  requireScope,
  rovingFocus,
  runInScope,
  Show,
  setAttr,
  trapTab,
  type ViewScope,
  watch,
} from '../dom/index.js';
import { classes } from './_generated/shell.module.css.js';
import {
  TopGroup,
  DefaultColumnSize,
  DefaultGroupSize,
  DefaultWindowSize,
  findGroup,
  MinGroupSize,
  MinPaneSize,
  findView,
  gridPeers,
  inRegion,
  type GroupNode,
  type LayoutIndex,
  type LayoutNode,
  layoutIndex,
  mainNode,
  nodeOf,
  type Rect,
  railGroup,
  rightColumn,
  type ShellLayout,
  type SplitNode,
  stacked,
  type ViewRef,
} from './layout.js';
import { fitRail, type RailEntry, railEntries } from './rail.js';
import {
  type PanePort,
  type RailItem,
  type ShellModel,
  type ShellViewDef,
  keyOf,
  type ViewPlacement,
  viewIdOf,
} from './shell-model.js';
import { viewBody, viewContext, ViewFrame } from './view-frame.js';
import { paneMenu, railMenu, viewMenu, type WindowAdjust } from './view-menu.js';

/**
 * What holds the app (shell §3.6, busy): a label, a progress 0–1 (none: indeterminate), a cancel. The overlay shows the
 * current value: a task replacing another changes its label, its progress and its Cancel (shown while it has one).
 */
export interface ShellBusy {
  readonly label: string;
  readonly progress?: number;
  cancel?(): void;
}

export interface ShellViewProps {
  readonly shell: ShellModel;
  /**
   * The app's main view (a canvas, a document): rendered once, kept across modes. `top` is the views overlaid on its top
   * (shell §2.1), for the app to place (a stage's band slot: `<Stage top={top()} …/>`).
   */
  readonly main: (ctx: {
    top: () => Child;
    /** The app's content as a pane (shell §2.1): its status, its insets, its title. */
    pane: PanePort;
  }) => Child;
  /**
   * Chrome over the whole main area, once (shell §2.1): its `status` at the bottom left, its `corner` at the bottom right,
   * over every pane and under the windows (a toolbar shared by the panes follows `shell.activePane`). Each stage keeps
   * its own slots too: which chrome goes where is the app's.
   */
  readonly over?: { readonly status?: () => Child; readonly corner?: () => Child };
  /** While set, a scrim and its card hold the app. */
  readonly busy?: ReadonlyState<ShellBusy | undefined>;
  /** Files dropped anywhere on the shell; without it, no drop overlay (it says `labels.dropFiles`). */
  readonly onDropFiles?: (files: File[]) => void;
  /** What can be dropped, a hint under the overlay's words (a list of formats). */
  readonly dropFormats?: string;
  /**
   * How far Expand reaches (shell §3.9): `shell` (default: the main area takes the shell's box) or `page` (the shell
   * itself fills the page while expanded: a viewer embedded in a page).
   */
  readonly fill?: 'shell' | 'page';
  /**
   * The rail can be hidden (default on): a Hide button at its foot, and a button in the main area's corner that shows
   * it again, with `railToggleIcon` (the app's mark; default the Mol* star).
   */
  readonly railToggle?: boolean;
  readonly railToggleIcon?: IconName;
  /**
   * The shell's own mark of the active pane while several show (an accent line along its top edge; default on). Off:
   * the panes mark themselves from their status (`ctx.pane`), or nothing does.
   */
  readonly paneMark?: boolean;
  /**
   * A ⋯ button in each view's header for its layout menu (Float, Dock, Split below, …). Off by default: the same menu
   * is a right-click (or Shift+F10) on a view's header, its tab or its rail button; ✕ closes.
   */
  readonly viewMenu?: boolean;
  /**
   * The shell's tooltips (default on): one layer for everything in it, the app's views included (`data-tip`, and a
   * row's card beside its view, `data-tip-row`). Off when the app installs its own around the shell.
   */
  readonly tooltips?: boolean;
  readonly class?: string;
}

/**
 * The shell (shell §1–4): the app's main area (its region of panes) with the docked layout around it (the right column
 * of groups, the rail at the edge), windows over the main area, a card in the right column's place, a page over the
 * main area, a modal over everything; compact below its compact width (the views in a bottom sheet). Mode by its own
 * width, the rail's steps by its height. Its accessible name is `labels.workspace`.
 */
export function ShellView(props: ShellViewProps) {
  const { shell } = props;
  const t = shell.labels;
  const owner = requireScope('ShellView');
  const mainArea = (
    <div
      class={classes.main}
      data-msui-shell-main=""
      data-msui-region=""
      tabIndex={-1}
      aria-label={t.mainArea}
      role="region"
    />
  ) as HTMLElement;
  const bodies = new ViewBodies(shell, owner, {
    menuButton: !!props.viewMenu,
    railToggle: props.railToggle !== false,
    main: mainArea,
  });
  onCleanup(() => bodies.dispose());
  // The app's main view: made once, moved between the desktop and compact layouts (its state, a canvas, kept); its
  // top overlay is the app's to place (shell §2.1).
  const mainContent = (
    <div class={classes.mainContent}>
      {props.main({
        top: () => <TopOverlay shell={shell} bodies={bodies} />,
        pane: shell.panePort('content'),
      })}
    </div>
  ) as HTMLElement;
  // The main area's region: the content and the views docked beside it (shell §2.1); the app's chrome over it.
  mainArea.append(Region({ shell, bodies, content: mainContent, mark: props.paneMark !== false }));
  if (props.over) {
    const status = (<div class={classes.overStatus}>{props.over.status?.() ?? null}</div>) as HTMLElement;
    const corner = (<div class={classes.overCorner}>{props.over.corner?.() ?? null}</div>) as HTMLElement;
    mainArea.append(
      (
        <div class={classes.over}>
          {status}
          {corner}
        </div>
      ) as HTMLElement,
    );
    // What the shared chrome covers of the panes in the main area's bottom corners (shell §2.1): its widths, and the
    // gap a pane keeps from it (the bar's padding, a unit each side), told to the panes through their ports.
    const covered = (el: HTMLElement) => (el.offsetWidth > 0 ? el.offsetWidth + 16 : 0);
    const measure = () => shell.setOverChrome({ bottomLeft: covered(status), bottomRight: covered(corner) });
    observeResize(status, measure);
    observeResize(corner, measure);
  }
  observeResize(mainArea, () => shell.setArea(mainArea.clientWidth, mainArea.clientHeight));
  // Under a page, the app's main view takes no pointer or focus (the page covers it).
  watch([shell.page], () => setAttr(mainContent, 'inert', !!shell.page.value));

  const root = (
    <div
      class={props.class ? `${classes.shell} ${props.class}` : classes.shell}
      role="application"
      aria-roledescription={t.workspaceRole}
      aria-label={t.workspace}
      data-msui-shell=""
    />
  ) as HTMLElement;
  observeResize(root, () => shell.setWidth(root.clientWidth));
  if (props.tooltips !== false) onCleanup(installTooltips(root));
  watch([shell.mode], () => setAttr(root, 'data-mode', shell.mode.value));
  // Expand (shell §3.9): the main area alone; with `fill: 'page'`, the shell over the page.
  watch([shell.expanded], () => {
    setAttr(root, 'data-expanded', shell.expanded.value);
    setAttr(root, 'data-fill', shell.expanded.value && props.fill === 'page' ? 'page' : false);
  });
  // Escape steps back, after whatever is inside took it (a menu, a card: they prevent its default): out of Expand, then
  // out of a maximized pane.
  on<KeyboardEvent>(root, 'keydown', (e) => {
    if (e.key !== 'Escape' || e.defaultPrevented) return;
    if (shell.expanded.value) shell.toggleExpanded(false);
    else if (shell.maximized.value) shell.maximize(undefined);
    else return;
    e.preventDefault();
  });
  // Fullscreen: the model asks, the view does it on the shell's root and writes back what the browser did (its own Esc
  // too); where the page or an iframe forbids it, the model says why.
  const doc = root.ownerDocument;
  shell.setFullscreenAvailable(
    typeof root.requestFullscreen !== 'function'
      ? t.fullscreenNoApi
      : doc.fullscreenEnabled === false
        ? t.fullscreenNotAllowed
        : true,
  );
  effect(shell.fullscreenRequest, (on) => {
    if (on === undefined) return;
    const done = () => shell.setFullscreen(doc.fullscreenElement === root);
    const ask = on ? root.requestFullscreen?.() : doc.fullscreenElement ? doc.exitFullscreen?.() : undefined;
    if (ask && typeof ask.then === 'function') ask.then(done, done);
    else done();
  });
  on(doc, 'fullscreenchange', () => shell.setFullscreen(doc.fullscreenElement === root));
  on<KeyboardEvent>(root, 'keydown', (e) => {
    if (e.key === 'F6' && !e.defaultPrevented) {
      e.preventDefault();
      cycleRegions(root, e.shiftKey ? -1 : 1);
    }
  });

  const compact = derive(shell.mode, (m) => m === 'compact');
  // Everything but the layers that hold the app (a modal, the busy overlay): inert while one is up (the rest of the
  // shell takes no pointer, focus or assistive tech), decided as it changes, whatever is in it then.
  const layout = (
    <div class={classes.layout}>
      <Show when={compact} fallback={() => <Desktop shell={shell} bodies={bodies} main={mainArea} props={props} />}>
        {() => <Compact shell={shell} bodies={bodies} main={mainArea} />}
      </Show>
    </div>
  ) as HTMLElement;
  const holding = derive([shell.modal, props.busy ?? new State(undefined)], (m, b) => !!m || !!b);
  watch([holding], () => setAttr(layout, 'inert', holding.value));
  root.append(
    (
      <>
        {layout}
        <ModalHost shell={shell} bodies={bodies} />
        {props.busy ? <BusyOverlay busy={props.busy} cancel={t.cancel} /> : null}
        {props.onDropFiles ? (
          <FileDrop onDrop={props.onDropFiles} label={t.dropFiles} formats={props.dropFormats} />
        ) : null}
      </>
    ) as Node,
  );
  // The layers over the main area (in it in both modes): windows (desktop only), a window's maximized pane, the page,
  // the rail toggle.
  mainArea.append(
    (
      <>
        {/* Covered by a page, the windows are inert and paused (hidden UI does no work, shell §1). */}
        <Show when={derive(compact, (c) => !c)}>
          {() =>
            activity(
              // Covered by a page: inert and paused (away while expanded: each window pauses its own frame, below).
              derive(shell.page, (p) => !p),
              () => <WindowsLayer shell={shell} bodies={bodies} />,
            )
          }
        </Show>
        <MaximizedHost shell={shell} bodies={bodies} />
        <PageHost shell={shell} bodies={bodies} />
        {props.railToggle === false ? null : (
          <RailToggle shell={shell} compact={compact} icon={props.railToggleIcon ?? 'molstar'} />
        )}
      </>
    ) as Node,
  );
  return root;
}

// ---------------------------------------------------------------------------------------------------------------------
// View bodies: kept per view while it moves

/**
 * Each placed view's body (shell §3.1): made once, in a scope of its own, and handed to whichever frame shows the view
 * now (a group's, a window's, the sheet's); paused while no frame shows it, so a view keeps its state (a scroll, a
 * draft, a plot's zoom) when it is popped out, docked back or folded. Disposed when the view leaves the layout.
 */
class ViewBodies {
  readonly #bodies = new Map<
    ViewRef,
    { el: HTMLElement; scope: ViewScope; owner: object | undefined; def: ShellViewDef }
  >();
  readonly #placements = new Map<ViewRef, ViewPlacement>();

  constructor(
    readonly shell: ShellModel,
    readonly owner: ViewScope,
    /** What the shell's frames show (its props): the ⋯ button, the rail's Hide; the shell's main area. */
    readonly ui: { readonly menuButton: boolean; readonly railToggle: boolean; readonly main: HTMLElement },
  ) {
    // A view that leaves the layout and every layer, or is unregistered: its body goes. Followed by what is placed,
    // registered and layered (compared, not keyed by a string: not every layout write, a drag's sizes, sweeps).
    effect(
      derive(
        [shell.layout, shell.views, shell.card, shell.page, shell.modal],
        (l, v, c, p, m) => ({ l, v, c: c?.ref, p: p?.ref, m: m?.ref }),
        (a, b) => a.v === b.v && a.c === b.c && a.p === b.p && a.m === b.m && samePlaced(a.l, b.l),
      ),
      () => this.#sweep(),
    );
  }

  /** The body of `ref`, now shown by `owner` (a frame): active, `placement` where it shows. */
  take(ref: ViewRef, owner: object, placement: ViewPlacement): HTMLElement | undefined {
    const def = this.shell.view(ref);
    if (!def) return undefined;
    this.#placements.set(ref, placement);
    let b = this.#bodies.get(ref);
    // Registered anew (an extension reloaded): the old body was rendered by the old definition.
    if (b && b.def !== def) {
      this.#drop(ref);
      b = undefined;
    }
    if (!b) {
      const scope = this.owner.child();
      const ctx = viewContext(this.shell, ref, () => this.#placements.get(ref) ?? placement);
      const el = runInScope(scope, () => viewBody(def, ctx));
      b = { el, scope, owner: undefined, def };
      this.#bodies.set(ref, b);
    }
    b.owner = owner;
    b.scope.setActive(true);
    return b.el;
  }

  /** The frame `owner` showing it is paused or resumed (a folded group, a hidden column): its body with it. */
  setActive(ref: ViewRef, owner: object, active: boolean): void {
    const b = this.#bodies.get(ref);
    if (b && b.owner === owner) b.scope.setActive(active);
  }

  /** `owner` no longer shows it: paused (unless another frame took it meanwhile), disposed when it left for good. */
  release(ref: ViewRef, owner: object): void {
    const b = this.#bodies.get(ref);
    if (!b || b.owner !== owner) return;
    b.owner = undefined;
    b.scope.setActive(false);
    b.el.remove();
    // The sweep that followed the change ran while this frame still held it.
    if (this.#gone(ref)) this.#drop(ref);
  }

  /** Whether `ref` left the layout and every layer, or is no longer registered (or registered anew). */
  #gone(ref: ViewRef): boolean {
    const b = this.#bodies.get(ref);
    if (b && this.shell.view(ref) !== b.def) return true;
    const layered = [this.shell.card.value, this.shell.page.value, this.shell.modal.value].some((x) => x?.ref === ref);
    return !layered && !findView(this.shell.layout.value, ref);
  }

  #drop(ref: ViewRef): void {
    const b = this.#bodies.get(ref);
    if (!b) return;
    b.scope.dispose();
    b.el.remove();
    this.#bodies.delete(ref);
    this.#placements.delete(ref);
  }

  #sweep(): void {
    for (const [ref, b] of [...this.#bodies]) {
      if (!b.owner && this.#gone(ref)) this.#drop(ref);
    }
  }

  dispose(): void {
    for (const b of this.#bodies.values()) b.scope.dispose();
    this.#bodies.clear();
  }
}

/** Whether two lists of rail entries are the same buttons (ids, and the def or item each shows). */
function sameEntries(a: readonly RailEntry[], b: readonly RailEntry[]): boolean {
  return (
    a.length === b.length &&
    a.every((e, i) => e.id === b[i].id && e.view?.def === b[i].view?.def && e.item === b[i].item)
  );
}

/** Whether two layouts place the same views (wherever). */
function samePlaced(a: ShellLayout, b: ShellLayout): boolean {
  if (a === b) return true;
  const x = layoutIndex(a).views;
  const y = layoutIndex(b).views;
  if (x.size !== y.size) return false;
  for (const ref of x.keys()) if (!y.has(ref)) return false;
  return true;
}

function sameList(a: readonly string[], b: readonly string[]): boolean {
  if (a === b) return true;
  if (a.length !== b.length) return false;
  for (let i = 0; i < a.length; i++) if (a[i] !== b[i]) return false;
  return true;
}

/** Whether two values of a group read the same (its views, the one shown, its fold, its size). */
function sameGroup(a: GroupNode | undefined, b: GroupNode | undefined): boolean {
  return (
    a === b ||
    (!!a &&
      !!b &&
      a.id === b.id &&
      a.active === b.active &&
      !!a.collapsed === !!b.collapsed &&
      a.size === b.size &&
      sameList(a.views, b.views))
  );
}

/**
 * Where the focus goes when a layer opens (shell §4, WAI-ARIA dialog): the first control of its content, never its
 * header's (its ✕ would be pressed by the next Enter; a modal's content puts its least destructive action first). With
 * `content: false` (a card opened by a press, not from the keyboard) or no control in the content, the layer itself
 * (no ring: nothing was chosen).
 */
function focusInitial(layer: HTMLElement, content: boolean): void {
  const body = layer.querySelector('[data-msui-view] > :not(header)');
  const first = content && body ? focusables(body)[0] : undefined;
  if (first) {
    first.focus({ preventScroll: true });
    return;
  }
  if (!layer.hasAttribute('tabindex')) layer.tabIndex = -1;
  layer.focus({ preventScroll: true });
}

/**
 * Marks what a region shows anew (not its first render) to fade in (styles §7.1: a swap crossfades in place): a view
 * replacing another in a group, a window, the sheet.
 */
function entering(): <T>(node: T) => T {
  let first = true;
  return (node) => {
    if (!first && node instanceof Element) node.setAttribute('data-reveal', 'fade');
    first = false;
    return node;
  };
}

/** A frame showing `ref`'s kept body (the body taken while it renders, given back when it goes). */
function KeptFrame(props: {
  shell: ShellModel;
  bodies: ViewBodies;
  ref: ViewRef;
  placement: ViewPlacement;
  head?: 'rail' | 'header' | 'none';
  tabs?: Child;
  menu?: () => ReturnType<typeof viewMenu>;
  onClose?: () => void;
  headerRef?: (el: HTMLElement) => void;
  /** While false, the body is hidden and paused, the header stays (a folded tab group). */
  open?: ReadonlyState<boolean>;
}) {
  const def = props.shell.view(props.ref);
  if (!def) return null;
  const token = {};
  const body = props.bodies.take(props.ref, token, props.placement);
  // The body pauses with the frame (the column hidden with the rail) and while folded, not only when it leaves.
  let frameActive = true;
  const sync = () => props.bodies.setActive(props.ref, token, frameActive && (props.open?.value ?? true));
  onActive((active) => {
    frameActive = active;
    sync();
  });
  if (props.open) watch([props.open], sync);
  onCleanup(() => props.bodies.release(props.ref, token));
  return (
    <ViewFrame
      shell={props.shell}
      ref={props.ref}
      def={def}
      placement={props.placement}
      head={props.head}
      tabs={props.tabs}
      menu={props.menu}
      onClose={props.onClose}
      menuButton={props.bodies.ui.menuButton}
      headerRef={props.headerRef}
      body={body}
    />
  );
}

// ---------------------------------------------------------------------------------------------------------------------
// Desktop: the docked tree and the rail

/** The views of a group that show (registered, not hidden, their flag not false). */
function shownViews(shell: ShellModel, g: GroupNode | undefined): ViewRef[] {
  return g ? g.views.filter((v) => shell.shows(v)) : [];
}

/** Whether a group has a view that shows (no list made: asked per node while the tree's state is written). */
function hasShown(shell: ShellModel, g: GroupNode): boolean {
  for (const v of g.views) if (shell.shows(v)) return true;
  return false;
}

/** The view a group shows: its active one if it shows, else its first that does. */
function activeOf(shell: ShellModel, g: GroupNode | undefined): ViewRef | undefined {
  const shown = shownViews(shell, g);
  return g?.active && shown.includes(g.active) ? g.active : shown[0];
}

function Desktop(props: { shell: ShellModel; bodies: ViewBodies; main: HTMLElement; props: ShellViewProps }) {
  const { shell, bodies } = props;
  const dock = (<div class={classes.dock} />) as HTMLElement;
  // The card's place in the dock, kept apart from the tree (which is replaced at each layout change).
  const cardSlot = (<div class={classes.cardSlot} />) as HTMLElement;
  const tree = DockTree({
    shell,
    bodies,
    leaf: 'main',
    leafEl: props.main,
    rootOf: (l) => l.root,
    // Covered by the card (it takes the right column's place): paused and inert, as hidden UI is.
    covered: (l, n) => !!shell.card.value && n.kind === 'group' && l.rail && layoutIndex(l).rightGroups.has(n.id),
    // Expanded (shell §3.9): the main area alone; every group away, paused.
    away: (_, n) => shell.expanded.value && n.kind === 'group',
    sources: [shell.card, shell.expanded],
  });
  dock.append(tree, cardSlot);
  // The card: over the right column (the rail beside it), or at the main area's right edge.
  const cardWidth = derive(shell.layout, (l) => {
    const col = rightColumn(l);
    return col && l.rail && (col.kind === 'group' || col.kind === 'split')
      ? (col.size ?? DefaultColumnSize)
      : DefaultColumnSize;
  });
  watch([cardWidth], () => dock.style.setProperty('--msui-shell-card', `${cardWidth.value}px`));
  dock.style.position = 'relative';
  return (
    <>
      {dock}
      <CardHost shell={shell} bodies={bodies} into={cardSlot} />
      <Show
        when={derive(
          [shell.layout, shell.absent, shell.views, shell.expanded],
          (l, _a, _v, e) => !!railGroup(l) && l.rail && !e,
        )}
        reveal="side"
      >
        {() => <Rail shell={shell} hideable={bodies.ui.railToggle} />}
      </Show>
    </>
  );
}

/**
 * The main area's region (shell §2.1): the app's content and the views docked beside it, one tree of the same groups,
 * splits and splitters as the shell's. Its panes: the pane a person works in is the shell's active pane (a press, the
 * focus), marked while there are several; a maximized pane shows alone, the rest away and paused, Escape back.
 */
function Region(props: { shell: ShellModel; bodies: ViewBodies; content: HTMLElement; mark: boolean }) {
  const { shell, bodies } = props;
  // The region's node shown alone: the maximized pane's, when it is the region's (a window's comes out over it
  // instead).
  const alone = (): string | undefined => {
    const m = shell.maximized.value;
    return m && shell.paneNode(m);
  };
  // Several panes shown (the model's: of the panes that show; a maximized one shows alone, no mark of being active).
  const several = derive(shell.paneStatus('content'), (s) => s.several);
  const marked = derive([shell.activePane, shell.layout], (p) => shell.paneNode(p));
  // The shell's own mark of the active pane (several panes); what a pane shows on its status is the pane's (`ctx.pane`).
  // Written as either changes and after each rebuild (the panes' wrappers are new then).
  const mark = (host: HTMLElement) => {
    for (const el of host.querySelectorAll('[data-pane]')) {
      setAttr(el, 'data-active', props.mark && several.value && el.getAttribute('data-pane') === marked.value);
    }
  };
  const tree = DockTree({
    shell,
    bodies,
    leaf: 'content',
    leafEl: props.content,
    rootOf: (l) => mainNode(l).root ?? { kind: 'content', id: 'content' },
    // Maximized (shell §2.1): only that pane and what holds it.
    away: (l, n) => {
      const m = alone();
      return !!m && !!mainNode(l).root && !holds(layoutIndex(l), n, m);
    },
    sources: [shell.maximized],
    motion: false,
    onBuild: mark,
  });
  tree.classList.add(classes.region);
  watch([several], () => setAttr(tree, 'data-panes', several.value));
  // The pane under a press or the focus is the active one.
  const paneOf = (target: EventTarget | null): string | undefined => {
    // A splitter is the layout's, not the pane's it sits in: pressing, focusing or hovering it changes no pane.
    if ((target as Element | null)?.closest?.('[role="separator"]')) return undefined;
    const node = (target as Element | null)?.closest?.('[data-pane]');
    if (!node || !tree.contains(node)) return undefined;
    const id = node.getAttribute('data-pane')!;
    if (id === 'content') return id;
    return activeOf(shell, findGroup(shell.layout.value, id)?.group);
  };
  const activate = (e: Event) => {
    const pane = paneOf(e.target);
    if (pane) shell.setActivePane(pane);
  };
  on(tree, 'pointerdown', activate, { capture: true });
  // The pane's status (shell §2.1): the view feeds the model where the pointer and the focus are.
  on(tree, 'focusin', (e: Event) => {
    activate(e);
    shell.setPaneFocused(paneOf(e.target));
  });
  on(tree, 'focusout', (e: FocusEvent) => {
    if (!paneOf(e.relatedTarget)) shell.setPaneFocused(undefined);
  });
  on(tree, 'pointerover', (e: Event) => shell.setPaneHovered(paneOf(e.target)));
  on(tree, 'pointerleave', () => shell.setPaneHovered(undefined));
  watch([marked, several], () => mark(tree));
  return tree;
}

/** Whether a node is a pane (`content` or a group's id) or holds it (up the pane's parents: no walk of the subtree). */
function holds(ix: LayoutIndex, n: LayoutNode, pane: string): boolean {
  for (let id: string | undefined = pane; id !== undefined; id = ix.parents.get(id)?.id) {
    if (id === n.id) return true;
  }
  return false;
}

interface DockTreeOptions {
  shell: ShellModel;
  bodies: ViewBodies;
  /** The level's leaf (shell §2.1): the main area at the shell's level, the app's content in the main area's region. */
  leaf: 'main' | 'content';
  leafEl: HTMLElement;
  rootOf: (l: ShellLayout) => LayoutNode;
  /** A node covered (inert and paused, still shown): the right column under a card. */
  covered?: (l: ShellLayout, n: LayoutNode) => boolean;
  /** A node away (not shown, paused): the groups while expanded, the other panes while one is maximized. */
  away?: (l: ShellLayout, n: LayoutNode) => boolean;
  /** What `covered` and `away` read besides the layout. */
  sources?: readonly ReadonlyState<unknown>[];
  /**
   * Groups slide in and out as they show and hide (default on: a panel folding beside the main area). Off for the main
   * area's panes: a pane maximized or back shows at once (sliding each pane out in turn drew every layout between).
   */
  motion?: boolean;
  /** After each rebuild (its nodes are new): what is written on them from outside (the active pane's mark). */
  onBuild?: (host: HTMLElement) => void;
}

/**
 * A docked tree (shell §3): groups and splits around the level's leaf. Built only when its structure changes
 * (`sameShape`); everything else is written in place (`applyState`): which nodes show and which are folded, which takes
 * the rest, the sizes, the splitters' sides. Why: rebuilding moves the leaf (the app's canvas), which drops the focus in
 * it and ends a splitter's drag; folding the right group or hiding the rail must not. Each group's element is kept by id
 * while the tree around it changes (its view stays mounted).
 */
function DockTree(o: DockTreeOptions): HTMLElement {
  const { shell, bodies } = o;
  const owner = requireScope('DockTree');
  const host = (<div class={classes.tree} />) as HTMLElement;
  const groups = new Map<string, { el: HTMLElement; scope: ViewScope }>();
  const nodes = new Map<string, HTMLElement>();
  // A group that opens slides in, one that closes slides out before it goes (styles §7.1, `side`); not at the first
  // write, nor when the tree was just rebuilt (its nodes are new).
  let settled = false;
  const closing = new WeakMap<HTMLElement, Animation>();
  const showNode = (el: HTMLElement, n: LayoutNode, shows: boolean) => {
    const was = el.getAttribute('data-shown') !== 'false';
    const pending = closing.get(el);
    if (pending) {
      closing.delete(el);
      pending.cancel();
    }
    if (shows === was || !settled || n.kind !== 'group' || o.motion === false)
      return setAttr(el, 'data-shown', shows ? false : 'false');
    if (shows) {
      setAttr(el, 'data-shown', false);
      el.removeAttribute('data-reveal');
      void el.offsetWidth;
      el.setAttribute('data-reveal', 'side');
      return;
    }
    const anim = Motion.close?.(el, 'side');
    if (!anim) return setAttr(el, 'data-shown', 'false');
    closing.set(el, anim);
    anim.finished.then(
      () => {
        if (closing.get(el) !== anim) return;
        closing.delete(el);
        setAttr(el, 'data-shown', 'false');
        anim.cancel();
      },
      () => {},
    );
  };
  // The last size written to each node (a write only when it changed: a drag's frames touch one node).
  const sizes = new Map<HTMLElement, number>();
  // Whether each node shows, asked once a node per write (a folded group asks its column's, a split its children's).
  const shown = new Map<LayoutNode, boolean>();
  const applyState = (l: ShellLayout) => {
    const ix = layoutIndex(l);
    shown.clear();
    const away = (n: LayoutNode) => o.away?.(l, n) ?? false;
    const shows = (n: LayoutNode): boolean => {
      let v = shown.get(n);
      if (v === undefined) shown.set(n, (v = nodeShows(shell, l, ix, n, shows)));
      return v;
    };
    const walk = (n: LayoutNode, flex: boolean, at: 'start' | 'end' | undefined) => {
      const el = nodes.get(n.id);
      const covered = o.covered?.(l, n) ?? false;
      const show = shows(n) && !away(n);
      if (el) {
        showNode(el, n, show);
        setAttr(el, 'data-folded', n.kind === 'group' && !!n.collapsed && ix.rail?.id !== n.id);
        setAttr(el, 'data-flex', flex);
        setAttr(el, 'data-at', at ?? false);
        // The node that takes the rest has no splitter (hidden for assistive tech too, not only by its style).
        const splitter = el.lastElementChild;
        if (splitter?.getAttribute('role') === 'separator') setAttr(splitter, 'hidden', flex);
        setAttr(el, 'inert', covered);
        if (!flex && (n.kind === 'group' || n.kind === 'split')) {
          const size = n.size ?? DefaultGroupSize;
          if (sizes.get(el) !== size) {
            sizes.set(el, size);
            el.style.setProperty('--msui-shell-size', `${size}px`);
          }
        }
      }
      if (n.kind === 'group') groups.get(n.id)?.scope.setActive(show && !covered);
      if (n.kind !== 'split') return;
      const f = flexIndex(ix, n, o.leaf, shows, away);
      for (let i = 0; i < n.children.length; i++) {
        walk(n.children[i], i === f, i === f ? undefined : i > f ? 'start' : 'end');
      }
    };
    walk(o.rootOf(l), true, undefined);
  };
  const sources = o.sources ?? [];
  // Written again only when what it reads changes (shell §3.1): the tree's nodes, sizes, folds and views, what shows,
  // the rail, and the options' sources; a window's move or the other tree's drag skips it.
  effect(
    derive(
      [shell.layout, shell.absent, shell.views, shell.hiddenViews, ...sources],
      (l, a, v, h, ...rest) => ({ l, a, v, h, rest }),
      (x, y) =>
        x.a === y.a &&
        x.v === y.v &&
        x.h === y.h &&
        x.l.rail === y.l.rail &&
        x.rest.every((r, i) => r === y.rest[i]) &&
        sameTree(o.rootOf(x.l), o.rootOf(y.l)),
    ),
    ({ l }) => {
      applyState(l);
      settled = true;
    },
  );
  let buildScope: ViewScope | undefined;
  // Rebuilt only when its structure changes (compared, not keyed by a string made at each write).
  const shape = derive(
    shell.layout,
    (l) => l,
    (a, b) => railGroup(a)?.id === railGroup(b)?.id && sameShape(o.rootOf(a), o.rootOf(b)),
  );
  effect(shape, () => {
    const l = shell.layout.value;
    buildScope?.dispose();
    buildScope = owner.child();
    const used = new Set<string>();
    nodes.clear();
    sizes.clear();
    const tree = runInScope(buildScope, () =>
      buildNode(shell, o.rootOf(l), undefined, {
        groups,
        used,
        owner,
        bodies,
        leaf: o.leaf,
        leafEl: o.leafEl,
        nodes,
        region: o.leaf === 'content',
      }),
    );
    host.replaceChildren(tree);
    const was = settled;
    settled = false;
    applyState(l);
    settled = was;
    o.onBuild?.(host);
    for (const [id, g] of groups) {
      if (used.has(id)) continue;
      g.scope.dispose();
      g.el.remove();
      groups.delete(id);
    }
  });
  onCleanup(() => {
    for (const g of groups.values()) g.scope.dispose();
  });
  return host;
}

/**
 * Whether a node shows: a group with a view to show (a rail group folded does not: its rail stays), a split with one.
 * `shows` asks it of another node (memoized by the caller: once a node per write).
 */
function nodeShows(
  shell: ShellModel,
  l: ShellLayout,
  ix: LayoutIndex,
  n: LayoutNode,
  shows: (n: LayoutNode) => boolean,
): boolean {
  if (n.kind === 'main' || n.kind === 'content') return true;
  if (n.kind === 'group') {
    if (!hasShown(shell, n)) return false;
    if (n.collapsed && ix.rail?.id === n.id) return false;
    if (!l.rail && ix.rightGroups.has(n.id)) return false;
    // A folded group shows its header only beside an open group of its column: with every group of the column folded,
    // the column closes to the rail (no header left alone over an empty column).
    if (n.collapsed) {
      const column = ix.columns.get(n.id);
      const open = (c: LayoutNode): boolean =>
        c.kind === 'group'
          ? c.id !== n.id && !c.collapsed && hasShown(shell, c)
          : c.kind === 'split' && c.children.some(open);
      return !!column && column.children.some(open);
    }
    return true;
  }
  for (const c of n.children) if (shows(c)) return true;
  return false;
}

/** Whether two trees read the same to `applyState`: nodes, order, direction, sizes, folds, views (not the active one). */
function sameTree(a: LayoutNode, b: LayoutNode): boolean {
  if (a === b) return true;
  if (a.kind !== b.kind || a.id !== b.id) return false;
  if (a.kind === 'group') {
    const g = b as GroupNode;
    return a.size === g.size && !!a.collapsed === !!g.collapsed && sameList(a.views, g.views);
  }
  if (a.kind === 'split') {
    const s = b as SplitNode;
    if (a.dir !== s.dir || a.size !== s.size || a.children.length !== s.children.length) return false;
    for (let i = 0; i < a.children.length; i++) if (!sameTree(a.children[i], s.children[i])) return false;
  }
  // The main area is a leaf of the shell's tree (its region is the other tree's).
  return true;
}

/** Whether two trees have one structure (`buildNode`'s): nodes, order, direction. */
function sameShape(a: LayoutNode, b: LayoutNode): boolean {
  if (a === b) return true;
  if (a.kind !== b.kind || a.id !== b.id) return false;
  if (a.kind !== 'split') return true;
  const s = b as SplitNode;
  if (a.dir !== s.dir || a.children.length !== s.children.length) return false;
  for (let i = 0; i < a.children.length; i++) if (!sameShape(a.children[i], s.children[i])) return false;
  return true;
}

/**
 * The child of a split that takes the rest: the one holding the main area, else the first that shows open with no
 * size, else the first that shows open (a folded group gives its height to the other).
 */
function flexIndex(
  ix: LayoutIndex,
  n: SplitNode,
  leaf: 'main' | 'content',
  shows: (n: LayoutNode) => boolean,
  away: (n: LayoutNode) => boolean,
): number {
  const open = (c: LayoutNode) => shows(c) && !(c.kind === 'group' && c.collapsed) && !away(c);
  // The leaf takes the rest, unless it is away (another pane maximized): then the pane that shows does.
  let flex = n.children.findIndex((c) => holds(ix, c, leaf) && !away(c));
  if (flex < 0) flex = n.children.findIndex((c) => (c as GroupNode).size === undefined && open(c));
  if (flex < 0) flex = n.children.findIndex(open);
  return Math.max(0, flex);
}

interface BuildContext {
  readonly groups: Map<string, { el: HTMLElement; scope: ViewScope }>;
  readonly used: Set<string>;
  readonly owner: ViewScope;
  readonly bodies: ViewBodies;
  readonly leaf: 'main' | 'content';
  readonly leafEl: HTMLElement;
  readonly nodes: Map<string, HTMLElement>;
  /** The main area's region: its groups and its content are panes. */
  readonly region: boolean;
}

/** A node's element (and its children's); what changes without the structure is `applyState`'s. */
function buildNode(
  shell: ShellModel,
  n: LayoutNode,
  parentDir: 'row' | 'column' | undefined,
  ctx: BuildContext,
): HTMLElement {
  const { groups, used, owner, bodies } = ctx;
  const wrap = (<div class={classes.node} data-node={n.id} />) as HTMLElement;
  ctx.nodes.set(n.id, wrap);
  if (ctx.region && (n.kind === 'group' || n.kind === 'content')) wrap.setAttribute('data-pane', n.id);
  if (n.kind === 'main' || n.kind === 'content') {
    // The leaf takes the rest: no splitter of its own (a sized neighbour's is on the neighbour).
    if (n.kind === ctx.leaf) wrap.append(ctx.leafEl);
    return wrap;
  }
  if (n.kind === 'group') {
    used.add(n.id);
    let g = groups.get(n.id);
    if (!g) {
      const scope = owner.child();
      const el = runInScope(scope, () =>
        ctx.region
          ? ((<PaneGroupView shell={shell} id={n.id} bodies={bodies} />) as HTMLElement)
          : ((<GroupView shell={shell} id={n.id} bodies={bodies} />) as HTMLElement),
      );
      g = { el, scope };
      groups.set(n.id, g);
    }
    wrap.append(g.el);
  } else {
    const split = (<div class={classes.split} data-dir={n.dir} />) as HTMLElement;
    for (const c of n.children) split.append(buildNode(shell, c, n.dir, ctx));
    wrap.append(split);
  }
  // Every node in a split has its splitter; it shows on a sized node (`data-at`: the side toward the one that takes
  // the rest), written in place.
  if (parentDir)
    wrap.append(
      <Splitter shell={shell} id={n.id} dir={parentDir} wrap={wrap} min={ctx.region ? MinPaneSize : MinGroupSize} />,
    );
  return wrap;
}

/** A splitter on a sized node's edge toward the node that takes the rest: drags or arrows resize it. */
function Splitter(props: { shell: ShellModel; id: string; dir: 'row' | 'column'; wrap: HTMLElement; min: number }) {
  const { shell, id, dir } = props;
  // Its side (toward the node that takes the rest) is its node's, written in place: read when it is used.
  const sign = () => (props.wrap.dataset.at === 'start' ? -1 : 1);
  const size = () => {
    const r = props.wrap.getBoundingClientRect();
    return dir === 'row' ? r.width : r.height;
  };
  const el = (
    <div
      class={classes.splitter}
      data-dir={dir}
      role="separator"
      tabIndex={0}
      aria-orientation={dir === 'row' ? 'vertical' : 'horizontal'}
      aria-label={shell.labels.resize}
      aria-valuemin={String(props.min)}
      onKeyDown={(e: KeyboardEvent) => {
        const step = e.shiftKey ? 64 : 16;
        const keys = dir === 'row' ? ['ArrowLeft', 'ArrowRight'] : ['ArrowUp', 'ArrowDown'];
        const i = keys.indexOf(e.key);
        if (i < 0) return;
        e.preventDefault();
        shell.resize(id, size() + (i === 0 ? -1 : 1) * step * sign());
      }}
    />
  ) as HTMLElement;
  const next = frameState(0);
  watch([next], () => {
    if (next.value > 0) shell.resize(id, next.value);
  });
  // Its value is the node's size (a focusable separator says where it is).
  const value = derive(shell.layout, (l) => {
    const n = nodeOf(l, id);
    return n?.kind === 'group' || n?.kind === 'split' ? n.size : undefined;
  });
  watch([value], () => setAttr(el, 'aria-valuenow', value.value === undefined ? false : String(value.value)));
  // In a lined-up grid, the cells at its place in the other lines move with it: their splitters light with it (hovered,
  // focused, dragged), so the line across reads as one.
  let lit: Element[] = [];
  const light = (on: boolean) => {
    for (const p of lit) p.removeAttribute('data-linked');
    lit = [];
    if (!on) return;
    const tree = props.wrap.closest(`.${classes.tree}`);
    for (const peer of gridPeers(shell.layout.value, id)) {
      const s = tree?.querySelector(`[data-node="${peer}"] > [role="separator"]`);
      if (s) {
        s.setAttribute('data-linked', '');
        lit.push(s);
      }
    }
  };
  let hovered = false;
  let focused = false;
  let dragging = false;
  const relight = () => light(hovered || focused || dragging);
  on(el, 'pointerenter', () => {
    hovered = true;
    relight();
  });
  on(el, 'pointerleave', () => {
    hovered = false;
    relight();
  });
  on(el, 'focus', () => {
    focused = true;
    relight();
  });
  on(el, 'blur', () => {
    focused = false;
    relight();
  });
  onCleanup(() => light(false));
  drag(el, {
    start: (e) => {
      e.preventDefault();
      dragging = true;
      relight();
      el.setAttribute('data-dragging', '');
      return { from: dir === 'row' ? e.clientX : e.clientY, size: size(), sign: sign() };
    },
    move: (e, s) => next.offer(s.size + ((dir === 'row' ? e.clientX : e.clientY) - s.from) * s.sign),
    end: () => {
      el.removeAttribute('data-dragging');
      dragging = false;
      relight();
    },
  });
  return el;
}

/** A docked group: its view's frame, with the group's tabs in its header when it is not the rail's. */
function GroupView(props: { shell: ShellModel; id: string; bodies: ViewBodies }) {
  const { shell, id, bodies } = props;
  const node = groupNode(shell, (l) => findGroup(l, id)?.group);
  const isRail = derive(shell.layout, (l) => railGroup(l)?.id === id);
  const state = derive(
    [node, isRail, shell.absent, shell.views, shell.hiddenViews],
    (g, rail) => ({ active: activeOf(shell, g), rail, views: shownViews(shell, g) }),
    (a, b) => a.active === b.active && a.rail === b.rail && sameList(a.views, b.views),
  );
  const key = derive(state, (s) => `${s.active}|${s.rail}|${s.views.join(',')}`);
  const open = derive(node, (g) => !g?.collapsed);
  const enter = entering();
  return (
    // A row's tooltip card goes outside the group, toward the roomier side (the main area).
    <div class={classes.group} data-group={id} data-tip-edge="" data-msui-region="">
      <Show when={key}>
        {() => {
          const s = state.value;
          if (!s.active) return null;
          const ref = s.active;
          const tabs =
            !s.rail && s.views.length > 1 ? (
              <Tabs shell={shell} views={s.views} active={ref} menu={(r) => viewMenu(shell, r)} />
            ) : undefined;
          let header!: HTMLElement;
          return enter(
            <KeptFrame
              shell={shell}
              bodies={bodies}
              ref={ref}
              placement="docked"
              head={s.rail ? 'rail' : 'header'}
              tabs={tabs}
              open={open}
              headerRef={(el) => (header = el)}
              // ✕: the rail's group folds (its rail stays to reopen it); a group under it closes the split, its views
              // back into the rail's group (none folded out of reach).
              onClose={() => {
                const rail = railGroup(shell.layout.value);
                if (!rail || rail.id === id) shell.collapse(id);
                else shell.closeGroup(id);
              }}
              menu={() => viewMenu(shell, ref, { anchorRect: () => anchorIn(header, shell, bodies.ui.main) })}
            />,
          );
        }}
      </Show>
    </div>
  );
}

/**
 * A group in the main area's region (shell §2.1): a pane. A canvas view alone has no header (no chrome above a canvas;
 * its pane's menu is the canvas's right-click, the app's); several views show their tabs; a panel shows its header. Its
 * ✕ closes its views (they leave the region; shown again they go home).
 */
function PaneGroupView(props: { shell: ShellModel; id: string; bodies: ViewBodies }) {
  const { shell, id, bodies } = props;
  const state = shownState(
    shell,
    groupNode(shell, (l) => findGroup(l, id)?.group),
  );
  const key = derive(state, (s) => `${s.active}|${s.views.join(',')}`);
  const enter = entering();
  return (
    <div class={classes.group} data-group={id} data-tip-edge="" data-msui-region="">
      <Show when={key}>
        {() => {
          const s = state.value;
          if (!s.active) return null;
          const ref = s.active;
          const def = shell.view(ref);
          const tabs =
            s.views.length > 1 ? (
              <Tabs shell={shell} views={s.views} active={ref} menu={(r) => paneMenu(shell, r)} />
            ) : undefined;
          let header!: HTMLElement;
          return enter(
            <KeptFrame
              shell={shell}
              bodies={bodies}
              ref={ref}
              placement="stage"
              head={def?.canvas && !tabs ? 'none' : 'header'}
              tabs={tabs}
              headerRef={(el) => (header = el)}
              onClose={() => {
                for (const v of findGroup(shell.layout.value, id)?.group.views ?? []) shell.close(v);
              }}
              menu={() => paneMenu(shell, ref, { anchorRect: () => anchorIn(header, shell, bodies.ui.main) })}
            />,
          );
        }}
      </Show>
    </div>
  );
}

/**
 * The top overlay (shell §2.1): the views the layout puts on the top of the main area's content, for the app to place
 * (`main`'s `top`: a stage's band slot). The active one's body, no header (its tabs when there are several).
 */
function TopOverlay(props: { shell: ShellModel; bodies: ViewBodies }) {
  const { shell, bodies } = props;
  const state = shownState(
    shell,
    groupNode(shell, (l) => mainNode(l).top),
  );
  const key = derive(state, (s) => `${s.active}|${s.views.join(',')}`);
  return (
    <div class={classes.topOverlay} data-group={TopGroup} data-msui-region="">
      <Show when={key}>
        {() => {
          const s = state.value;
          if (!s.active) return null;
          const tabs =
            s.views.length > 1 ? (
              <Tabs shell={shell} views={s.views} active={s.active} menu={(r) => paneMenu(shell, r)} />
            ) : undefined;
          return (
            <KeptFrame
              shell={shell}
              bodies={bodies}
              ref={s.active}
              placement="top"
              head={tabs ? 'header' : 'none'}
              tabs={tabs}
              onClose={() => shell.close(s.active!)}
              menu={() => paneMenu(shell, s.active!)}
            />
          );
        }}
      </Show>
    </div>
  );
}

/** A rectangle of `el` relative to the shell's main area (where a window opens beside it). */
function anchorIn(el: Element | undefined, shell: ShellModel, main: HTMLElement): Rect | undefined {
  if (!el) return undefined;
  const a = main.getBoundingClientRect();
  const r = el.getBoundingClientRect();
  const area = shell.area.value;
  return { x: Math.min(r.left - a.left, area.width), y: Math.max(0, r.top - a.top), width: r.width, height: r.height };
}

/** A group's value, a new one only when it reads differently (a drag's frames and a window's move write none). */
function groupNode(
  shell: ShellModel,
  of: (l: ShellLayout) => GroupNode | undefined,
): ReadonlyState<GroupNode | undefined> {
  return derive(shell.layout, of, sameGroup);
}

/** The views of a group that show and the one it shows, a new value only when they change. */
function shownState(
  shell: ShellModel,
  node: ReadonlyState<GroupNode | undefined>,
): ReadonlyState<{ active: ViewRef | undefined; views: ViewRef[] }> {
  return derive(
    [node, shell.absent, shell.views, shell.hiddenViews],
    (g) => ({ active: activeOf(shell, g), views: shownViews(shell, g) }),
    (a, b) => a.active === b.active && sameList(a.views, b.views),
  );
}

/**
 * A group's labelled tabs (a group not at the rail, a window of several views): a tab list, the arrows between. A tab's
 * context menu (or Shift+F10) is `menu`'s: the view's where it is docked or in a window, the pane's in the main area.
 */
function Tabs(props: {
  shell: ShellModel;
  views: readonly ViewRef[];
  active: ViewRef;
  menu: (ref: ViewRef) => ContextMenuItem[];
}) {
  const { shell } = props;
  const list = (
    <div class={classes.tabs} role="tablist">
      {props.views.map((ref) => {
        const def = shell.view(ref)!;
        return (
          <button
            type="button"
            role="tab"
            class={classes.tab}
            aria-selected={String(ref === props.active)}
            tabIndex={ref === props.active ? 0 : -1}
            onClick={() => shell.select(ref)}
            onContextMenu={(e: MouseEvent) => {
              e.preventDefault();
              openContextMenu(menuRoot(e.currentTarget as Element), { x: e.clientX, y: e.clientY }, props.menu(ref));
            }}
            onKeyDown={(e: KeyboardEvent) => {
              if ((e.key === 'F10' && e.shiftKey) || e.key === 'ContextMenu') {
                e.preventDefault();
                const r = (e.currentTarget as Element).getBoundingClientRect();
                openContextMenu(menuRoot(e.currentTarget as Element), { x: r.left, y: r.bottom }, props.menu(ref));
              }
            }}
          >
            <Icon name={def.icon} />
            {def.heading && typeof def.heading === 'string' ? def.heading : shell.paneTitle(ref)}
          </button>
        );
      })}
    </div>
  ) as HTMLElement;
  rovingFocus(list, {
    items: '[role="tab"]',
    axis: 'horizontal',
    preferred: (all) => all.find((t) => t.ariaSelected === 'true'),
  });
  return list;
}

function menuRoot(el: Element): HTMLElement {
  return el.closest<HTMLElement>('.msui-root') ?? (el.ownerDocument.body as HTMLElement);
}

/**
 * The rail (shell §3.2): the rail group's views as a vertical tab list (head, strip, foot), the app's rail items by group,
 * at the shell's right edge, full height. Its buttons step down with the shell's height; what does not fit goes under
 * More. A press on the open view's button folds the column; a button's context menu (or Shift+F10) hides it, shows a
 * hidden one, resets the layout.
 */
function Rail(props: { shell: ShellModel; hideable: boolean }) {
  const { shell } = props;
  const t = shell.labels;
  const height = new State(0);
  const railNode = groupNode(shell, railGroup);
  // The rail's buttons, a new list only when they change (not at each layout write: a drag's sizes, a window's move).
  const entries = derive(
    [railNode, shell.absent, shell.views, shell.railItems, shell.hiddenViews],
    (g, absent, _, items) => {
      const views = shownViews(shell, g).map((ref) => ({ ref, def: shell.view(ref)! }));
      const tools = items.filter((i) => !absent.has(`item:${i.id}`));
      return railEntries(views, tools);
    },
    sameEntries,
  );
  // A button is its entry's id and what it shows (a view's def, an item): a view registered anew under its id (another
  // title or icon) is a new button.
  const versions = new WeakMap<object, number>();
  let version = 0;
  const entryKey = (e: RailEntry) => {
    const o = e.view?.def ?? e.item!;
    let v = versions.get(o);
    if (v === undefined) versions.set(o, (v = ++version));
    return `${e.id}|${v}`;
  };
  const fit = derive(
    [entries, height],
    (e, h) => fitRail(e, h),
    (a, b) => a.step === b.step && a.hidden.size === b.hidden.size && [...a.hidden].every((id) => b.hidden.has(id)),
  );
  const active = derive([railNode, shell.absent, shell.views, shell.hiddenViews], (g) =>
    g && !g.collapsed ? activeOf(shell, g) : undefined,
  );
  const rail = (<nav class={classes.rail} aria-label={t.rail} data-msui-region="" />) as HTMLElement;
  observeResize(rail, () => {
    const h = Math.round(rail.clientHeight);
    if (h !== height.value) height.set(h);
  });
  rovingFocus(rail, { items: '[data-rail-button]', axis: 'vertical', attributes: ['aria-selected'] });
  const section = (name: string) => (es: readonly RailEntry[]) => es.filter((e) => e.section === name);
  // A rail item's press: its `run`, else its menu beside the button (aria-expanded while open).
  const press = (item: RailItem, anchor: Element) => {
    if (item.run) return item.run(anchor);
    if (!item.menu) return;
    anchor.setAttribute('aria-expanded', 'true');
    openMenuBeside(menuRoot(rail), anchor, item.menu(), 'left', () => anchor.setAttribute('aria-expanded', 'false'));
  };
  const button = (e: RailEntry) => {
    if (e.view) {
      const { ref, def } = e.view;
      const el = (
        <button
          type="button"
          role="tab"
          class={classes.railButton}
          data-rail-button=""
          data-tip-side="left"
          aria-label={def.title}
          data-tip={def.title}
          onClick={() => shell.select(ref)}
          onContextMenu={(ev: MouseEvent) => {
            ev.preventDefault();
            openContextMenu(menuRoot(rail), { x: ev.clientX, y: ev.clientY }, railMenu(shell, ref));
          }}
          onKeyDown={(ev: KeyboardEvent) => {
            if ((ev.key === 'F10' && ev.shiftKey) || ev.key === 'ContextMenu') {
              ev.preventDefault();
              const r = (ev.currentTarget as Element).getBoundingClientRect();
              openContextMenu(menuRoot(rail), { x: r.left, y: r.bottom }, railMenu(shell, ref));
            }
          }}
        >
          <Icon name={def.icon} />
        </button>
      ) as HTMLElement;
      watch([active, fit], () => {
        setAttr(el, 'aria-selected', String(active.value === ref));
        setAttr(el, 'hidden', fit.value.hidden.has(e.id));
      });
      return el;
    }
    const item = e.item!;
    const el = (
      <button
        type="button"
        class={classes.railButton}
        data-rail-button=""
        data-tip-side="left"
        data-captioned={item.caption ? '' : undefined}
        aria-label={item.title}
        data-tip={item.title}
        aria-haspopup={item.menu && !item.run ? 'menu' : undefined}
        onClick={(ev: MouseEvent) => press(item, ev.currentTarget as Element)}
      >
        <Icon name={item.icon} />
        {item.caption ? <span class={classes.railCaption}>{item.caption}</span> : null}
      </button>
    ) as HTMLElement;
    watch([fit], () => setAttr(el, 'hidden', fit.value.hidden.has(e.id)));
    if (item.pressed) watch([item.pressed], () => setAttr(el, 'aria-pressed', String(!!item.pressed!.value)));
    return el;
  };
  const run = (name: string, tablist = false) => {
    const list = derive(entries, section(name), sameEntries);
    // A section whose buttons are all under More draws nothing (no hairline of its own).
    const any = derive([list, fit], (l, f) => l.some((e) => !f.hidden.has(e.id)));
    return (
      <Show when={any}>
        {() => (
          <div
            class={classes.railSection}
            role={tablist ? 'tablist' : 'group'}
            data-rail-section={name}
            aria-orientation={tablist ? 'vertical' : undefined}
            // Named for people (its views, or tools), not by its internal section (`head`, `group:app`).
            aria-label={tablist ? shell.labels.railViews : shell.labels.railTools}
          >
            <For each={list} by={entryKey}>
              {(e) => button(e)}
            </For>
          </div>
        )}
      </Show>
    );
  };
  const groupsOf = derive(entries, (es) => [
    ...new Set(es.filter((e) => e.section.startsWith('group:')).map((e) => e.section)),
  ]);
  let moreButton!: HTMLElement;
  const overflow = derive([entries, fit], (es, f) => es.filter((e) => f.hidden.has(e.id)));
  rail.append(
    (
      <>
        {run('head', true)}
        {run('strip', true)}
        <For each={groupsOf} by={(g) => g}>
          {(g) => run(g)}
        </For>
        <div class={classes.railFoot}>
          <Show when={derive(overflow, (o) => o.length > 0)}>
            {() => (
              <div class={classes.railSection}>
                <button
                  type="button"
                  class={classes.railButton}
                  data-rail-button=""
                  data-tip-side="left"
                  aria-label={t.more}
                  data-tip={t.more}
                  aria-haspopup="menu"
                  ref={(el: HTMLElement) => (moreButton = el)}
                  onClick={() => {
                    openMenuBeside(
                      menuRoot(rail),
                      moreButton,
                      overflow.value.map((e) => ({
                        id: e.id,
                        label: e.view ? e.view.def.title : e.item!.title,
                        icon: e.view ? e.view.def.icon : e.item!.icon,
                        run: () => (e.view ? shell.select(e.view.ref) : press(e.item!, moreButton)),
                      })),
                    );
                  }}
                >
                  <Icon name="ellipsis" />
                </button>
              </div>
            )}
          </Show>
          {run('foot', true)}
          {run('tools')}
          {props.hideable ? (
            <div class={classes.railSection}>
              <button
                type="button"
                class={classes.railButton}
                data-rail-button=""
                data-tip-side="left"
                aria-label={t.hideViews}
                data-tip={t.hideViews}
                onClick={() => shell.toggleRail(false)}
              >
                <Icon name="chevrons-right" />
              </button>
            </div>
          ) : null}
          {run('corner')}
        </div>
      </>
    ) as Node,
  );
  return rail;
}

/** The rail's toggle (the app's mark): in the main area's corner, it shows the rail and the docked column again. */
function RailToggle(props: { shell: ShellModel; compact: ReadonlyState<boolean>; icon: IconName }) {
  const { shell } = props;
  const shown = derive([shell.railShown, props.compact, shell.expanded], (r, c, e) => !c && !r && !e);
  return (
    <Show when={shown}>
      {() => (
        <div class={classes.railToggle}>
          <Button
            variant="ghost"
            icon={<Icon name={props.icon} />}
            title={shell.labels.showViews}
            onClick={() => shell.toggleRail(true)}
          />
        </div>
      )}
    </Show>
  );
}

// ---------------------------------------------------------------------------------------------------------------------
// Windows

function WindowsLayer(props: { shell: ShellModel; bodies: ViewBodies }) {
  const { shell } = props;
  // In a stable order (by id, not by stacking): a raise must not move a window's node (moving it would drop the focus
  // in it, reload a frame in it); the stacking is each window's `z-index`.
  const ids = derive(
    shell.layout,
    (l) => l.windows.map((w) => w.group.id).sort(),
    (a, b) => a.length === b.length && a.every((x, i) => x === b[i]),
  );
  // Windows there at the first render stay still; one opened later floats in.
  let ready = false;
  queueMicrotask(() => (ready = true));
  // Away while expanded, or under a window's maximized pane: out of the focus order and the pointer's way (hidden by
  // the style while expanded; each frame pauses).
  const away = derive(
    [shell.expanded, shell.maximized, shell.layout],
    (e, m, l) => e || (!!m && !!findView(l, m)?.window),
  );
  return (
    <div class={classes.windows} ref={(el: HTMLElement) => watch([away], () => setAttr(el, 'inert', away.value))}>
      <For each={ids} by={(id) => id}>
        {(id) => <WindowView shell={shell} bodies={props.bodies} id={id} enter={ready} />}
      </For>
    </div>
  );
}

/**
 * A window (shell §3.4): its header drags it, its right and bottom edges and corner resize it (in the frame lane), kept
 * inside the main area; a press or focus in it brings it to the top. ✕ or Escape (the focused window's) closes its
 * views; a press outside does not. Its menu docks it, moves it into another window, or moves and resizes it from the
 * keyboard (arrows by 4px, Shift by 16, Enter or Escape ends).
 */
function WindowView(props: { shell: ShellModel; bodies: ViewBodies; id: string; enter?: boolean }) {
  const { shell, id } = props;
  // Its rectangle and place in the stack, a new value only when they change (not at another window's move).
  const win = derive(
    shell.layout,
    (l) => l.windows.find((w) => w.group.id === id),
    (a, b) =>
      a === b ||
      (!!a &&
        !!b &&
        a.z === b.z &&
        a.rect.x === b.rect.x &&
        a.rect.y === b.rect.y &&
        a.rect.width === b.rect.width &&
        a.rect.height === b.rect.height),
  );
  const view = derive(
    [
      shownState(
        shell,
        groupNode(shell, (l) => findGroup(l, id)?.group),
      ),
      shell.maximized,
    ],
    (s, maximized) => {
      // Its view maximized (shell §2.1): out over the main area (`MaximizedHost`), its body back here after.
      const out = !!s.active && s.active === maximized;
      return { active: out ? undefined : s.active, views: s.views, key: `${s.active}|${s.views.join(',')}|${out}` };
    },
  );
  const key = derive(view, (v) => v.key);
  const notExpanded = derive(
    [shell.expanded, shell.maximized, shell.layout],
    (e, m, l) => !e && !(m && findView(l, m)?.window),
  );
  const adjusting = new State<WindowAdjust | undefined>(undefined);
  const enterView = entering();
  const rect = frameState<Rect | undefined>(undefined);
  watch([rect], () => {
    if (rect.value) shell.setWindowRect(id, rect.value);
  });
  const el = (
    <div
      class={classes.window}
      role="dialog"
      aria-modal="false"
      data-window={id}
      data-tip-edge=""
      data-msui-region=""
      data-enter={props.enter ? '' : undefined}
    />
  ) as HTMLElement;
  watch([win, adjusting], () => {
    const w = win.value;
    if (!w) return;
    el.style.left = `${w.rect.x}px`;
    el.style.top = `${w.rect.y}px`;
    el.style.width = `${w.rect.width}px`;
    el.style.height = `${w.rect.height}px`;
    el.style.zIndex = String(w.z);
    setAttr(el, 'data-adjusting', adjusting.value ?? false);
  });
  // Its name, and nothing on screen when none of its views can show (unregistered, flagged off, hidden): not an empty
  // box with no header to close it by.
  // Named by its pane's title, followed as it changes (`setTitle`, shell §2.1).
  let stopTitle: (() => void) | undefined;
  watch([key], () => {
    const active = view.value.active;
    setAttr(el, 'hidden', !active);
    stopTitle?.();
    stopTitle = undefined;
    if (!active) return setAttr(el, 'aria-label', shell.labels.window);
    const title = shell.paneTitle(active);
    setAttr(el, 'aria-label', title.value);
    stopTitle = listen(title, () => setAttr(el, 'aria-label', title.value));
  });
  onCleanup(() => stopTitle?.());
  on(el, 'pointerdown', () => shell.raise(id), { capture: true });
  // A window's view is a pane too (a detached stage): pressed or focused, it is the active one.
  const activate = () => {
    const active = view.value.active;
    if (active) shell.setActivePane(active);
  };
  on(el, 'pointerdown', activate, { capture: true });
  on(el, 'focusin', activate);
  on(el, 'focusin', () => shell.raise(id));
  // Its pane's status (shell §2.1): hovered and focused, fed as the region's are.
  paneInputs(shell, el, () => view.value.active);
  const closeAll = () => {
    for (const v of findGroup(shell.layout.value, id)?.group.views ?? []) shell.close(v);
  };
  dismissable({
    layer: () => el,
    close: (reason) => {
      if (reason === 'escape') {
        if (adjusting.value) adjusting.set(undefined);
        else closeAll();
      }
    },
    outside: false,
    ownsKey: () => false,
    refocus: () => props.bodies.ui.main,
  });
  on<KeyboardEvent>(el, 'keydown', (e) => {
    const how = adjusting.value;
    const w = win.value;
    if (!how || !w) return;
    if (e.key === 'Enter') {
      e.preventDefault();
      adjusting.set(undefined);
      return;
    }
    const step = e.shiftKey ? 16 : 4;
    const d = { ArrowLeft: [-step, 0], ArrowRight: [step, 0], ArrowUp: [0, -step], ArrowDown: [0, step] }[e.key];
    if (!d) return;
    e.preventDefault();
    const r = w.rect;
    shell.setWindowRect(
      id,
      how === 'move'
        ? { ...r, x: r.x + d[0], y: r.y + d[1] }
        : { ...r, width: r.width + d[0], height: r.height + d[1] },
    );
  });
  const resized = (
    s: { x: number; y: number; r: Rect },
    e: PointerEvent,
    edge?: 'right' | 'bottom' | 'corner',
  ): Rect => ({
    ...s.r,
    width: edge === 'bottom' ? s.r.width : s.r.width + e.clientX - s.x,
    height: edge === 'right' ? s.r.height : s.r.height + e.clientY - s.y,
  });
  const grip = (edge: 'right' | 'bottom' | 'corner') => {
    const g = (<div class={classes.grip} data-edge={edge} aria-hidden="true" />) as HTMLElement;
    drag(g, {
      start: (e) => {
        const w = win.value;
        if (!w) return undefined;
        e.preventDefault();
        return { x: e.clientX, y: e.clientY, r: w.rect };
      },
      move: (e, s) => rect.offer(resized(s, e, edge)),
      end: (s, how, e) => {
        // The release's place, through the frame state (a commit beside it would be overwritten by its pending value).
        if (how === 'up' && e) rect.offer(resized(s, e, edge));
      },
    });
    return g;
  };
  el.append(
    (
      <>
        <Show when={key}>
          {() => {
            const v = view.value;
            if (!v.active) return null;

            const ref = v.active;
            let header!: HTMLElement;
            const frame = (
              <KeptFrame
                shell={shell}
                bodies={props.bodies}
                ref={ref}
                placement="window"
                tabs={
                  v.views.length > 1 ? (
                    <Tabs shell={shell} views={v.views} active={ref} menu={(r) => viewMenu(shell, r)} />
                  ) : undefined
                }
                headerRef={(h) => (header = h)}
                onClose={closeAll}
                menu={() =>
                  viewMenu(shell, ref, {
                    adjust: (how) => {
                      adjusting.set(how);
                      el.focus();
                    },
                    resetPosition: () => {
                      const def = shell.view(ref);
                      const area = shell.area.value;
                      const size = def?.size ?? DefaultWindowSize;
                      shell.setWindowRect(id, {
                        x: Math.max(0, (area.width - size.width) / 2),
                        y: Math.max(0, (area.height - size.height) / 2),
                        ...size,
                      });
                    },
                  })
                }
              />
            );
            drag(header, {
              start: (e) => {
                // A press on a button or a tab in the header is that control's.
                if ((e.target as Element).closest('button, [role="tab"], input')) return undefined;
                const w = win.value;
                if (!w) return undefined;
                // Not a text selection.
                e.preventDefault();
                return { x: e.clientX, y: e.clientY, r: w.rect };
              },
              move: (e, s) => rect.offer({ ...s.r, x: s.r.x + e.clientX - s.x, y: s.r.y + e.clientY - s.y }),
              end: (s, how, e) => {
                if (how === 'up' && e) rect.offer({ ...s.r, x: s.r.x + e.clientX - s.x, y: s.r.y + e.clientY - s.y });
              },
            });
            // Away while expanded (shell §3.9): hidden and paused; a window whose view is out (expanded) has dropped it.
            return enterView(activity(notExpanded, () => frame) as Node);
          }}
        </Show>
        {grip('right')}
        {grip('bottom')}
        {grip('corner')}
      </>
    ) as Node,
  );
  el.tabIndex = -1;
  return el;
}

/**
 * A pane's hover and focus outside the region (a window's view, a window's pane maximized), fed to the model as the
 * region's are (shell §2.1): the pointer over `el`, the focus in it.
 */
function paneInputs(shell: ShellModel, el: HTMLElement, pane: () => string | undefined): void {
  on(el, 'pointerover', () => shell.setPaneHovered(pane()));
  on(el, 'pointerleave', () => shell.setPaneHovered(undefined));
  on(el, 'focusin', () => shell.setPaneFocused(pane()));
  on(el, 'focusout', (e: FocusEvent) => {
    if (!el.contains(e.relatedTarget as Node | null)) shell.setPaneFocused(undefined);
  });
}

/**
 * A window's pane maximized (shell §2.1): its body out of the window, over the main area (the shell's chrome stays),
 * its frame without a header for a stage (a panel keeps its header, its ✕ shows every pane again); back in its window
 * after.
 */
function MaximizedHost(props: { shell: ShellModel; bodies: ViewBodies }) {
  const { shell, bodies } = props;
  const windowed = derive([shell.maximized, shell.layout], (m, l) => (m && findView(l, m)?.window ? m : undefined));
  return (
    <Show when={windowed}>
      {(ref) => (
        <div
          class={classes.maximizedPane}
          data-maximized-pane={ref}
          ref={(el: HTMLElement) => paneInputs(shell, el, () => ref)}
        >
          <KeptFrame
            shell={shell}
            bodies={bodies}
            ref={ref}
            placement="stage"
            head={shell.view(ref)?.canvas ? 'none' : 'header'}
            onClose={() => shell.maximize(undefined)}
          />
        </div>
      )}
    </Show>
  );
}

// ---------------------------------------------------------------------------------------------------------------------
// Layers: the card, the page, the modal

/**
 * The card (shell §3.6): one at a time, in the right column's place (over its groups, the rail beside it), full height; it
 * closes on its ✕, Escape, or a press outside it, the main area included (not a persistent card).
 */
function CardHost(props: { shell: ShellModel; bodies: ViewBodies; into: HTMLElement }) {
  const { shell } = props;
  const card = shell.card;
  props.into.append(
    (
      <Show when={card} reveal="side">
        {(open) => {
          let el!: HTMLElement;
          const restore = rememberFocus(() => open.anchor as HTMLElement | undefined);
          onCleanup(() => restore(el));
          dismissable({
            layer: () => el,
            anchor: open.anchor,
            close: () => shell.close(open.ref),
            outside: !open.def.persistent,
            // Escape anywhere is the card's, but not with the focus in a window or a dialog: that one's (shell §3.6).
            ownsKey: (e) => !(eventTarget(e) as Element | null)?.closest?.('[data-window], [role="dialog"]'),
            refocus: () => open.anchor as HTMLElement | undefined,
          });
          el = (
            <div class={classes.card} data-card={open.ref} data-tip-edge="" data-msui-region="">
              <KeptFrame
                shell={shell}
                bodies={props.bodies}
                ref={open.ref}
                placement="card"
                onClose={() => shell.close(open.ref)}
              />
            </div>
          ) as HTMLElement;
          // Into the card when it was opened from the keyboard (its first control); else the card itself.
          const byKeyboard = focusedByKeyboard(open.anchor);
          queueMicrotask(() => {
            if (el.isConnected && !el.contains(deepActiveElement() ?? null)) focusInitial(el, byKeyboard);
          });
          return el;
        }}
      </Show>
    ) as Node,
  );
  return null;
}

/** The page (shell §3.6): opaque over the main area, its own header; Escape or ✕ close it, the focus goes back. */
function PageHost(props: { shell: ShellModel; bodies: ViewBodies }) {
  const { shell } = props;
  return (
    <Show when={derive(shell.page, (p) => p)} reveal="fade">
      {(open) => {
        let el!: HTMLElement;
        const restore = rememberFocus(() => open.anchor as HTMLElement | undefined);
        onCleanup(() => restore(el));
        dismissable({ layer: () => el, close: () => shell.close(open.ref), outside: false, ownsKey: () => true });
        el = (
          <div class={classes.page} role="region" aria-label={open.def.title} tabIndex={-1}>
            <KeptFrame
              shell={shell}
              bodies={props.bodies}
              ref={open.ref}
              placement="page"
              onClose={() => shell.close(open.ref)}
            />
          </div>
        ) as HTMLElement;
        queueMicrotask(() => el.isConnected && el.focus({ preventScroll: true }));
        return el;
      }}
    </Show>
  );
}

/** The modal (shell §3.6): over the whole shell with a scrim; Tab stays in it, the rest is inert; Escape or ✕ close it. */
function ModalHost(props: { shell: ShellModel; bodies: ViewBodies }) {
  const { shell } = props;
  return (
    <Show when={shell.modal} reveal="fade">
      {(open) => {
        let el!: HTMLElement;
        const restore = rememberFocus(() => open.anchor as HTMLElement | undefined);
        onCleanup(() => restore(el));
        dismissable({ layer: () => el, close: () => shell.close(open.ref), outside: false, ownsKey: () => true });
        const scrim = (
          <div class={classes.scrim}>
            <div
              class={classes.modal}
              role="dialog"
              aria-modal="true"
              aria-label={open.def.title}
              ref={(m: HTMLElement) => (el = m)}
              onKeyDown={(e: KeyboardEvent) => trapTab(e, el)}
            >
              <KeptFrame
                shell={shell}
                bodies={props.bodies}
                ref={open.ref}
                placement="modal"
                onClose={() => shell.close(open.ref)}
              />
            </div>
          </div>
        ) as HTMLElement;
        // The rest of the shell is inert while it is open (the shell's layout wrapper, `ShellView`).
        queueMicrotask(() => el.isConnected && focusInitial(el, true));
        return scrim;
      }}
    </Show>
  );
}

// ---------------------------------------------------------------------------------------------------------------------
// Busy and file drop

/**
 * What holds the app (shell §3.6, busy): a scrim and one card, the rest inert (the shell's layout wrapper). A dialog:
 * the focus goes to it (its Cancel, else the card) and stays in it, and goes back when it ends. It follows the current
 * task: its label, its progress, and its Cancel (shown while the current one has `cancel`, and cancelling that one).
 */
function BusyOverlay(props: { busy: ReadonlyState<ShellBusy | undefined>; cancel: string }) {
  // A new task replacing one (another label) keeps the overlay: shown while there is one.
  const some = derive(props.busy, (b) => !!b);
  const cancellable = derive(props.busy, (b) => !!b?.cancel);
  return (
    <Show when={some} reveal="fade">
      {() => {
        let card!: HTMLElement;
        const restore = rememberFocus();
        onCleanup(() => restore(card));
        queueMicrotask(() => card?.isConnected && (focusables(card)[0] ?? card).focus({ preventScroll: true }));
        const bar = (<span />) as HTMLElement;
        const progress = (
          <div
            class={classes.progress}
            role="progressbar"
            aria-label={derive(props.busy, (x) => x?.label ?? '')}
            data-indeterminate={derive(props.busy, (x) => x?.progress === undefined)}
          >
            {bar}
          </div>
        ) as HTMLElement;
        watch([props.busy], () => {
          const p = props.busy.value?.progress;
          bar.style.width = p === undefined ? '' : `${Math.round(Math.max(0, Math.min(1, p)) * 100)}%`;
          setAttr(progress, 'aria-valuenow', p === undefined ? false : String(Math.round(p * 100)));
        });
        return (
          <div
            class={classes.busy}
            role="alertdialog"
            aria-modal="true"
            aria-label={derive(props.busy, (x) => x?.label ?? '')}
            onKeyDown={(e: KeyboardEvent) => trapTab(e, card)}
          >
            <div class={classes.busyCard} tabIndex={-1} ref={(el: HTMLElement) => (card = el)}>
              <div class={classes.busyRow}>
                <span style="flex: 1">{derive(props.busy, (x) => x?.label ?? '')}</span>
                <Show when={cancellable}>
                  {() => (
                    <Button size="sm" onClick={() => props.busy.value?.cancel?.()}>
                      {props.cancel}
                    </Button>
                  )}
                </Show>
              </div>
              {progress}
            </div>
          </div>
        );
      }}
    </Show>
  );
}

/** Files dragged over the shell: an overlay saying what a drop does; dropped, they go to `onDrop`. */
function FileDrop(props: { onDrop: (files: File[]) => void; label: string; formats?: string }) {
  const over = new State(false);
  let depth = 0;
  const hasFiles = (e: DragEvent) => !!e.dataTransfer && [...e.dataTransfer.types].includes('Files');
  const target = () => (anchor.parentElement ?? anchor) as HTMLElement;
  const anchor = (<span hidden />) as HTMLElement;
  queueMicrotask(() => {
    const t = target();
    t.addEventListener('dragenter', enter);
    t.addEventListener('dragover', overHandler);
    t.addEventListener('dragleave', leave);
    t.addEventListener('drop', drop);
  });
  const enter = (e: DragEvent) => {
    if (!hasFiles(e)) return;
    depth++;
    over.set(true);
  };
  const overHandler = (e: DragEvent) => {
    if (!hasFiles(e)) return;
    e.preventDefault();
    // A drop opens a copy of the files (the cursor says so), never moves them.
    if (e.dataTransfer) e.dataTransfer.dropEffect = 'copy';
  };
  const leave = () => {
    depth = Math.max(0, depth - 1);
    if (depth === 0) over.set(false);
  };
  const drop = (e: DragEvent) => {
    if (!hasFiles(e)) return;
    depth = 0;
    over.set(false);
    // A view under the drop that took it (prevented its default) handled it: not the shell's too.
    if (e.defaultPrevented) return;
    e.preventDefault();
    const files = [...(e.dataTransfer?.files ?? [])];
    if (files.length) props.onDrop(files);
  };
  onCleanup(() => {
    const t = target();
    t.removeEventListener('dragenter', enter);
    t.removeEventListener('dragover', overHandler);
    t.removeEventListener('dragleave', leave);
    t.removeEventListener('drop', drop);
  });
  return (
    <>
      {anchor}
      <Show when={over}>
        {() => (
          <div class={classes.drop}>
            <div class={classes.dropText}>
              <span>{props.label}</span>
              {props.formats ? <span class={classes.dropFormats}>{props.formats}</span> : null}
            </div>
          </div>
        )}
      </Show>
    </>
  );
}

// ---------------------------------------------------------------------------------------------------------------------
// Compact

/**
 * Compact (shell §3.7): the main area (its region's panes and top kept where they are) over a bottom sheet; the sheet
 * shows the card, else the view picked last: a window's (a window is a card here) or one docked beside the main area;
 * its bar holds those docked views marked `compactBar` and More (the rest).
 */
function Compact(props: { shell: ShellModel; bodies: ViewBodies; main: HTMLElement }) {
  const { shell } = props;
  // The groups docked beside the main area: the region's and its top's views stay in the main area (their bodies are
  // the panes').
  const sideGroups = (l: ShellLayout) => {
    const ix = layoutIndex(l);
    return ix.docked.filter((g) => !ix.regionGroups.has(g.id));
  };
  const all = derive([shell.layout, shell.absent, shell.views], (l) =>
    sideGroups(l).flatMap((g) => shownViews(shell, g)),
  );
  // What the sheet shows (besides a card): the view picked last, from the bar, More, a `show`, a window opened or
  // raised. Not an order of kinds: with a window open, a docked view picked from the bar still shows.
  const topWindow = derive(shell.layout, (l) => {
    const top = stacked(l).at(-1);
    return top ? activeOf(shell, top.group) : undefined;
  });
  const openDocked = derive(shell.layout, (l) =>
    sideGroups(l)
      .map((g) => (g.collapsed ? '' : (activeOf(shell, g) ?? '')))
      .join('|'),
  );
  const picked = new State<ViewRef | undefined>(topWindow.value ?? openDocked.value.split('|').find((r) => r));
  // A window opened or raised, a docked view shown (its group opened, or another of it activated): that is picked.
  effect(topWindow, (ref) => {
    if (ref) picked.set(ref);
  });
  let lastOpen = openDocked.value.split('|');
  effect(openDocked, (key) => {
    const now = key.split('|');
    const changed = now.find((r, i) => r && r !== lastOpen[i]);
    lastOpen = now;
    if (changed) picked.set(changed);
  });
  const shown = derive([shell.layout, shell.card, shell.absent, shell.views, picked], (l, card, _a, _v, ref) => {
    if (card) return { ref: card.ref, kind: 'card' as const };
    if (!ref || !shell.shows(ref)) return undefined;
    const place = findView(l, ref);
    // A pane (in the main area's region or on its top) shows there, not in the sheet.
    if (!place || (!place.window && inRegion(l, place.group.id))) return undefined;
    return { ref, kind: place.window ? ('window' as const) : ('docked' as const) };
  });
  const key = derive(shown, (s) => (s ? `${s.kind}:${s.ref}` : ''));
  const enterSheet = entering();
  let moreButton!: HTMLElement;
  const bar = derive(all, (refs) => ({
    buttons: refs.filter((r) => shell.view(r)?.compactBar),
    more: refs.filter((r) => !shell.view(r)?.compactBar),
  }));
  return (
    <div class={classes.compact}>
      {props.main}
      <div class={classes.sheet} data-open={derive(key, (k) => String(!!k))} data-tip-edge="">
        <div class={classes.sheetContent} hidden={derive(key, (k) => !k)}>
          <Show when={key}>
            {() => {
              const s = shown.value!;
              return enterSheet(
                <KeptFrame
                  shell={shell}
                  bodies={props.bodies}
                  ref={s.ref}
                  placement="sheet"
                  onClose={
                    s.kind === 'card'
                      ? () => shell.close(s.ref)
                      : s.kind === 'window'
                        ? () => {
                            picked.set(undefined);
                            shell.close(s.ref);
                          }
                        : () => picked.set(undefined)
                  }
                />,
              );
            }}
          </Show>
        </div>
        <div class={classes.sheetBar}>
          <div class={classes.sheetTabs} role="tablist" aria-label={shell.labels.compactBar} data-msui-region="">
            <For each={derive(bar, (b) => b.buttons)} by={(r) => r}>
              {(ref) => {
                const def = shell.view(ref)!;
                return (
                  <button
                    type="button"
                    role="tab"
                    class={classes.sheetButton}
                    aria-selected={derive(shown, (s) => String(s?.kind !== 'card' && s?.ref === ref))}
                    onClick={() => {
                      // The picked one again: the sheet closes; another: it shows.
                      if (picked.value === ref && shown.value?.kind !== 'card') return picked.set(undefined);
                      shell.show(viewIdOf(ref), { key: keyOf(ref) });
                      picked.set(ref);
                    }}
                  >
                    <Icon name={def.icon} />
                    <span>{def.title}</span>
                  </button>
                );
              }}
            </For>
          </div>
          <Show when={derive(bar, (b) => b.more.length > 0)}>
            {() => (
              <button
                type="button"
                class={classes.sheetButton}
                aria-haspopup="menu"
                ref={(el: HTMLElement) => (moreButton = el)}
                onClick={() => {
                  const r = moreButton.getBoundingClientRect();
                  openContextMenu(
                    menuRoot(moreButton),
                    { x: r.left, y: r.top },
                    bar.value.more.map((ref) => {
                      const def = shell.view(ref)!;
                      return {
                        id: ref,
                        label: def.title,
                        icon: def.icon,
                        run: () => {
                          shell.show(viewIdOf(ref), { key: keyOf(ref) });
                          picked.set(ref);
                        },
                      };
                    }),
                  );
                }}
              >
                <Icon name="ellipsis" />
                <span>{shell.labels.more}</span>
              </button>
            )}
          </Show>
        </div>
      </div>
    </div>
  );
}

// ---------------------------------------------------------------------------------------------------------------------
// F6

/**
 * The regions F6 moves between (shell §4): the main area, each docked group that shows, the rail, each window, the card,
 * the compact bar. The shell marks its own (`data-msui-region`): not found by a label (translated) or a tag (an app's
 * `nav` in a view), nor in a shell nested in a view.
 */
function cycleRegions(root: HTMLElement, dir: 1 | -1): void {
  const regions = [...root.querySelectorAll<HTMLElement>('[data-msui-region]')].filter(
    (r) =>
      r.closest('[data-msui-shell]') === root &&
      !r.closest('[hidden], [inert], [data-shown="false"]') &&
      r.getClientRects().length > 0,
  );
  if (!regions.length) return;
  const now = deepActiveElement();
  // The innermost region holding the focus (windows are in the main area: a later region in document order).
  let at = -1;
  for (let i = 0; i < regions.length; i++) if (now && regions[i].contains(now)) at = i;
  const next = regions[(at + dir + regions.length) % regions.length];
  const target = next.matches('[data-msui-shell-main]') ? next : (focusables(next)[0] ?? next);
  if (target === next && next.tabIndex < 0 && !next.hasAttribute('tabindex')) next.tabIndex = -1;
  target.focus({ preventScroll: true });
}
