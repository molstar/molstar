/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { derive, type ReadonlyState, State } from '@molstar/core/reactive';
import { Button } from '../components/button/button.js';
import { openContextMenu, openMenuBeside } from '../components/context-menu/context-menu.js';
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
  allGroups,
  TopGroup,
  containsLeaf,
  DefaultGroupSize,
  columnOf,
  dockedGroups,
  findGroup,
  MinGroupSize,
  MinPaneSize,
  findView,
  gridPeers,
  inRegion,
  type GroupNode,
  type LayoutNode,
  mainNode,
  type Rect,
  railGroup,
  rightColumn,
  type ShellLayout,
  stacked,
  type ViewRef,
} from './layout.js';
import { fitRail, type RailEntry, railEntries } from './rail.js';
import {
  type PanePort,
  type RailItem,
  type ShellModel,
  type ShellViewDef,
  type ViewPlacement,
  viewIdOf,
} from './shell-model.js';
import { viewBody, viewContext, ViewFrame } from './view-frame.js';
import { paneMenu, railMenu, viewMenu, type WindowAdjust } from './view-menu.js';

/** What holds the app (shell §3.6, busy): a label, a progress 0–1 (none: indeterminate), a cancel. */
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
  /** Files dropped anywhere on the shell; without it, no drop overlay. */
  readonly onDropFiles?: (files: File[]) => void;
  /** What the drop overlay says. */
  readonly dropLabel?: string;
  /** What can be dropped, a hint under the label ("PDB, mmCIF, sessions"). */
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
  /** The shell's accessible name (its region). */
  readonly label?: string;
  readonly class?: string;
}

/**
 * The shell (shell §1–4): the app's main area with the docked layout around it (the right column of groups, the rail
 * at the edge), windows over the main area, a card in the right column's place, a page over the main area, a modal
 * over everything; compact below its compact width (the views in a bottom sheet). Mode by its own width, the rail's
 * steps by its height.
 */
export function ShellView(props: ShellViewProps) {
  const { shell } = props;
  const owner = requireScope('ShellView');
  const mainArea = (
    <div class={classes.main} data-msui-shell-main="" tabIndex={-1} aria-label="Main area" role="region" />
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
      aria-roledescription="workspace"
      aria-label={props.label ?? 'Workspace'}
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
    if (shell.expanded.value) shell.collapseExpanded();
    else if (shell.maximized.value) shell.maximize(undefined);
    else return;
    e.preventDefault();
  });
  // Fullscreen: the model asks, the view does it on the shell's root and writes back what the browser did (its own Esc
  // too); where the page or an iframe forbids it, the model says why.
  const doc = root.ownerDocument;
  shell.setFullscreenAvailable(
    typeof root.requestFullscreen !== 'function'
      ? 'This browser has no fullscreen'
      : doc.fullscreenEnabled === false
        ? 'Fullscreen is not allowed here (the page embedding this does not allow it)'
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
        {props.busy ? <BusyOverlay busy={props.busy} /> : null}
        {props.onDropFiles ? (
          <FileDrop onDrop={props.onDropFiles} label={props.dropLabel} formats={props.dropFormats} />
        ) : null}
      </>
    ) as Node,
  );
  // The layers over the main area (in it in both modes): windows (desktop only), the page, the status slot.
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
    // registered and layered (a key: not every layout write, a drag's sizes, sweeps).
    effect(
      derive([shell.layout, shell.views, shell.card, shell.page, shell.modal], (l, v, c, p, m) => {
        let key = `${c?.ref}|${p?.ref}|${m?.ref}`;
        for (const g of allGroups(l)) key += `|${g.views.join(',')}`;
        for (const [id, def] of v) key += `|${id}:${defId(def)}`;
        return key;
      }),
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

/** The instance key of a ref (`id#key`), if it has one. */
function keyOf(ref: ViewRef): string | undefined {
  const i = ref.indexOf('#');
  return i < 0 ? undefined : ref.slice(i + 1);
}

/** A number for each view definition (a key that changes when a view is registered anew under its id). */
const defIds = new WeakMap<ShellViewDef, number>();
let nextDefId = 1;
function defId(def: ShellViewDef | undefined): number {
  if (!def) return 0;
  let id = defIds.get(def);
  if (id === undefined) defIds.set(def, (id = nextDefId++));
  return id;
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
    covered: (l, n) => !!shell.card.value && n.kind === 'group' && l.rail && isInRightColumn(l, n.id),
    // Expanded (shell §3.9): the main area alone; every group away, paused.
    away: (_, n) => shell.expanded.value && n.kind === 'group',
    sources: [shell.card, shell.expanded],
  });
  dock.append(tree, cardSlot);
  // The card: over the right column (the rail beside it), or at the main area's right edge.
  const cardWidth = derive(shell.layout, (l) => {
    const col = rightColumn(l);
    return col && l.rail && (col.kind === 'group' || col.kind === 'split') ? (col.size ?? 340) : 340;
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
  // The region's node shown alone: the maximized one, when it is the region's (a window's comes out over it instead).
  const alone = (l: ShellLayout): string | undefined => {
    const m = shell.maximized.value;
    return m && (m === 'content' || inRegion(l, m)) ? m : undefined;
  };
  const tree = DockTree({
    shell,
    bodies,
    leaf: 'content',
    leafEl: props.content,
    rootOf: (l) => mainNode(l).root ?? { kind: 'content', id: 'content' },
    // Maximized (shell §2.1): only that pane and what holds it.
    away: (l, n) => {
      const m = alone(l);
      return !!m && !!mainNode(l).root && !holds(n, m);
    },
    sources: [shell.maximized],
    motion: false,
  });
  tree.classList.add(classes.region);
  // Several panes: the active one is marked (a shared toolbar follows it, shell §2.1).
  // Several panes shown (a maximized one shows alone: no mark of being the active one).
  const several = derive([shell.layout, shell.maximized], (l, m) => !!mainNode(l).root && !m);
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
  const marked = derive([shell.activePane, shell.layout], (p) => shell.paneNode(p));
  // The shell's own mark of the active pane (several panes); what a pane shows on its status is the pane's (`ctx.pane`).
  watch([marked, several], () => {
    for (const el of tree.querySelectorAll('[data-pane]')) {
      setAttr(el, 'data-active', props.mark && several.value && el.getAttribute('data-pane') === marked.value);
    }
  });
  return tree;
}

/** Whether a node is a pane (`content` or a group's id) or holds it. */
function holds(n: LayoutNode, pane: string): boolean {
  if (n.id === pane) return true;
  return n.kind === 'split' && n.children.some((c) => holds(c, pane));
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
}

/**
 * A docked tree (shell §3): groups and splits around the level's leaf. Built only when its structure changes
 * (`shapeOf`); everything else is written in place (`applyState`): which nodes show and which are folded, which takes
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
  const applyState = (l: ShellLayout) => {
    const walk = (n: LayoutNode, flex: boolean, at: 'start' | 'end' | undefined) => {
      const el = nodes.get(n.id);
      const covered = o.covered?.(l, n) ?? false;
      const shows = nodeShows(shell, l, n) && !(o.away?.(l, n) ?? false);
      if (el) {
        showNode(el, n, shows);
        setAttr(el, 'data-folded', n.kind === 'group' && !!n.collapsed && railGroup(l)?.id !== n.id);
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
      if (n.kind === 'group') groups.get(n.id)?.scope.setActive(shows && !covered);
      if (n.kind !== 'split') return;
      const f = flexIndex(shell, l, n, o.leaf, o.away);
      n.children.forEach((c, i) => walk(c, i === f, i === f ? undefined : i > f ? 'start' : 'end'));
    };
    walk(o.rootOf(l), true, undefined);
  };
  effect(
    derive([shell.layout, shell.absent, shell.views, ...(o.sources ?? [])], (l, a, v, ...rest) => ({ l, a, v, rest })),
    ({ l }) => {
      applyState(l);
      settled = true;
    },
  );
  let buildScope: ViewScope | undefined;
  const shape = derive(shell.layout, (l) => shapeOf(l, o.rootOf(l)));
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

/** Whether a node shows: a group with a view to show (a rail group folded does not: its rail stays), a split with one. */
function nodeShows(shell: ShellModel, l: ShellLayout, n: LayoutNode): boolean {
  if (n.kind === 'main' || n.kind === 'content') return true;
  if (n.kind === 'group') {
    if (shownViews(shell, n).length === 0) return false;
    if (n.collapsed && railGroup(l)?.id === n.id) return false;
    if (!l.rail && isInRightColumn(l, n.id)) return false;
    // A folded group shows its header only beside an open group of its column: with every group of the column folded,
    // the column closes to the rail (no header left alone over an empty column).
    if (n.collapsed) {
      const column = columnOf(l, n.id);
      const open = (c: LayoutNode): boolean =>
        c.kind === 'group'
          ? c.id !== n.id && !c.collapsed && shownViews(shell, c).length > 0
          : c.kind === 'split' && c.children.some(open);
      return !!column && column.children.some(open);
    }
    return true;
  }
  return n.children.some((c) => nodeShows(shell, l, c));
}

function isInRightColumn(l: ShellLayout, id: string): boolean {
  const col = rightColumn(l);
  if (!col) return false;
  if (col.kind === 'group') return col.id === id;
  const walk = (n: LayoutNode): boolean =>
    n.kind === 'group' ? n.id === id : n.kind === 'split' && n.children.some(walk);
  return walk(col);
}

/**
 * The docked tree's structure: its nodes, their order and direction, which group is the rail's. Equal structures build
 * equal trees; the rest (shown, folded, which takes the rest, sizes) is written in place.
 */
function shapeOf(l: ShellLayout, root: LayoutNode): string {
  const rail = railGroup(l)?.id;
  const walk = (n: LayoutNode): string => {
    if (n.kind === 'main' || n.kind === 'content') return n.kind;
    if (n.kind === 'group') return `g:${n.id}${n.id === rail ? '@' : ''}`;
    return `${n.dir}:${n.id}(${n.children.map(walk).join(',')})`;
  };
  return walk(root);
}

/**
 * The child of a split that takes the rest: the one holding the main area, else the first that shows open with no
 * size, else the first that shows open (a folded group gives its height to the other).
 */
function flexIndex(
  shell: ShellModel,
  l: ShellLayout,
  n: LayoutNode & { kind: 'split' },
  leaf: 'main' | 'content',
  away?: (l: ShellLayout, n: LayoutNode) => boolean,
): number {
  const open = (c: LayoutNode) =>
    nodeShows(shell, l, c) && !(c.kind === 'group' && c.collapsed) && !(away?.(l, c) ?? false);
  // The leaf takes the rest, unless it is away (another pane maximized): then the pane that shows does.
  let flex = n.children.findIndex((c) => containsLeaf(c, leaf) && !(away?.(l, c) ?? false));
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

function sizeOf(n: LayoutNode, id: string): number | undefined {
  if (n.kind === 'main') return n.root ? sizeOf(n.root, id) : undefined;
  if (n.kind === 'content') return undefined;
  if (n.id === id) return n.size;
  if (n.kind === 'split')
    for (const c of n.children) {
      const s = sizeOf(c, id);
      if (s !== undefined) return s;
    }
  return undefined;
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
      aria-label="Resize"
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
  const value = derive(shell.layout, (l) => sizeOf(l.root, id));
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
  const state = derive([shell.layout, shell.absent, shell.views], (l) => {
    const g = findGroup(l, id)?.group;
    const rail = railGroup(l)?.id === id;
    const views = shownViews(shell, g);
    return { active: activeOf(shell, g), rail, views, key: `${activeOf(shell, g)}|${rail}|${views.join(',')}` };
  });
  const key = derive(state, (s) => s.key);
  const open = derive(shell.layout, (l) => !findGroup(l, id)?.group.collapsed);
  const enter = entering();
  return (
    // A row's tooltip card goes outside the group, toward the roomier side (the main area).
    <div class={classes.group} data-group={id} data-tip-edge="">
      <Show when={key}>
        {() => {
          const s = state.value;
          if (!s.active) return null;
          const ref = s.active;
          const tabs = !s.rail && s.views.length > 1 ? <Tabs shell={shell} views={s.views} active={ref} /> : undefined;
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
 * A group in the main area's region (shell §2.1): a pane. A stage alone has no header (a canvas has no chrome above it;
 * its pane's menu is the canvas's right-click, the app's); several views show their tabs; a panel shows its header. Its
 * ✕ closes its views (they leave the region; shown again they go home).
 */
function PaneGroupView(props: { shell: ShellModel; id: string; bodies: ViewBodies }) {
  const { shell, id, bodies } = props;
  const state = derive([shell.layout, shell.absent, shell.views], (l) => {
    const g = findGroup(l, id)?.group;
    const views = shownViews(shell, g);
    const active = activeOf(shell, g);
    return { active, views, key: `${active}|${views.join(',')}` };
  });
  const key = derive(state, (s) => s.key);
  const enter = entering();
  return (
    <div class={classes.group} data-group={id} data-tip-edge="">
      <Show when={key}>
        {() => {
          const s = state.value;
          if (!s.active) return null;
          const ref = s.active;
          const def = shell.view(ref);
          const tabs = s.views.length > 1 ? <Tabs shell={shell} views={s.views} active={ref} /> : undefined;
          let header!: HTMLElement;
          return enter(
            <KeptFrame
              shell={shell}
              bodies={bodies}
              ref={ref}
              placement="stage"
              head={def?.stage && !tabs ? 'none' : 'header'}
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
  const state = derive([shell.layout, shell.absent, shell.views], (l) => {
    const g = mainNode(l).top;
    const views = shownViews(shell, g);
    const active = activeOf(shell, g);
    return { active, views, key: `${active}|${views.join(',')}` };
  });
  const key = derive(state, (s) => s.key);
  return (
    <div class={classes.topOverlay} data-group={TopGroup}>
      <Show when={key}>
        {() => {
          const s = state.value;
          if (!s.active) return null;
          const tabs = s.views.length > 1 ? <Tabs shell={shell} views={s.views} active={s.active} /> : undefined;
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

/** A group's labelled tabs (a group not at the rail, a window of several views): a tab list, the arrows between. */
function Tabs(props: { shell: ShellModel; views: readonly ViewRef[]; active: ViewRef }) {
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
              openContextMenu(
                menuRoot(e.currentTarget as Element),
                { x: e.clientX, y: e.clientY },
                viewMenu(shell, ref),
              );
            }}
            onKeyDown={(e: KeyboardEvent) => {
              if ((e.key === 'F10' && e.shiftKey) || e.key === 'ContextMenu') {
                e.preventDefault();
                const r = (e.currentTarget as Element).getBoundingClientRect();
                openContextMenu(menuRoot(e.currentTarget as Element), { x: r.left, y: r.bottom }, viewMenu(shell, ref));
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
  const height = new State(0);
  // The rail's buttons, a new list only when they change (not at each layout write: a drag's sizes, a window's move).
  const entries = derive(
    [shell.layout, shell.absent, shell.views, shell.railItems],
    (l, absent, _, items) => {
      const g = railGroup(l);
      const views = shownViews(shell, g).map((ref) => ({ ref, def: shell.view(ref)! }));
      const tools = items.filter((i) => !absent.has(`item:${i.id}`));
      return railEntries(views, tools);
    },
    (a, b) =>
      a.length === b.length &&
      a.every((e, i) => e.id === b[i].id && e.view?.def === b[i].view?.def && e.item === b[i].item),
  );
  const fit = derive(
    [entries, height],
    (e, h) => fitRail(e, h),
    (a, b) => a.step === b.step && a.hidden.size === b.hidden.size && [...a.hidden].every((id) => b.hidden.has(id)),
  );
  const active = derive([shell.layout, shell.absent, shell.views], (l) => {
    const g = railGroup(l);
    return g && !g.collapsed ? activeOf(shell, g) : undefined;
  });
  const rail = (<nav class={classes.rail} aria-label="Views and tools" />) as HTMLElement;
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
    const list = derive(entries, section(name), (a, b) => a.length === b.length && a.every((x, i) => x.id === b[i].id));
    // A section whose buttons are all under More draws nothing (no hairline of its own).
    const any = derive([list, fit], (l, f) => l.some((e) => !f.hidden.has(e.id)));
    return (
      <Show when={any}>
        {() => (
          <div
            class={classes.railSection}
            role={tablist ? 'tablist' : 'group'}
            aria-orientation={tablist ? 'vertical' : undefined}
            aria-label={name}
          >
            <For each={list} by={(e) => e.id}>
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
                  aria-label="More"
                  data-tip="More"
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
                aria-label="Hide the views"
                data-tip="Hide the views"
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
            title="Show the views"
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
  const view = derive([shell.layout, shell.absent, shell.views, shell.maximized], (l, _a, _v, maximized) => {
    const g = findGroup(l, id)?.group;
    const views = shownViews(shell, g);
    const active = activeOf(shell, g);
    // Its view maximized (shell §2.1): out over the main area (`MaximizedHost`), its body back here after.
    const out = !!active && active === maximized;
    return { active: out ? undefined : active, views, key: `${active}|${views.join(',')}|${out}` };
  });
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
  watch([key], () => {
    const active = view.value.active;
    setAttr(el, 'hidden', !active);
    setAttr(el, 'aria-label', active ? shell.paneTitle(active).value : 'Window');
  });
  on(el, 'pointerdown', () => shell.raise(id), { capture: true });
  // A window's view is a pane too (a detached stage): pressed or focused, it is the active one.
  const activate = () => {
    const active = view.value.active;
    if (active) shell.setActivePane(active);
  };
  on(el, 'pointerdown', activate, { capture: true });
  on(el, 'focusin', activate);
  on(el, 'focusin', () => shell.raise(id));
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
                tabs={v.views.length > 1 ? <Tabs shell={shell} views={v.views} active={ref} /> : undefined}
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
                      const size = def?.size ?? { width: 360, height: 320 };
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
        <div class={classes.maximizedPane} data-maximized-pane={ref}>
          <KeptFrame
            shell={shell}
            bodies={bodies}
            ref={ref}
            placement="stage"
            head={shell.view(ref)?.stage ? 'none' : 'header'}
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
            close: () => shell.closeCard(),
            outside: !open.def.persistent,
            // Escape anywhere is the card's, but not with the focus in a window or a dialog: that one's (shell §3.6).
            ownsKey: (e) => !(eventTarget(e) as Element | null)?.closest?.('[data-window], [role="dialog"]'),
            refocus: () => open.anchor as HTMLElement | undefined,
          });
          el = (
            <div class={classes.card} data-card={open.ref} data-tip-edge="">
              <KeptFrame
                shell={shell}
                bodies={props.bodies}
                ref={open.ref}
                placement="card"
                onClose={() => shell.closeCard()}
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
        dismissable({ layer: () => el, close: () => shell.closePage(), outside: false, ownsKey: () => true });
        el = (
          <div class={classes.page} role="region" aria-label={open.def.title} tabIndex={-1}>
            <KeptFrame
              shell={shell}
              bodies={props.bodies}
              ref={open.ref}
              placement="page"
              onClose={() => shell.closePage()}
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
        dismissable({ layer: () => el, close: () => shell.closeModal(), outside: false, ownsKey: () => true });
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
                onClose={() => shell.closeModal()}
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
 * the focus goes to it (its Cancel, else the card) and stays in it, and goes back when it ends.
 */
function BusyOverlay(props: { busy: ReadonlyState<ShellBusy | undefined> }) {
  // A new task replacing one (another label) keeps the overlay: shown while there is one.
  const some = derive(props.busy, (b) => !!b);
  return (
    <Show when={some} reveal="fade">
      {() => {
        const b = props.busy.value!;
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
                {b?.cancel ? (
                  <Button size="sm" onClick={() => props.busy.value?.cancel?.()}>
                    Cancel
                  </Button>
                ) : null}
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
function FileDrop(props: { onDrop: (files: File[]) => void; label?: string; formats?: string }) {
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
    e.preventDefault();
    depth = 0;
    over.set(false);
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
              <span>{props.label ?? 'Drop files to open them'}</span>
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
 * Compact (shell §3.7): the main area over a bottom sheet; the sheet shows the card, else the top window's view (a window
 * opens as a card here), else the open docked view; its bar holds the views marked `compactBar` and More (the rest).
 */
function Compact(props: { shell: ShellModel; bodies: ViewBodies; main: HTMLElement }) {
  const { shell } = props;
  const all = derive([shell.layout, shell.absent, shell.views], (l) =>
    dockedGroups(l).flatMap((g) => shownViews(shell, g)),
  );
  // What the sheet shows (besides a card): the view picked last, from the bar, More, a `show`, a window opened or
  // raised. Not an order of kinds: with a window open, a docked view picked from the bar still shows.
  const topWindow = derive(shell.layout, (l) => {
    const top = stacked(l).at(-1);
    return top ? activeOf(shell, top.group) : undefined;
  });
  const openDocked = derive(shell.layout, (l) =>
    dockedGroups(l)
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
    if (!place) return undefined;
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
                      ? () => shell.closeCard()
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
          <div class={classes.sheetTabs} role="tablist" aria-label="Views">
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
                <span>More</span>
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

/** The regions F6 moves between (shell §4): the main area, each docked group that shows, the rail, each window. */
function cycleRegions(root: HTMLElement, dir: 1 | -1): void {
  const regions = [
    ...root.querySelectorAll<HTMLElement>(
      '[data-msui-shell-main], [data-group], nav, [data-window], [data-card], [role="tablist"][aria-label="Views"]',
    ),
  ].filter((r) => !r.closest('[hidden], [inert], [data-shown="false"]') && r.getClientRects().length > 0);
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
