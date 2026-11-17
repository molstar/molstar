/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { derive, type ReadonlyState, State } from '@molstar/core/reactive';
import { Button } from '../components/button/button.js';
import { openContextMenu } from '../components/context-menu/context-menu.js';
import { Icon } from '../components/icon/icon.js';
import { installTooltips } from '../components/tooltip/tooltip.js';
import { dismissable } from '../components/popover/dismiss.js';
import {
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
  DefaultGroupSize,
  columnOf,
  findGroup,
  MinGroupSize,
  findView,
  type GroupNode,
  type LayoutNode,
  type Rect,
  railGroup,
  rightColumn,
  type ShellLayout,
  stacked,
  type ViewRef,
} from './layout.js';
import { fitRail, type RailEntry, railEntries } from './rail.js';
import type { ShellModel, ViewPlacement } from './shell-model.js';
import { viewBody, viewContext, ViewFrame } from './view-frame.js';
import { railMenu, viewMenu, type WindowAdjust } from './view-menu.js';

/** What holds the app (shell §3.6, busy): a label, a progress 0–1 (none: indeterminate), a cancel. */
export interface ShellBusy {
  readonly label: string;
  readonly progress?: number;
  cancel?(): void;
}

export interface ShellViewProps {
  readonly shell: ShellModel;
  /** The app's main view (a canvas, a document): rendered once, kept across modes. */
  readonly main: () => Child;
  /** The status slot (the main area's bottom left): the app's. */
  readonly status?: () => Child;
  /** While set, a scrim and its card hold the app. */
  readonly busy?: ReadonlyState<ShellBusy | undefined>;
  /** Files dropped anywhere on the shell; without it, no drop overlay. */
  readonly onDropFiles?: (files: File[]) => void;
  /** What the drop overlay says. */
  readonly dropLabel?: string;
  /**
   * The rail can be hidden (default on): a Hide button at its foot, and the Mol* star in the main area's corner that
   * shows it again.
   */
  readonly railToggle?: boolean;
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
  const bodies = new ViewBodies(shell, owner, { menuButton: !!props.viewMenu, railToggle: props.railToggle !== false });
  onCleanup(() => bodies.dispose());

  // The app's main view: made once, moved between the desktop and compact layouts (its state, a canvas, kept).
  const mainContent = (<div class={classes.mainContent}>{props.main()}</div>) as HTMLElement;
  const mainArea = (
    <div class={classes.main} data-msui-shell-main="" tabIndex={-1} aria-label="Main area" role="region" />
  ) as HTMLElement;
  mainArea.append(mainContent);
  observeResize(mainArea, () => shell.setArea(mainArea.clientWidth, mainArea.clientHeight));

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
  on<KeyboardEvent>(root, 'keydown', (e) => {
    if (e.key === 'F6' && !e.defaultPrevented) {
      e.preventDefault();
      cycleRegions(root, e.shiftKey ? -1 : 1);
    }
  });

  const compact = derive(shell.mode, (m) => m === 'compact');
  root.append(
    (
      <>
        <Show when={compact} fallback={() => <Desktop shell={shell} bodies={bodies} main={mainArea} props={props} />}>
          {() => <Compact shell={shell} bodies={bodies} main={mainArea} />}
        </Show>
        <ModalHost shell={shell} bodies={bodies} />
        {props.busy ? <BusyOverlay busy={props.busy} /> : null}
        {props.onDropFiles ? <FileDrop onDrop={props.onDropFiles} label={props.dropLabel} /> : null}
      </>
    ) as Node,
  );
  // The layers over the main area (in it in both modes): windows (desktop only), the page, the status slot.
  mainArea.append(
    (
      <>
        <Show when={derive(compact, (c) => !c)}>{() => <WindowsLayer shell={shell} bodies={bodies} />}</Show>
        <PageHost shell={shell} bodies={bodies} />
        {props.status ? <div class={classes.status}>{props.status()}</div> : null}
        {props.railToggle === false ? null : <RailToggle shell={shell} compact={compact} />}
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
  readonly #bodies = new Map<ViewRef, { el: HTMLElement; scope: ViewScope; owner: object | undefined }>();
  readonly #placements = new Map<ViewRef, ViewPlacement>();

  constructor(
    readonly shell: ShellModel,
    readonly owner: ViewScope,
    /** What the shell's frames show (its props): the ⋯ button, the rail's Hide. */
    readonly ui: { readonly menuButton: boolean; readonly railToggle: boolean },
  ) {
    // A view that leaves the layout and every layer: its body goes.
    effect(
      derive([shell.layout, shell.views, shell.card, shell.page, shell.modal], (l) => l),
      () => this.#sweep(),
    );
  }

  /** The body of `ref`, now shown by `owner` (a frame): active, `placement` where it shows. */
  take(ref: ViewRef, owner: object, placement: ViewPlacement): HTMLElement | undefined {
    const def = this.shell.view(ref);
    if (!def) return undefined;
    this.#placements.set(ref, placement);
    let b = this.#bodies.get(ref);
    if (!b) {
      const scope = this.owner.child();
      const ctx = viewContext(this.shell, ref, () => this.#placements.get(ref) ?? placement);
      const el = runInScope(scope, () => viewBody(def, ctx));
      b = { el, scope, owner: undefined };
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

  /** Whether `ref` left the layout and every layer (or is no longer registered). */
  #gone(ref: ViewRef): boolean {
    const layered = [this.shell.card.value, this.shell.page.value, this.shell.modal.value].some((x) => x?.ref === ref);
    return !layered && !(findView(this.shell.layout.value, ref) && this.shell.view(ref));
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
  head?: 'rail' | 'header';
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
  const owner = requireScope('Desktop');
  const dock = (<div class={classes.dock} />) as HTMLElement;
  // The card's place in the dock, kept apart from the tree (which is replaced at each layout change).
  const cardSlot = (<div class={classes.cardSlot} />) as HTMLElement;
  // Each group's element, kept by id while the tree around it changes (its view stays mounted).
  const groups = new Map<string, { el: HTMLElement; scope: ViewScope }>();
  // The sized nodes' elements by id: a size changes in place (a splitter's drag goes on, the main area stays put).
  const sized = new Map<string, HTMLElement>();
  // Every node's element by id: which show and which are folded are written in place (folding the right group or
  // hiding the rail does not rebuild the tree: the main area stays put, with the focus in it).
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
    if (shows === was || !settled || n.kind !== 'group') return setAttr(el, 'data-shown', shows ? false : 'false');
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
  const applyState = (l: ShellLayout) => {
    const walk = (n: LayoutNode) => {
      const el = nodes.get(n.id);
      const shows = nodeShows(shell, l, n);
      if (el) {
        showNode(el, n, shows);
        setAttr(el, 'data-folded', n.kind === 'group' && !!n.collapsed && railGroup(l)?.id !== n.id);
      }
      // Hidden UI is paused (shell §3.3): a folded rail group, the column hidden with the rail.
      if (n.kind === 'group') groups.get(n.id)?.scope.setActive(shows);
      if (n.kind === 'split') for (const c of n.children) walk(c);
    };
    walk(l.root);
  };
  effect(
    derive([shell.layout, shell.absent, shell.views], (l) => l),
    (l) => {
      applyState(l);
      settled = true;
    },
  );
  let buildScope: ViewScope | undefined;
  // Rebuilt only when the tree's shape changes (not on a size, nor on a window's move or raise): rebuilding moves the
  // main area, which drops the focus in it (a window's) and ends a splitter's drag.
  const shape = derive([shell.layout, shell.absent, shell.views], (l) => shapeOf(shell, l));
  effect(shape, () => {
    const l = shell.layout.value;
    buildScope?.dispose();
    buildScope = owner.child();
    const used = new Set<string>();
    sized.clear();
    nodes.clear();
    const tree = runInScope(buildScope, () =>
      buildNode(shell, l, l.root, undefined, { groups, used, owner, bodies, main: props.main, sized, nodes }),
    );
    dock.replaceChildren(tree, cardSlot);
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
  effect(
    derive(shell.layout, (l) => l.root),
    (root) => {
      const walk = (n: LayoutNode) => {
        if (n.kind === 'split') for (const c of n.children) walk(c);
        if (n.kind === 'main') return;
        sized.get(n.id)?.style.setProperty('--msui-shell-size', `${n.size ?? DefaultGroupSize}px`);
      };
      walk(root);
    },
  );
  onCleanup(() => {
    for (const g of groups.values()) g.scope.dispose();
  });
  // The card: over the right column (the rail beside it), or at the main area's right edge.
  const cardWidth = derive(shell.layout, (l) => {
    const col = rightColumn(l);
    return col && l.rail && col.kind !== 'main' ? (col.size ?? 340) : 340;
  });
  watch([cardWidth], () => dock.style.setProperty('--msui-shell-card', `${cardWidth.value}px`));
  dock.style.position = 'relative';
  return (
    <>
      {dock}
      <CardHost shell={shell} bodies={bodies} into={cardSlot} />
      <Show when={derive([shell.layout, shell.absent, shell.views], (l) => !!railGroup(l) && l.rail)} reveal="side">
        {() => <Rail shell={shell} hideable={bodies.ui.railToggle} />}
      </Show>
    </>
  );
}

/** Whether a node shows: a group with a view to show (a rail group folded does not: its rail stays), a split with one. */
function nodeShows(shell: ShellModel, l: ShellLayout, n: LayoutNode): boolean {
  if (n.kind === 'main') return true;
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
 * What the docked tree's elements follow, but its sizes: the nodes, their order and direction, which show, which are
 * folded, which take the rest. Equal shapes build equal trees.
 */
function shapeOf(shell: ShellModel, l: ShellLayout): string {
  const rail = railGroup(l)?.id;
  // Which show and which are folded are written in place (`applyState`), except where they choose the node that takes
  // the rest: in a split without the main area (a folded group gives its height to the other).
  const walk = (n: LayoutNode, flexByOpen: boolean): string => {
    const state = flexByOpen
      ? `${nodeShows(shell, l, n) ? '' : '-'}${n.kind === 'group' && n.collapsed ? '~' : ''}`
      : '';
    if (n.kind === 'main') return 'main';
    if (n.kind === 'group') return `g:${n.id}${state}${n.size === undefined ? '*' : ''}${n.id === rail ? '@' : ''}`;
    const byOpen = !n.children.some(containsMainNode);
    return `${n.dir}:${n.id}${state}${n.size === undefined ? '*' : ''}(${n.children.map((c) => walk(c, byOpen)).join(',')})`;
  };
  return walk(l.root, false);
}

interface BuildContext {
  readonly groups: Map<string, { el: HTMLElement; scope: ViewScope }>;
  readonly used: Set<string>;
  readonly owner: ViewScope;
  readonly bodies: ViewBodies;
  readonly main: HTMLElement;
  readonly sized: Map<string, HTMLElement>;
  readonly nodes: Map<string, HTMLElement>;
}

function buildNode(
  shell: ShellModel,
  l: ShellLayout,
  n: LayoutNode,
  parent: { dir: 'row' | 'column'; flex: boolean; at: 'start' | 'end' | undefined } | undefined,
  ctx: BuildContext,
): HTMLElement {
  const { groups, used, owner, bodies, main } = ctx;
  const wrap = (<div class={classes.node} data-node={n.id} />) as HTMLElement;
  if (!parent || parent.flex) wrap.setAttribute('data-flex', '');
  else if (n.kind !== 'main') {
    wrap.style.setProperty('--msui-shell-size', `${n.size ?? DefaultGroupSize}px`);
    ctx.sized.set(n.id, wrap);
  }
  ctx.nodes.set(n.id, wrap);
  if (n.kind === 'main') {
    wrap.append(main);
    return wrap;
  }
  if (n.kind === 'group') {
    used.add(n.id);
    let g = groups.get(n.id);
    if (!g) {
      const scope = owner.child();
      const el = runInScope(scope, () => (<GroupView shell={shell} id={n.id} bodies={bodies} />) as HTMLElement);
      g = { el, scope };
      groups.set(n.id, g);
    }
    wrap.append(g.el);
  } else {
    const split = (<div class={classes.split} data-dir={n.dir} />) as HTMLElement;
    // The child that takes the rest: the one holding the main area, else the first that shows open (a folded group
    // above a open one gives it its height).
    const open = (c: LayoutNode) => nodeShows(shell, l, c) && !(c.kind === 'group' && c.collapsed);
    let flex = n.children.findIndex((c) => c.kind === 'main' || (c.kind === 'split' && containsMainNode(c)));
    if (flex < 0) flex = n.children.findIndex((c) => (c as GroupNode).size === undefined && open(c));
    if (flex < 0) flex = n.children.findIndex(open);
    n.children.forEach((c, i) => {
      const at = i === flex ? undefined : i > flex ? 'start' : 'end';
      split.append(buildNode(shell, l, c, { dir: n.dir, flex: i === flex, at }, ctx));
    });
    wrap.append(split);
  }
  if (parent && !parent.flex && parent.at) {
    wrap.append(<Splitter shell={shell} id={n.id} dir={parent.dir} at={parent.at} wrap={wrap} />);
  }
  return wrap;
}

function containsMainNode(n: LayoutNode): boolean {
  return n.kind === 'main' || (n.kind === 'split' && n.children.some(containsMainNode));
}

function sizeOf(n: LayoutNode, id: string): number | undefined {
  if (n.kind === 'main') return undefined;
  if (n.id === id) return n.size;
  if (n.kind === 'split')
    for (const c of n.children) {
      const s = sizeOf(c, id);
      if (s !== undefined) return s;
    }
  return undefined;
}

/** A splitter on a sized node's edge toward the node that takes the rest: drags or arrows resize it. */
function Splitter(props: {
  shell: ShellModel;
  id: string;
  dir: 'row' | 'column';
  at: 'start' | 'end';
  wrap: HTMLElement;
}) {
  const { shell, id, dir, at } = props;
  const sign = at === 'start' ? -1 : 1;
  const size = () => {
    const r = props.wrap.getBoundingClientRect();
    return dir === 'row' ? r.width : r.height;
  };
  const el = (
    <div
      class={classes.splitter}
      data-dir={dir}
      data-at={at}
      role="separator"
      tabIndex={0}
      aria-orientation={dir === 'row' ? 'vertical' : 'horizontal'}
      aria-label="Resize"
      aria-valuemin={String(MinGroupSize)}
      onKeyDown={(e: KeyboardEvent) => {
        const step = e.shiftKey ? 64 : 16;
        const keys = dir === 'row' ? ['ArrowLeft', 'ArrowRight'] : ['ArrowUp', 'ArrowDown'];
        const i = keys.indexOf(e.key);
        if (i < 0) return;
        e.preventDefault();
        shell.resize(id, size() + (i === 0 ? -1 : 1) * step * sign);
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
  drag(el, {
    start: (e) => {
      e.preventDefault();
      el.setAttribute('data-dragging', '');
      return { from: dir === 'row' ? e.clientX : e.clientY, size: size() };
    },
    move: (e, s) => next.offer(s.size + ((dir === 'row' ? e.clientX : e.clientY) - s.from) * sign),
    end: () => el.removeAttribute('data-dragging'),
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
              menu={() => viewMenu(shell, ref, { anchorRect: () => anchorIn(header, shell) })}
            />,
          );
        }}
      </Show>
    </div>
  );
}

/** A rectangle of `el` relative to the main area (where a window opens beside it). */
function anchorIn(el: Element | undefined, shell: ShellModel): Rect | undefined {
  const main = el?.ownerDocument.querySelector('[data-msui-shell-main]');
  if (!el || !main) return undefined;
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
            {def.heading && typeof def.heading === 'string' ? def.heading : def.title}
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
  const items = derive(shell.railItems, (all) => all);
  const entries = derive([shell.layout, shell.absent, shell.views, items], (l) => {
    const g = railGroup(l);
    const views = shownViews(shell, g).map((ref) => ({ ref, def: shell.view(ref)! }));
    const tools = shell.railItems.value.filter((i) => !i.visible || i.visible.value);
    return railEntries(views, tools);
  });
  const fit = derive([entries, height], (e, h) => fitRail(e, h));
  const active = derive(shell.layout, (l) => {
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
        onClick={(ev: MouseEvent) => item.run(ev.currentTarget as Element)}
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
                    const r = moreButton.getBoundingClientRect();
                    openContextMenu(
                      menuRoot(rail),
                      { x: r.left, y: r.top },
                      overflow.value.map((e) => ({
                        id: e.id,
                        label: e.view ? e.view.def.title : e.item!.title,
                        icon: e.view ? e.view.def.icon : e.item!.icon,
                        run: () => (e.view ? shell.select(e.view.ref) : e.item!.run(moreButton)),
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
        </div>
      </>
    ) as Node,
  );
  return rail;
}

/** The rail's toggle (the Mol* star): in the main area's corner, it shows the rail and the docked column, or hides them. */
function RailToggle(props: { shell: ShellModel; compact: ReadonlyState<boolean> }) {
  const { shell } = props;
  const shown = derive([shell.railShown, props.compact], (r, c) => !c && !r);
  return (
    <Show when={shown}>
      {() => (
        <div class={classes.railToggle}>
          <Button
            variant="ghost"
            icon={<Icon name="molstar" />}
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
  const ids = derive(
    shell.layout,
    (l) => l.windows.map((w) => w.group.id),
    (a, b) => a.length === b.length && a.every((x, i) => x === b[i]),
  );
  // Windows there at the first render stay still; one opened later floats in.
  let ready = false;
  queueMicrotask(() => (ready = true));
  return (
    <div class={classes.windows}>
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
  const win = derive(shell.layout, (l) => l.windows.find((w) => w.group.id === id));
  const view = derive([shell.layout, shell.absent, shell.views], (l) => {
    const g = findGroup(l, id)?.group;
    const views = shownViews(shell, g);
    return { active: activeOf(shell, g), views, key: `${activeOf(shell, g)}|${views.join(',')}` };
  });
  const key = derive(view, (v) => v.key);
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
    setAttr(el, 'aria-label', shell.view(view.value.active ?? '')?.title ?? 'Window');
  });
  on(el, 'pointerdown', () => shell.raise(id), { capture: true });
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
    refocus: () => el.ownerDocument.querySelector<HTMLElement>('[data-msui-shell-main]') ?? undefined,
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
        if (how === 'up' && e) shell.setWindowRect(id, resized(s, e, edge));
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
                if (how === 'up' && e)
                  shell.setWindowRect(id, { ...s.r, x: s.r.x + e.clientX - s.x, y: s.r.y + e.clientY - s.y });
              },
            });
            return enterView(frame);
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

// ---------------------------------------------------------------------------------------------------------------------
// Layers: the card, the page, the modal

/**
 * The card (shell §3.6): one at a time, in the right column's place (over its groups, the rail beside it), full height; it
 * closes on its ✕, Escape, or a press outside it, the main area included (not a persistent card).
 */
function CardHost(props: { shell: ShellModel; bodies: ViewBodies; into: HTMLElement }) {
  const { shell } = props;
  const card = derive(shell.card, (c) => c);
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
    <Show when={derive(shell.modal, (m) => m)} reveal="fade">
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
        // The rest of the shell is inert while it is open.
        const siblings = () => [...(scrim.parentElement?.children ?? [])].filter((c) => c !== scrim);
        queueMicrotask(() => {
          for (const s of siblings()) s.setAttribute('inert', '');
          focusInitial(el, true);
        });
        onCleanup(() => {
          for (const s of siblings()) s.removeAttribute('inert');
        });
        return scrim;
      }}
    </Show>
  );
}

// ---------------------------------------------------------------------------------------------------------------------
// Busy and file drop

function BusyOverlay(props: { busy: ReadonlyState<ShellBusy | undefined> }) {
  return (
    <Show when={derive(props.busy, (b) => b)} reveal="fade">
      {(b) => {
        const bar = (<span />) as HTMLElement;
        const progress = (
          <div
            class={classes.progress}
            role="progressbar"
            aria-label={b.label}
            data-indeterminate={b.progress === undefined ? '' : undefined}
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
          <div class={classes.busy} role="alertdialog" aria-modal="true" aria-label={b.label}>
            <div class={classes.busyCard}>
              <div class={classes.busyRow}>
                <span style="flex: 1">{derive(props.busy, (x) => x?.label ?? '')}</span>
                {b.cancel ? (
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
function FileDrop(props: { onDrop: (files: File[]) => void; label?: string }) {
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
    if (hasFiles(e)) e.preventDefault();
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
      <Show when={over}>{() => <div class={classes.drop}>{props.label ?? 'Drop files to open them'}</div>}</Show>
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
    [...dockedGroupsOf(l)].flatMap((g) => shownViews(shell, g)),
  );
  const shown = derive([shell.layout, shell.card, shell.absent, shell.views], (l, card) => {
    if (card) return { ref: card.ref, kind: 'card' as const };
    const top = stacked(l).at(-1);
    if (top) {
      const a = activeOf(shell, top.group);
      if (a) return { ref: a, kind: 'window' as const };
    }
    for (const g of dockedGroupsOf(l)) {
      if (g.collapsed) continue;
      const a = activeOf(shell, g);
      if (a) return { ref: a, kind: 'docked' as const };
    }
    return undefined;
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
                        ? () => shell.close(s.ref)
                        : () => {
                            const g = findView(shell.layout.value, s.ref)?.group;
                            if (g) shell.collapse(g.id);
                          }
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
                    aria-selected={derive(shown, (s) => String(s?.kind === 'docked' && s.ref === ref))}
                    onClick={() => shell.select(ref)}
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
                      return { id: ref, label: def.title, icon: def.icon, run: () => shell.select(ref) };
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

function dockedGroupsOf(l: ShellLayout): GroupNode[] {
  const out: GroupNode[] = [];
  const walk = (n: LayoutNode) => {
    if (n.kind === 'group') out.push(n);
    else if (n.kind === 'split') for (const c of n.children) walk(c);
  };
  walk(l.root);
  return out;
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
