/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The shell's layout (shell §2): a value, and pure functions on it. A tree whose leaves are groups of views and the one
// main area, its inner nodes splits; beside it the floating windows (each a group with a rectangle). The main area is
// itself a region (shell §2.1): without a tree it is the app's content alone; with one, a tree of the same groups and
// splits around the one `content` leaf (the app's own main view), so views dock beside the app's stage, a stage is a
// view, and the stage splits. Every operation returns a new layout that keeps the invariants (`normalize`), so the model
// only ever holds a valid one, and a stored one comes back through `validateLayout`.

/** A view in the layout: its id, or `id#key` for one instance of a view that has several (shell §2). */
export type ViewRef = string;

export const LayoutVersion = 1;

/** A group: views shown one at a time, a strip to switch them. */
export interface GroupNode {
  readonly kind: 'group';
  readonly id: string;
  readonly views: readonly ViewRef[];
  /** The view it shows; undefined only when it has none to show. */
  readonly active?: ViewRef;
  /** Only its strip shows (a rail's panel closed, a tab row's content folded). */
  readonly collapsed?: boolean;
  /** Its size along its parent split, px; undefined: it takes the rest. */
  readonly size?: number;
}

/**
 * The app's main area: exactly one in a layout, in the tree. Its region (shell §2.1): `root`, a tree of groups and
 * splits around the one `content` leaf (none: the content alone); `top`, the views overlaid on the content's top (the
 * app places them in its content, e.g. a stage's top slot, shell §2.1).
 */
export interface MainNode {
  readonly kind: 'main';
  readonly id: 'main';
  readonly root?: LayoutNode;
  readonly top?: GroupNode;
}

/** The app's own main view in its region: exactly one there when the region has a tree; it takes the rest. */
export interface ContentNode {
  readonly kind: 'content';
  readonly id: 'content';
}

export interface SplitNode {
  readonly kind: 'split';
  readonly id: string;
  /** `row`: children side by side; `column`: stacked. */
  readonly dir: 'row' | 'column';
  readonly children: readonly LayoutNode[];
  readonly size?: number;
  /**
   * A grid (shell §2.1): its children are splits of the other direction with as many children each, lined up: resizing
   * a cell resizes the cells of its row (or column) in every child, so a line across the grid moves as one. Kept only
   * while the shape holds (a pane split again, one closed: no grid).
   */
  readonly grid?: true;
}

export type LayoutNode = GroupNode | MainNode | SplitNode | ContentNode;

export interface Rect {
  readonly x: number;
  readonly y: number;
  readonly width: number;
  readonly height: number;
}

/** A floating window: a group with its rectangle (relative to the main area) and its place in the stack. */
export interface WindowNode {
  readonly group: GroupNode;
  readonly rect: Rect;
  /** Higher is on top. */
  readonly z: number;
}

export interface ShellLayout {
  readonly version: typeof LayoutVersion;
  readonly root: LayoutNode;
  readonly windows: readonly WindowNode[];
  /** Views out of view by the person's choice (`hide`); they keep their place. */
  readonly hidden: readonly ViewRef[];
  /** Whether the rail (and the docked column with it) shows. */
  readonly rail: boolean;
  /** The next id for a split or a group the layout makes. */
  readonly seq: number;
}

/** The right group's id: the default home of a view. */
export const RightGroup = 'right';
/** The id of the group overlaid on the top of the main area's content (shell §2.1). */
export const TopGroup = 'top';

/** Smallest size a docked group takes (px). */
export const MinGroupSize = 180;
/**
 * Smallest size a pane of the main area's region takes (px, shell §2.1): a strip docked above the content is a few rows,
 * not a panel.
 */
export const MinPaneSize = 64;
/** The least width the main area keeps when docked sizes grow (px). */
export const MinMainSize = 240;
/** The rail's width (px; `--msui-shell-rail`). */
export const RailWidth = 40;
/** A docked group's size when nothing gave it one (px). */
export const DefaultGroupSize = 280;
/** The right column's width when it is made (px; `ShellOptions.size`), and the card's without a column. */
export const DefaultColumnSize = 340;
/** A group split below another (Split below, px). */
export const DefaultSplitSize = 240;
/** A `top` view's height in a split above the content (px: its header and a few rows; shell §2.1). */
export const DefaultTopSize = 96;
/** A window's size when its view gives none (px). */
export const DefaultWindowSize = { width: 360, height: 320 } as const;
/** Smallest size a window takes (px). */
export const MinWindow = { width: 240, height: 160 } as const;

/** The default layout: the main area alone (the right group is made as the first view is placed in it). */
export function defaultLayout(): ShellLayout {
  return { version: LayoutVersion, root: { kind: 'main', id: 'main' }, windows: [], hidden: [], rail: true, seq: 1 };
}

// ---------------------------------------------------------------------------------------------------------------------
// Reading

/**
 * What the layout's readers ask of it, looked up once per layout value (`layoutIndex`, kept by the value: a layout is
 * immutable, so every reader of one value, the model's derived states and the view's trees, shares one index). Why: a
 * drag writes a layout a frame, and each tree, group, splitter and port asked again by walking it (O(n²) a frame).
 */
export interface LayoutIndex {
  /** The main area's node (the first in the tree; a bare one when it has none). */
  readonly main: MainNode;
  /** Every group in the tree, in document order (`dockedGroups`). */
  readonly docked: readonly GroupNode[];
  /** Every node of the tree (the main area's region and its top group included) by id; the first of an id wins. */
  readonly nodes: ReadonlyMap<string, LayoutNode>;
  /** Each node's parent split (or the main area, for its region's root and its top group). */
  readonly parents: ReadonlyMap<string, LayoutNode>;
  /** Where each view is (docked first, then the windows). */
  readonly views: ReadonlyMap<ViewRef, ViewPlace>;
  /** Each group by id (docked first, then the windows). */
  readonly groups: ReadonlyMap<string, ViewPlace>;
  readonly rail: GroupNode | undefined;
  readonly rightColumn: LayoutNode | undefined;
  /** The groups of the right column. */
  readonly rightGroups: ReadonlySet<string>;
  /** The groups of the main area's region and its top group. */
  readonly regionGroups: ReadonlySet<string>;
  /** Every node of the main area's region (panes and splits). */
  readonly regionNodes: ReadonlySet<string>;
  /** The column split each group is a direct child of. */
  readonly columns: ReadonlyMap<string, SplitNode>;
  /** The cells at the same place in the other lines of their grid (`gridPeers`), for the cells of a grid only. */
  readonly gridPeers: ReadonlyMap<string, readonly string[]>;
}

const indexes = new WeakMap<ShellLayout, LayoutIndex>();

/** The layout's index (`LayoutIndex`), built at the first ask of a layout value and kept with it. */
export function layoutIndex(layout: ShellLayout): LayoutIndex {
  let ix = indexes.get(layout);
  if (!ix) indexes.set(layout, (ix = buildIndex(layout)));
  return ix;
}

function buildIndex(layout: ShellLayout): LayoutIndex {
  const docked: GroupNode[] = [];
  const nodes = new Map<string, LayoutNode>();
  const parents = new Map<string, LayoutNode>();
  const columns = new Map<string, SplitNode>();
  const peers = new Map<string, string[]>();
  let main: MainNode | undefined;
  const walk = (n: LayoutNode, parent: LayoutNode | undefined, nested: boolean) => {
    if (!nodes.has(n.id)) {
      nodes.set(n.id, n);
      if (parent) parents.set(n.id, parent);
    }
    if (n.kind === 'group') docked.push(n);
    else if (n.kind === 'split') {
      if (n.grid) gridLines(n, peers);
      for (const c of n.children) {
        if (n.dir === 'column' && c.kind === 'group') columns.set(c.id, n);
        walk(c, n, nested);
      }
    } else if (n.kind === 'main') {
      if (!main && !nested) main = n;
      // The main area's region, then its top overlay (shell §2.1).
      if (n.root) walk(n.root, n, true);
      if (n.top) walk(n.top, n, true);
    }
  };
  walk(layout.root, undefined, false);
  const groups = new Map<string, ViewPlace>();
  const views = new Map<ViewRef, ViewPlace>();
  const add = (place: ViewPlace) => {
    if (!groups.has(place.group.id)) groups.set(place.group.id, place);
    for (const v of place.group.views) if (!views.has(v)) views.set(v, place);
  };
  for (const g of docked) add({ group: g });
  for (const w of layout.windows) add({ group: w.group, window: w });
  const m = main ?? { kind: 'main', id: 'main' };
  const regionNodes = new Set<string>();
  const regionGroups = new Set<string>();
  const inRegion = (n: LayoutNode) => {
    regionNodes.add(n.id);
    if (n.kind === 'group') regionGroups.add(n.id);
    else if (n.kind === 'split') for (const c of n.children) inRegion(c);
  };
  if (m.root) inRegion(m.root);
  if (m.top) regionGroups.add(m.top.id);
  const column = findRightColumn(layout);
  const rightGroups = new Set<string>();
  const inColumn = (n: LayoutNode) => {
    if (n.kind === 'group') rightGroups.add(n.id);
    else if (n.kind === 'split') for (const c of n.children) inColumn(c);
  };
  if (column) inColumn(column);
  return {
    main: m,
    docked,
    nodes,
    parents,
    views,
    groups,
    rail: findRailGroup(layout),
    rightColumn: column,
    rightGroups,
    regionGroups,
    regionNodes,
    columns,
    gridPeers: peers,
  };
}

/** Each cell of a grid's lines: the cells at its place in the other lines (`LayoutIndex.gridPeers`). */
function gridLines(grid: SplitNode, peers: Map<string, string[]>): void {
  for (const line of grid.children) {
    if (line.kind !== 'split') continue;
    for (let i = 0; i < line.children.length; i++) {
      const list: string[] = [];
      for (const other of grid.children) {
        if (other !== line && other.kind === 'split' && other.children[i]) list.push(other.children[i].id);
      }
      peers.set(line.children[i].id, list);
    }
  }
}

/** Where a view is: its group, and whether that group is a window's. */
export interface ViewPlace {
  readonly group: GroupNode;
  readonly window?: WindowNode;
}

/** Every group in the tree, in document order (left to right, top to bottom). */
export function dockedGroups(layout: ShellLayout): readonly GroupNode[] {
  return layoutIndex(layout).docked;
}

/** The main area's node. */
export function mainNode(layout: ShellLayout): MainNode {
  return layoutIndex(layout).main;
}

/** A node of the tree by id (the main area's region and its top group included). */
export function nodeOf(layout: ShellLayout, id: string): LayoutNode | undefined {
  return layoutIndex(layout).nodes.get(id);
}

/** Whether a node of any kind (a pane, a split) is in the main area's region. */
export function regionHas(layout: ShellLayout, id: string): boolean {
  return layoutIndex(layout).regionNodes.has(id);
}

/** Whether a docked group is in the main area's region (beside the app's content, shell §2.1) or its top overlay. */
export function inRegion(layout: ShellLayout, groupId: string): boolean {
  return layoutIndex(layout).regionGroups.has(groupId);
}

/** Every group, docked then windows (in stacking order). */
export function allGroups(layout: ShellLayout): GroupNode[] {
  return [...dockedGroups(layout), ...stacked(layout).map((w) => w.group)];
}

/** The windows bottom to top. */
export function stacked(layout: ShellLayout): WindowNode[] {
  return stacked0(layout.windows);
}

export function findView(layout: ShellLayout, ref: ViewRef): ViewPlace | undefined {
  return layoutIndex(layout).views.get(ref);
}

export function findGroup(layout: ShellLayout, id: string): ViewPlace | undefined {
  return layoutIndex(layout).groups.get(id);
}

/** Every view in the layout (docked and in windows). */
export function placedViews(layout: ShellLayout): ViewRef[] {
  return allGroups(layout).flatMap((g) => g.views);
}

/**
 * The group whose strip is the rail (shell §2): the top group of the docked column at the shell's right edge (the last
 * child of the root row, or the root's right column). Undefined when no group is docked at the right edge.
 */
export function railGroup(layout: ShellLayout): GroupNode | undefined {
  return layoutIndex(layout).rail;
}

function findRailGroup(layout: ShellLayout): GroupNode | undefined {
  let n: LayoutNode = layout.root;
  for (;;) {
    if (n.kind === 'group') return n;
    if (n.kind === 'main' || n.kind === 'content') return undefined;
    const last: LayoutNode = n.dir === 'row' ? n.children[n.children.length - 1] : n.children[0];
    if (!last || containsMain(last)) return undefined;
    n = last;
  }
}

/** The docked column at the shell's right edge (the rail's group, or the column split it is the top of). */
export function rightColumn(layout: ShellLayout): LayoutNode | undefined {
  return layoutIndex(layout).rightColumn;
}

function findRightColumn(layout: ShellLayout): LayoutNode | undefined {
  const root = layout.root;
  if (root.kind !== 'split' || root.dir !== 'row') return undefined;
  const last = root.children[root.children.length - 1];
  return last && !containsMain(last) ? last : undefined;
}

export function containsMain(n: LayoutNode): boolean {
  return containsLeaf(n, 'main');
}

// ---------------------------------------------------------------------------------------------------------------------
// Invariants

/**
 * The layout as it must be (every operation ends here): exactly one main area, in the tree; no empty group (a group
 * whose views are all gone leaves; one holding only hidden views stays, to come back to); no split with fewer than two
 * children (one child takes its place, with the split's size); no split directly in a split of its direction (its
 * children join the parent); a group's active view one of its views (else its first that is not hidden); each view in one
 * place; windows' rectangles at least the minimum size; ids unique. A node that needs no change is kept as it is (the
 * same object), and so is every node above it whose children are all kept: a write that moves one node (a window's
 * move, a drag's size) rebuilds only the path to it, and readers compare nodes by identity.
 */
export function normalize(layout: ShellLayout): ShellLayout {
  const seen = new Set<ViewRef>();
  const ids = new Set<string>();
  const hidden = layout.hidden.length ? new Set(layout.hidden) : undefined;
  const isHidden = (v: ViewRef) => !!hidden?.has(v);
  let seq = layout.seq;
  const freshId = (prefix: string) => {
    let id: string;
    do id = `${prefix}${seq++}`;
    while (ids.has(id));
    return id;
  };
  const group = (g: GroupNode): GroupNode | undefined => {
    // A copy of the views only when one is dropped (another place has it, or it is no ref).
    let copy: ViewRef[] | undefined;
    for (let i = 0; i < g.views.length; i++) {
      const v = g.views[i];
      if (typeof v === 'string' && !seen.has(v)) {
        seen.add(v);
        copy?.push(v);
      } else if (!copy) copy = g.views.slice(0, i);
    }
    const views = copy ?? g.views;
    if (views.length === 0) return undefined;
    let active = g.active;
    if (active === undefined || isHidden(active) || !views.includes(active)) active = views.find((v) => !isHidden(v));
    const id = ids.has(g.id) ? freshId('g') : g.id;
    ids.add(id);
    const collapsed = g.collapsed || undefined;
    if (views === g.views && active === g.active && id === g.id && collapsed === g.collapsed && !hasUndefined(g)) {
      return g;
    }
    return clean({ ...g, id, views, active, collapsed });
  };
  // The tree's one leaf that takes the rest: the main area at the top, the app's content in the main area's region
  // (shell §2.1). Each level keeps exactly one of its leaf and none of the other's.
  const counts = { main: 0, content: 0 };
  const node = (n: LayoutNode, leaf: 'main' | 'content'): LayoutNode | undefined => {
    if (n.kind === 'main' || n.kind === 'content') {
      if (n.kind !== leaf || counts[leaf]++ > 0) return undefined;
      ids.add(n.id);
      if (n.kind === 'main') return mainArea(n);
      return hasOnly(n, 2) ? n : { kind: 'content', id: 'content' };
    }
    if (n.kind === 'group') return group(n);
    let children: LayoutNode[] | undefined;
    for (let i = 0; i < n.children.length; i++) {
      const c = n.children[i];
      const k = node(c, leaf);
      // A split of the same direction joins this one.
      const joins = k?.kind === 'split' && k.dir === n.dir;
      if ((k !== c || joins) && !children) children = n.children.slice(0, i) as LayoutNode[];
      if (!children || !k) continue;
      if (joins) children.push(...(k as SplitNode).children);
      else children.push(k);
    }
    const list = children ?? n.children;
    if (list.length === 0) return undefined;
    if (list.length === 1) {
      const only = list[0];
      if (only.kind === 'main' || only.kind === 'content') return only;
      const size = n.size ?? only.size;
      return size === only.size ? only : clean({ ...only, size });
    }
    const id = ids.has(n.id) ? freshId('s') : n.id;
    ids.add(id);
    // The child holding the level's leaf (or, without one, the first) takes the rest: it has no size. Every other
    // child has a size, at least the level's least (a group made without one gets the default).
    const flex = list.findIndex((c) => containsLeaf(c, leaf));
    const at = flex >= 0 ? flex : 0;
    const min = leaf === 'content' ? MinPaneSize : MinGroupSize;
    let sized: LayoutNode[] | undefined;
    for (let i = 0; i < list.length; i++) {
      const c = list[i];
      let k: LayoutNode = c;
      if (i === at || c.kind === 'main' || c.kind === 'content') k = withoutSize(c);
      else {
        const size = Math.max(min, Math.round(c.size ?? DefaultGroupSize));
        if (size !== c.size) k = { ...c, size };
      }
      if (k !== c && !sized) sized = list.slice(0, i);
      sized?.push(k);
    }
    const out = sized ?? list;
    // A grid stays one while its children line up: splits of the other direction, as many children each.
    const grid = n.grid && isGridShape(n.dir, out) ? (true as const) : undefined;
    if (out === n.children && id === n.id && grid === n.grid && !hasUndefined(n)) return n;
    return clean({ ...n, id, children: out, grid });
  };
  /** The main area with its region normalized: the content kept in it (at the left when it was lost), or no tree. */
  const mainArea = (m: MainNode): MainNode => {
    const outer = counts.content;
    counts.content = 0;
    let root = m.root ? node(m.root, 'content') : undefined;
    if (root && counts.content === 0) root = withLeaf(root, { kind: 'content', id: 'content' }, ids, freshId);
    counts.content = outer;
    // The content alone is no tree.
    if (root?.kind === 'content') root = undefined;
    // The top group: its own id, no size.
    const t =
      m.top &&
      (m.top.id === TopGroup && m.top.size === undefined ? m.top : { ...m.top, id: TopGroup, size: undefined });
    const top = t && group(t);
    if (root === m.root && top === m.top && !hasUndefined(m)) return m;
    return clean({ kind: 'main' as const, id: 'main' as const, root, top });
  };
  let root = node(layout.root, 'main');
  if (counts.main === 0) root = withLeaf(root, mainArea({ kind: 'main', id: 'main' }), ids, freshId);
  let windows: WindowNode[] | undefined;
  const order = stackOrder(layout.windows);
  for (let i = 0; i < order.length; i++) {
    const w = order[i];
    let g = group(w.group);
    if (g && (g.size !== undefined || g.collapsed !== undefined))
      g = clean({ ...g, size: undefined, collapsed: undefined });
    const rect = g && minRect(w.rect);
    const z = (windows?.length ?? i) + 1;
    const same = g === w.group && rect === w.rect && z === w.z && order === layout.windows;
    if (!same && !windows) windows = order.slice(0, i) as WindowNode[];
    if (windows && g && rect) windows.push(same ? w : { group: g, rect, z });
  }
  let hiddenList: readonly ViewRef[] = layout.hidden;
  for (let i = 0; i < layout.hidden.length; i++) {
    const v = layout.hidden[i];
    if (!seen.has(v) || layout.hidden.indexOf(v) !== i) {
      hiddenList = layout.hidden.filter((x, j, a) => seen.has(x) && a.indexOf(x) === j);
      break;
    }
  }
  if (
    root === layout.root &&
    !windows &&
    hiddenList === layout.hidden &&
    seq === layout.seq &&
    layout.version === LayoutVersion &&
    hasOnly(layout, 6)
  ) {
    return layout;
  }
  return {
    version: LayoutVersion,
    root: root!,
    windows: windows ?? layout.windows,
    hidden: hiddenList,
    rail: layout.rail,
    seq,
  };
}

/** The windows bottom to top: the same array when it is in order already. */
function stackOrder(windows: readonly WindowNode[]): readonly WindowNode[] {
  for (let i = 1; i < windows.length; i++) if (windows[i - 1].z > windows[i].z) return stacked0(windows);
  return windows;
}

function stacked0(windows: readonly WindowNode[]): WindowNode[] {
  return [...windows].sort((a, b) => a.z - b.z);
}

/** Whether an object has a field set to undefined (a layout is compared and serialized as plain data: none). */
function hasUndefined(o: object): boolean {
  for (const k in o) if ((o as Record<string, unknown>)[k] === undefined) return true;
  return false;
}

/** Whether an object has at most `n` fields. */
function hasOnly(o: object, n: number): boolean {
  let i = 0;
  for (const _ in o) if (++i > n) return false;
  return true;
}

/** A level's leaf put back at the left of what is left of its tree (the leaf alone when nothing is). */
function withLeaf(
  root: LayoutNode | undefined,
  leaf: MainNode | ContentNode,
  ids: Set<string>,
  freshId: (prefix: string) => string,
): LayoutNode {
  if (!root) return leaf;
  const id = leaf.kind === 'main' && !ids.has('root') ? 'root' : freshId('s');
  ids.add(id);
  return normalizeShallow({ kind: 'split', id, dir: 'row', children: [leaf, root] }, leaf.kind);
}

function normalizeShallow(split: SplitNode, leaf: 'main' | 'content'): SplitNode {
  const flex = split.children.findIndex((c) => containsLeaf(c, leaf));
  const min = leaf === 'content' ? MinPaneSize : MinGroupSize;
  return {
    ...split,
    children: split.children.map((c, i) => (i === flex ? withoutSize(c) : withLeastSize(c, min))),
  };
}

/** A sized node at least `min` (the level's least size; a stored size is clamped here, not when read). */
function withLeastSize(n: LayoutNode, min: number): LayoutNode {
  return (n.kind === 'group' || n.kind === 'split') && n.size !== undefined && n.size < min ? { ...n, size: min } : n;
}

/** Whether a node holds a level's leaf (the main area, or the app's content in the main area's region). */
export function containsLeaf(n: LayoutNode, leaf: 'main' | 'content'): boolean {
  if (n.kind === leaf) return true;
  return n.kind === 'split' && n.children.some((c) => containsLeaf(c, leaf));
}

function withoutSize<N extends LayoutNode>(n: N): N {
  if (n.kind === 'main' || n.kind === 'content' || (n as GroupNode | SplitNode).size === undefined) return n;
  const { size: _, ...rest } = n as GroupNode | SplitNode;
  return rest as N;
}

/** Drops undefined fields (a layout is compared and serialized as plain data). */
function clean<T extends object>(o: T): T {
  const out: any = {};
  for (const [k, v] of Object.entries(o)) if (v !== undefined) out[k] = v;
  return out;
}

/** A window's rectangle in whole px, at least `MinWindow` (the same object when it is already). */
function minRect(r: Rect): Rect {
  const x = Math.round(r.x);
  const y = Math.round(r.y);
  const width = Math.max(MinWindow.width, Math.round(r.width));
  const height = Math.max(MinWindow.height, Math.round(r.height));
  return x === r.x && y === r.y && width === r.width && height === r.height ? r : { x, y, width, height };
}

// ---------------------------------------------------------------------------------------------------------------------
// Operations

/** Where a view goes. */
export type LayoutTarget =
  /** Into a group (at `index` of its views; default the end). */
  | { readonly kind: 'group'; readonly id: string; readonly index?: number }
  /**
   * Beside a node (a group, the main area, the app's content, a split): a new group on that side, `size` px (default
   * `DefaultGroupSize`; at least the level's least, `normalize`).
   */
  | {
      readonly kind: 'edge';
      readonly of: string;
      readonly side: 'left' | 'right' | 'top' | 'bottom';
      readonly size?: number;
      readonly id?: string;
    }
  /** A new window. */
  | { readonly kind: 'float'; readonly rect: Rect; readonly id?: string }
  /** Overlaid on the top of the main area's content (shell §2.1). */
  | { readonly kind: 'top' };

/**
 * Maps every node of the tree (bottom up); `fn` returns the node or undefined to remove it. A split whose children all
 * come back as they were is handed to `fn` as it is (the same object): an operation that changes one node rebuilds only
 * the path to it.
 */
function mapTree(n: LayoutNode, fn: (n: LayoutNode) => LayoutNode | undefined): LayoutNode | undefined {
  if (n.kind === 'split') {
    let children: LayoutNode[] | undefined;
    for (let i = 0; i < n.children.length; i++) {
      const c = n.children[i];
      const k = mapTree(c, fn);
      if (k !== c && !children) children = n.children.slice(0, i) as LayoutNode[];
      if (children && k) children.push(k);
    }
    return fn(children ? { ...n, children } : n);
  }
  if (n.kind === 'main' && (n.root || n.top)) {
    // Into the main area's region and top overlay (shell §2.1); its group removed, no overlay.
    const root = n.root && mapTree(n.root, fn);
    const t = n.top && fn(n.top);
    const top = t?.kind === 'group' ? t : undefined;
    return fn(root === n.root && top === n.top ? n : clean({ ...n, root, top }));
  }
  return fn(n);
}

/** Maps every group, docked and in windows; the layout as it is when `fn` returns every group as it was. */
function mapGroups(layout: ShellLayout, fn: (g: GroupNode) => GroupNode | undefined): ShellLayout {
  const root = mapTree(layout.root, (n) => (n.kind === 'group' ? fn(n) : n)) ?? { kind: 'main', id: 'main' };
  let windows: WindowNode[] | undefined;
  for (let i = 0; i < layout.windows.length; i++) {
    const w = layout.windows[i];
    const g = fn(w.group);
    if (g !== w.group && !windows) windows = layout.windows.slice(0, i) as WindowNode[];
    if (windows && g) windows.push(g === w.group ? w : { ...w, group: g });
  }
  if (root === layout.root && !windows) return layout;
  return { ...layout, root, windows: windows ?? layout.windows };
}

/** The layout without `ref` in any group (a group left empty goes, a split left with one child collapses). */
export function removeView(layout: ShellLayout, ref: ViewRef): ShellLayout {
  return normalize(
    mapGroups(layout, (g) =>
      g.views.includes(ref)
        ? { ...g, views: g.views.filter((v) => v !== ref), active: g.active === ref ? undefined : g.active }
        : g,
    ),
  );
}

/** `ref` placed at `target` (moved, if it is somewhere), and made its group's active view. */
export function moveView(layout: ShellLayout, ref: ViewRef, target: LayoutTarget): ShellLayout {
  // Into its own group, at an index: a reorder.
  const from = findView(layout, ref);
  let next = from ? withoutViewKeepingGroups(layout, ref) : layout;
  let seq = next.seq;
  if (target.kind === 'group') {
    const place = findGroup(next, target.id);
    if (!place) return layout;
    next = mapGroups(next, (g) => {
      if (g.id !== target.id) return g;
      const views = [...g.views];
      const at = target.index === undefined ? views.length : Math.max(0, Math.min(views.length, target.index));
      views.splice(at, 0, ref);
      return { ...g, views, active: ref, collapsed: false };
    });
  } else if (target.kind === 'top') {
    const main = mainNode(next);
    const top: GroupNode = main.top
      ? { ...main.top, views: [...main.top.views, ref], active: ref, collapsed: false }
      : { kind: 'group', id: TopGroup, views: [ref], active: ref };
    next = { ...next, root: mapTree(next.root, (n) => (n.kind === 'main' ? { ...n, top } : n))! };
  } else if (target.kind === 'float') {
    const id = target.id ?? `w${seq++}`;
    const top = layout.windows.reduce((m, w) => Math.max(m, w.z), 0);
    next = {
      ...next,
      seq,
      windows: [
        ...next.windows,
        { group: { kind: 'group', id, views: [ref], active: ref }, rect: target.rect, z: top + 1 },
      ],
    };
  } else {
    const id = target.id ?? `g${seq++}`;
    const group: GroupNode = { kind: 'group', id, views: [ref], active: ref, size: target.size };
    const dir = target.side === 'left' || target.side === 'right' ? 'row' : 'column';
    const before = target.side === 'left' || target.side === 'top';
    // Beside the app's content in a region with no tree yet: the tree starts as the content alone.
    if (target.of === 'content' && !mainNode(next).root) {
      next = {
        ...next,
        root: mapTree(next.root, (n) => (n.kind === 'main' ? { ...n, root: { kind: 'content', id: 'content' } } : n))!,
      };
    }
    let placed = false;
    const root = mapTree(next.root, (n) => {
      if (n.id !== target.of || placed) return n;
      placed = true;
      const children = before ? [group, withoutSize(n)] : [withoutSize(n), group];
      const size = n.kind === 'group' || n.kind === 'split' ? n.size : undefined;
      // The main area's node holds its region: beside it is beside the whole main area (its region goes with it).
      return { kind: 'split', id: `s${seq++}`, dir, children, size };
    });
    if (!placed) return layout;
    next = { ...next, root: root!, seq };
  }
  return normalize({ ...next, seq: Math.max(seq, next.seq), hidden: next.hidden.filter((v) => v !== ref) });
}

/** `ref` out of its group, the group kept even when empty (a move into it, or next to it, still finds it). */
function withoutViewKeepingGroups(layout: ShellLayout, ref: ViewRef): ShellLayout {
  return mapGroups(layout, (g) =>
    g.views.includes(ref)
      ? { ...g, views: g.views.filter((v) => v !== ref), active: g.active === ref ? undefined : g.active }
      : g,
  );
}

/** `ref` its group's active view (its group opened, the view shown again if it was hidden). */
export function selectView(layout: ShellLayout, ref: ViewRef): ShellLayout {
  if (!findView(layout, ref)) return layout;
  const hidden = layout.hidden.filter((v) => v !== ref);
  return normalize({
    ...mapGroups({ ...layout, hidden }, (g) => (g.views.includes(ref) ? { ...g, active: ref, collapsed: false } : g)),
    hidden,
  });
}

/** A group's active view and fold set as they are given (placing a view is not opening it). */
export function setGroupState(
  layout: ShellLayout,
  groupId: string,
  active: ViewRef | undefined,
  collapsed: boolean,
): ShellLayout {
  return normalize(mapGroups(layout, (g) => (g.id === groupId ? { ...g, active, collapsed } : g)));
}

export function setCollapsed(layout: ShellLayout, groupId: string, collapsed: boolean): ShellLayout {
  return normalize(mapGroups(layout, (g) => (g.id === groupId ? { ...g, collapsed } : g)));
}

/** A docked node's size (px, at least `MinGroupSize`; the node that takes the rest has none). */
export function resizeNode(layout: ShellLayout, id: string, size: number): ShellLayout {
  const s = Math.max(regionHas(layout, id) ? MinPaneSize : MinGroupSize, Math.round(size));
  // In a grid's line: the cell at the same place in every line resizes with it (the line across moves as one).
  const peers = gridPeers(layout, id);
  const root = mapTree(layout.root, (n) =>
    (n.id === id || peers.includes(n.id)) &&
    (n.kind === 'group' || n.kind === 'split') &&
    n.size !== undefined &&
    n.size !== s
      ? { ...n, size: s }
      : n,
  );
  return root === layout.root ? layout : normalize({ ...layout, root: root! });
}

/** Several docked nodes' sizes at once (px; a node that takes the rest keeps none; the level's least size kept). */
export function setSizes(layout: ShellLayout, sizes: ReadonlyMap<string, number>): ShellLayout {
  const root = mapTree(layout.root, (n) => {
    const s = sizes.get(n.id);
    return s !== undefined && (n.kind === 'group' || n.kind === 'split') && n.size !== undefined && n.size !== s
      ? { ...n, size: s }
      : n;
  });
  return root === layout.root ? layout : normalize({ ...layout, root: root! });
}

function isGridShape(dir: 'row' | 'column', children: readonly LayoutNode[]): boolean {
  const first = children[0];
  if (children.length < 2 || first?.kind !== 'split' || first.dir === dir) return false;
  const n = first.children.length;
  return children.every((c) => c.kind === 'split' && c.dir !== dir && c.children.length === n);
}

const NoPeers: readonly string[] = [];

/** The ids of the cells at the same place as `id` in the other lines of its grid (none outside a grid; indexed). */
export function gridPeers(layout: ShellLayout, id: string): readonly string[] {
  return layoutIndex(layout).gridPeers.get(id) ?? NoPeers;
}

/**
 * The main area's region as a grid (shell §2.1): its panes that show (the app's content first, then the docked groups
 * in order; `shows`, default every group) in `columns` columns (default the square root, rounded up), filled top to
 * bottom; equal sizes over `area` (px). A full grid is lined up (`grid`); one with a short last column is not (its lines
 * would not match). A group that does not show (its views flagged off) keeps a place after the grid, taking no room.
 */
export function arrangeGrid(
  layout: ShellLayout,
  area: { width: number; height: number },
  columns?: number,
  shows: (g: GroupNode) => boolean = () => true,
): ShellLayout {
  const main = mainNode(layout);
  if (!main.root) return layout;
  const panes: LayoutNode[] = [];
  const rest: LayoutNode[] = [];
  const walk = (n: LayoutNode) => {
    if (n.kind === 'split') for (const c of n.children) walk(c);
    else if (n.kind === 'content') panes.unshift(withoutSize(n));
    else if (n.kind === 'group') (shows(n) ? panes : rest).push(withoutSize(n));
  };
  walk(main.root);
  const cols = Math.max(1, Math.min(panes.length, Math.round(columns ?? Math.ceil(Math.sqrt(panes.length)))));
  const rows = Math.ceil(panes.length / cols);
  const width = Math.round(area.width / cols) || DefaultGroupSize;
  const height = Math.round(area.height / rows) || DefaultGroupSize;
  let seq = layout.seq;
  const lines: LayoutNode[] = [];
  for (let c = 0; c < cols; c++) {
    const cells = panes.slice(c * rows, (c + 1) * rows).map((p) => ({ ...p, size: height }) as LayoutNode);
    if (!cells.length) continue;
    lines.push(
      cells.length === 1
        ? ({ ...cells[0], size: width } as LayoutNode)
        : { kind: 'split', id: `s${seq++}`, dir: 'column', children: cells, size: width },
    );
  }
  lines.push(...rest);
  const root: LayoutNode =
    lines.length === 1 ? lines[0] : { kind: 'split', id: `s${seq++}`, dir: 'row', children: lines, grid: true };
  const next = mapTree(layout.root, (n) => (n.kind === 'main' ? { ...n, root } : n))!;
  return normalize({ ...layout, root: next, seq });
}

export function setWindowRect(layout: ShellLayout, groupId: string, rect: Rect): ShellLayout {
  return normalize({ ...layout, windows: layout.windows.map((w) => (w.group.id === groupId ? { ...w, rect } : w)) });
}

/** The window on top of the others. */
export function raiseWindow(layout: ShellLayout, groupId: string): ShellLayout {
  const top = layout.windows.reduce((m, w) => Math.max(m, w.z), 0);
  const w = layout.windows.find((x) => x.group.id === groupId);
  if (!w || w.z === top) return layout;
  return normalize({ ...layout, windows: layout.windows.map((x) => (x === w ? { ...x, z: top + 1 } : x)) });
}

/** A window's views merged into another group (a window docked into the right group, into another window). */
export function mergeGroup(layout: ShellLayout, fromId: string, intoId: string): ShellLayout {
  const from = findGroup(layout, fromId);
  if (!from || fromId === intoId || !findGroup(layout, intoId)) return layout;
  let next = layout;
  for (const v of from.group.views) next = moveView(next, v, { kind: 'group', id: intoId });
  // Views the person hid stay hidden (a move shows a view: not these).
  next = keepHidden(next, layout.hidden);
  return from.group.active && !layout.hidden.includes(from.group.active) ? selectView(next, from.group.active) : next;
}

/** `layout` with the views of `hidden` that it still places hidden again (after moves, which show what they move). */
export function keepHidden(layout: ShellLayout, hidden: readonly ViewRef[]): ShellLayout {
  const placed = new Set(placedViews(layout));
  const again = hidden.filter((v) => placed.has(v) && !layout.hidden.includes(v));
  return again.length ? normalize({ ...layout, hidden: [...layout.hidden, ...again] }) : layout;
}

/** Whether a docked group is in the right column (the rail's group and those under it). */
export function inRightColumn(layout: ShellLayout, groupId: string): boolean {
  return layoutIndex(layout).rightGroups.has(groupId);
}

/** Moves a view one place up or down among the docked groups of its column (Move up, Move down). */
export function moveInColumn(layout: ShellLayout, ref: ViewRef, dir: -1 | 1): ShellLayout {
  const place = findView(layout, ref);
  if (!place || place.window) return layout;
  const column = columnOf(layout, place.group.id);
  if (!column) return layout;
  const groups = column.children.filter((c): c is GroupNode => c.kind === 'group');
  const i = groups.findIndex((g) => g.id === place.group.id);
  const to = groups[i + dir];
  if (!to) return layout;
  return moveView(layout, ref, { kind: 'group', id: to.id });
}

/** The column split a docked group is in, if it is in one. */
export function columnOf(layout: ShellLayout, groupId: string): SplitNode | undefined {
  return layoutIndex(layout).columns.get(groupId);
}

export function setHidden(layout: ShellLayout, ref: ViewRef, hidden: boolean): ShellLayout {
  const list = layout.hidden.filter((v) => v !== ref);
  return normalize({ ...layout, hidden: hidden ? [...list, ref] : list });
}

export function setRail(layout: ShellLayout, shown: boolean): ShellLayout {
  return layout.rail === shown ? layout : { ...layout, rail: shown };
}

// ---------------------------------------------------------------------------------------------------------------------
// Stored layouts

/**
 * A stored layout back as a layout, or undefined when it is not one (another version, or any field of a wrong shape or
 * missing: no field is filled in, a stored layout is the value `save` gave or nothing), then normalized (shell §2,
 * persistence). Unknown view refs are kept: a view that registers later lands where it was left.
 */
export function validateLayout(value: unknown): ShellLayout | undefined {
  if (!isObject(value) || value.version !== LayoutVersion) return undefined;
  const { windows, hidden, rail, seq } = value;
  if (!Array.isArray(windows) || !isStrings(hidden) || typeof rail !== 'boolean' || !isNumber(seq)) return undefined;
  const root = validNode(value.root, 0);
  if (!root) return undefined;
  const wins: WindowNode[] = [];
  for (const w of windows) {
    if (!isObject(w) || !isNumber(w.z) || !isObject(w.rect)) return undefined;
    const g = validNode(w.group, 0);
    const { x, y, width, height } = w.rect;
    if (g?.kind !== 'group' || !isNumber(x) || !isNumber(y) || !isNumber(width) || !isNumber(height)) return undefined;
    wins.push({ group: g, rect: { x, y, width, height }, z: w.z });
  }
  return normalize({ version: LayoutVersion, root, windows: wins, hidden, rail, seq: Math.max(1, Math.floor(seq)) });
}

/** A stored node, or undefined when it or anything in it is not of its shape. */
function validNode(v: unknown, depth: number): LayoutNode | undefined {
  if (!isObject(v) || depth > 16) return undefined;
  // The least size is the level's (`normalize`: a pane's is `MinPaneSize`, a docked group's `MinGroupSize`).
  const size = v.size;
  if (size !== undefined && !isNumber(size)) return undefined;
  if (v.kind === 'main') {
    const root = v.root === undefined ? undefined : validNode(v.root, depth + 1);
    const top = v.top === undefined ? undefined : validNode(v.top, depth + 1);
    if ((v.root !== undefined && !root) || (v.top !== undefined && top?.kind !== 'group')) return undefined;
    return clean({ kind: 'main' as const, id: 'main' as const, root, top: top as GroupNode | undefined });
  }
  if (v.kind === 'content') return { kind: 'content', id: 'content' };
  if (typeof v.id !== 'string') return undefined;
  if (v.kind === 'group') {
    const { views, active, collapsed } = v;
    if (!isStrings(views) || (active !== undefined && typeof active !== 'string')) return undefined;
    if (collapsed !== undefined && typeof collapsed !== 'boolean') return undefined;
    return clean({ kind: 'group' as const, id: v.id, views, active, collapsed: collapsed || undefined, size });
  }
  if (v.kind === 'split' && (v.dir === 'row' || v.dir === 'column') && Array.isArray(v.children)) {
    if (v.grid !== undefined && typeof v.grid !== 'boolean') return undefined;
    const children: LayoutNode[] = [];
    for (const c of v.children) {
      const n = validNode(c, depth + 1);
      if (!n) return undefined;
      children.push(n);
    }
    const grid = v.grid === true ? (true as const) : undefined;
    return clean({ kind: 'split' as const, id: v.id, dir: v.dir as 'row' | 'column', children, size, grid });
  }
  return undefined;
}

function isNumber(v: unknown): v is number {
  return typeof v === 'number' && Number.isFinite(v);
}

function isStrings(v: unknown): v is string[] {
  return Array.isArray(v) && v.every((x) => typeof x === 'string');
}

function isObject(v: unknown): v is Record<string, unknown> {
  return typeof v === 'object' && v !== null && !Array.isArray(v);
}

// ---------------------------------------------------------------------------------------------------------------------
// Window placement

/**
 * A window's rectangle kept inside the area it floats over (shell §2): at least `MinWindow`, no larger than the area
 * (unless the area is smaller than `MinWindow`), and placed so that all of it is in (its top-left at least, when the
 * area is smaller): its header stays reachable.
 */
export function clampRect(rect: Rect, area: { width: number; height: number }): Rect {
  const width = Math.min(Math.max(MinWindow.width, rect.width), Math.max(MinWindow.width, area.width));
  const height = Math.min(Math.max(MinWindow.height, rect.height), Math.max(MinWindow.height, area.height));
  const x = Math.min(Math.max(0, rect.x), Math.max(0, area.width - width));
  const y = Math.min(Math.max(0, rect.y), Math.max(0, area.height - height));
  return { x: Math.round(x), y: Math.round(y), width: Math.round(width), height: Math.round(height) };
}

/** Px a new window is offset from the one before it (the cascade). */
export const CascadeStep = 24;
/** Px a window keeps from the area's edge when it is placed. */
export const WindowMargin = 16;

/**
 * Where a new window opens (shell §2): beside `anchor` (a button's rectangle, relative to the area: to its left, its top
 * level with the anchor's, else below it), else cascaded from the area's top-left past the windows already there;
 * clamped into the area.
 */
export function placeWindow(
  size: { width: number; height: number },
  area: { width: number; height: number },
  existing: readonly Rect[],
  anchor?: Rect,
): Rect {
  if (anchor) {
    const left = anchor.x - size.width - WindowMargin / 2;
    const x = left >= 0 ? left : anchor.x;
    const y = left >= 0 ? anchor.y : anchor.y + anchor.height + WindowMargin / 2;
    return clampRect({ x, y, ...size }, area);
  }
  let x = WindowMargin;
  let y = WindowMargin;
  for (let i = 0; i < 64 && existing.some((r) => Math.abs(r.x - x) < 4 && Math.abs(r.y - y) < 4); i++) {
    x += CascadeStep;
    y += CascadeStep;
    if (x + size.width > area.width || y + size.height > area.height) {
      x = WindowMargin + ((i + 1) % 4) * 8;
      y = WindowMargin;
    }
  }
  return clampRect({ x, y, ...size }, area);
}
