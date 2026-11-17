/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The shell's layout (shell §2): a value, and pure functions on it. A tree whose leaves are groups of views and the one
// main area, its inner nodes splits; beside it the floating windows (each a group with a rectangle). Every operation
// returns a new layout that keeps the invariants (`normalize`), so the model only ever holds a valid one, and a stored
// one comes back through `validateLayout`.

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

/** The app's main area: exactly one in a layout, in the tree. */
export interface MainNode {
  readonly kind: 'main';
  readonly id: 'main';
}

export interface SplitNode {
  readonly kind: 'split';
  readonly id: string;
  /** `row`: children side by side; `column`: stacked. */
  readonly dir: 'row' | 'column';
  readonly children: readonly LayoutNode[];
  readonly size?: number;
}

export type LayoutNode = GroupNode | MainNode | SplitNode;

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

/** Smallest sizes a docked group or a window takes (px). */
export const MinGroupSize = 180;
/** The least width the main area keeps when docked sizes grow (px). */
export const MinMainSize = 240;
/** The rail's width (px; `--msui-shell-rail`). */
export const RailWidth = 40;
/** A docked group's size when nothing gave it one (px). */
export const DefaultGroupSize = 280;
export const MinWindow = { width: 240, height: 160 } as const;

/** The default layout: the main area, the right group beside it (340 wide). */
export function defaultLayout(
  rightViews: readonly ViewRef[] = [],
  options?: { size?: number; active?: ViewRef },
): ShellLayout {
  return normalize({
    version: LayoutVersion,
    root: {
      kind: 'split',
      id: 'root',
      dir: 'row',
      children: [
        { kind: 'main', id: 'main' },
        { kind: 'group', id: RightGroup, views: rightViews, active: options?.active, size: options?.size ?? 340 },
      ],
    },
    windows: [],
    hidden: [],
    rail: true,
    seq: 1,
  });
}

// ---------------------------------------------------------------------------------------------------------------------
// Reading

/** Where a view is: its group, and whether that group is a window's. */
export interface ViewPlace {
  readonly group: GroupNode;
  readonly window?: WindowNode;
}

/** Every group in the tree, in document order (left to right, top to bottom). */
export function dockedGroups(layout: ShellLayout): GroupNode[] {
  const out: GroupNode[] = [];
  const walk = (n: LayoutNode) => {
    if (n.kind === 'group') out.push(n);
    else if (n.kind === 'split') for (const c of n.children) walk(c);
  };
  walk(layout.root);
  return out;
}

/** Every group, docked then windows (in stacking order). */
export function allGroups(layout: ShellLayout): GroupNode[] {
  return [...dockedGroups(layout), ...stacked(layout).map((w) => w.group)];
}

/** The windows bottom to top. */
export function stacked(layout: ShellLayout): WindowNode[] {
  return [...layout.windows].sort((a, b) => a.z - b.z);
}

export function findView(layout: ShellLayout, ref: ViewRef): ViewPlace | undefined {
  for (const g of dockedGroups(layout)) if (g.views.includes(ref)) return { group: g };
  for (const w of layout.windows) if (w.group.views.includes(ref)) return { group: w.group, window: w };
  return undefined;
}

export function findGroup(layout: ShellLayout, id: string): ViewPlace | undefined {
  for (const g of dockedGroups(layout)) if (g.id === id) return { group: g };
  for (const w of layout.windows) if (w.group.id === id) return { group: w.group, window: w };
  return undefined;
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
  let n: LayoutNode = layout.root;
  for (;;) {
    if (n.kind === 'group') return n;
    if (n.kind === 'main') return undefined;
    const last: LayoutNode = n.dir === 'row' ? n.children[n.children.length - 1] : n.children[0];
    if (!last || containsMain(last)) return undefined;
    n = last;
  }
}

/** The docked column at the shell's right edge (the rail's group, or the column split it is the top of). */
export function rightColumn(layout: ShellLayout): LayoutNode | undefined {
  const root = layout.root;
  if (root.kind !== 'split' || root.dir !== 'row') return undefined;
  const last = root.children[root.children.length - 1];
  return last && !containsMain(last) ? last : undefined;
}

export function containsMain(n: LayoutNode): boolean {
  if (n.kind === 'main') return true;
  if (n.kind === 'group') return false;
  return n.children.some(containsMain);
}

// ---------------------------------------------------------------------------------------------------------------------
// Invariants

/**
 * The layout as it must be (every operation ends here): exactly one main area, in the tree; no empty group (a group
 * whose views are all gone leaves; one holding only hidden views stays, to come back to); no split with fewer than two
 * children (one child takes its place, with the split's size); no split directly in a split of its direction (its
 * children join the parent); a group's active view one of its views (else its first that is not hidden); each view in one
 * place; windows' rectangles at least the minimum size; ids unique.
 */
export function normalize(layout: ShellLayout): ShellLayout {
  const seen = new Set<ViewRef>();
  const ids = new Set<string>();
  let seq = layout.seq;
  const freshId = (prefix: string) => {
    let id: string;
    do id = `${prefix}${seq++}`;
    while (ids.has(id));
    return id;
  };
  const group = (g: GroupNode): GroupNode | undefined => {
    const views = g.views.filter((v) => typeof v === 'string' && !seen.has(v));
    for (const v of views) seen.add(v);
    if (views.length === 0) return undefined;
    const shown = views.filter((v) => !layout.hidden.includes(v));
    const active = g.active && shown.includes(g.active) ? g.active : shown[0];
    const id = ids.has(g.id) ? freshId('g') : g.id;
    ids.add(id);
    return clean({ ...g, id, views, active, collapsed: g.collapsed || undefined });
  };
  let mains = 0;
  const node = (n: LayoutNode): LayoutNode | undefined => {
    if (n.kind === 'main') {
      if (mains++ > 0) return undefined;
      ids.add('main');
      return { kind: 'main', id: 'main' };
    }
    if (n.kind === 'group') return group(n);
    const children: LayoutNode[] = [];
    for (const c of n.children) {
      const k = node(c);
      if (!k) continue;
      // A split of the same direction joins this one.
      if (k.kind === 'split' && k.dir === n.dir) children.push(...k.children);
      else children.push(k);
    }
    if (children.length === 0) return undefined;
    if (children.length === 1) {
      const only = children[0];
      return only.kind === 'main' ? only : clean({ ...only, size: n.size ?? only.size });
    }
    const id = ids.has(n.id) ? freshId('s') : n.id;
    ids.add(id);
    // The child holding the main area (or, without one, the first) takes the rest: it has no size.
    const flex = children.findIndex(containsMain);
    const at = flex >= 0 ? flex : 0;
    // Every other child has a size, at least the minimum (a group made without one gets the default).
    const sized = children.map((c, i) =>
      i === at || c.kind === 'main'
        ? withoutSize(c)
        : { ...c, size: Math.max(MinGroupSize, Math.round(c.size ?? DefaultGroupSize)) },
    );
    return clean({ ...n, id, children: sized });
  };
  let root = node(layout.root);
  if (mains === 0) {
    // The main area is always there: put it back at the left.
    root = root
      ? {
          kind: 'split',
          id: ids.has('root') ? freshId('s') : 'root',
          dir: 'row',
          children: [{ kind: 'main', id: 'main' }, root],
        }
      : { kind: 'main', id: 'main' };
    if (root.kind === 'split') root = normalizeShallow(root);
  }
  const windows: WindowNode[] = [];
  for (const w of stacked(layout)) {
    const g = group(w.group);
    if (!g) continue;
    windows.push({
      group: { ...g, size: undefined, collapsed: undefined },
      rect: minRect(w.rect),
      z: windows.length + 1,
    });
  }
  return {
    version: LayoutVersion,
    root: root!,
    windows,
    hidden: layout.hidden.filter((v, i, a) => seen.has(v) && a.indexOf(v) === i),
    rail: layout.rail,
    seq,
  };
}

function normalizeShallow(split: SplitNode): SplitNode {
  const flex = split.children.findIndex(containsMain);
  return { ...split, children: split.children.map((c, i) => (i === flex ? withoutSize(c) : c)) };
}

function withoutSize<N extends LayoutNode>(n: N): N {
  if (n.kind === 'main' || (n as GroupNode | SplitNode).size === undefined) return n;
  const { size: _, ...rest } = n as GroupNode | SplitNode;
  return rest as N;
}

/** Drops undefined fields (a layout is compared and serialized as plain data). */
function clean<T extends object>(o: T): T {
  const out: any = {};
  for (const [k, v] of Object.entries(o)) if (v !== undefined) out[k] = v;
  return out;
}

function minRect(r: Rect): Rect {
  return {
    x: Math.round(r.x),
    y: Math.round(r.y),
    width: Math.max(MinWindow.width, Math.round(r.width)),
    height: Math.max(MinWindow.height, Math.round(r.height)),
  };
}

// ---------------------------------------------------------------------------------------------------------------------
// Operations

/** Where a view goes. */
export type LayoutTarget =
  /** Into a group (at `index` of its views; default the end). */
  | { readonly kind: 'group'; readonly id: string; readonly index?: number }
  /**
   * Beside a group or the main area: a new group on that side, `size` px (default half the room, at least
   * `MinGroupSize`; the caller passes the room it knows).
   */
  | {
      readonly kind: 'edge';
      readonly of: string;
      readonly side: 'left' | 'right' | 'top' | 'bottom';
      readonly size?: number;
      readonly id?: string;
    }
  /** A new window. */
  | { readonly kind: 'float'; readonly rect: Rect; readonly id?: string };

/** Maps every node of the tree (bottom up); `fn` returns the node or undefined to remove it. */
function mapTree(n: LayoutNode, fn: (n: LayoutNode) => LayoutNode | undefined): LayoutNode | undefined {
  if (n.kind === 'split') {
    const children = n.children.map((c) => mapTree(c, fn)).filter((c): c is LayoutNode => !!c);
    return fn({ ...n, children });
  }
  return fn(n);
}

function mapGroups(layout: ShellLayout, fn: (g: GroupNode) => GroupNode | undefined): ShellLayout {
  const root = mapTree(layout.root, (n) => (n.kind === 'group' ? fn(n) : n)) ?? { kind: 'main', id: 'main' };
  const windows = layout.windows
    .map((w) => {
      const g = fn(w.group);
      return g ? { ...w, group: g } : undefined;
    })
    .filter((w): w is WindowNode => !!w);
  return { ...layout, root, windows };
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
    let placed = false;
    const root = mapTree(next.root, (n) => {
      if (n.id !== target.of || placed) return n;
      placed = true;
      const children = before ? [group, withoutSize(n)] : [withoutSize(n), group];
      return { kind: 'split', id: `s${seq++}`, dir, children, size: (n as GroupNode | SplitNode).size };
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
  const s = Math.max(MinGroupSize, Math.round(size));
  const root = mapTree(layout.root, (n) =>
    n.id === id && n.kind !== 'main' && n.size !== undefined ? { ...n, size: s } : n,
  );
  return normalize({ ...layout, root: root! });
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
  const col = rightColumn(layout);
  if (!col) return false;
  const walk = (n: LayoutNode): boolean =>
    n.kind === 'group' ? n.id === groupId : n.kind === 'split' && n.children.some(walk);
  return walk(col);
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
  let found: SplitNode | undefined;
  const walk = (n: LayoutNode) => {
    if (n.kind !== 'split') return;
    if (n.dir === 'column' && n.children.some((c) => c.kind === 'group' && c.id === groupId)) found = n;
    for (const c of n.children) walk(c);
  };
  walk(layout.root);
  return found;
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
 * A stored layout back as a layout, or undefined when it is not one (another version, a wrong shape): checked field by
 * field, then normalized (shell §2, persistence). Unknown view refs are kept: a view that registers later lands where
 * it was left.
 */
export function validateLayout(value: unknown): ShellLayout | undefined {
  if (!isObject(value) || value.version !== LayoutVersion) return undefined;
  const root = validNode(value.root, 0);
  if (!root) return undefined;
  const windows: WindowNode[] = [];
  if (Array.isArray(value.windows)) {
    for (const w of value.windows) {
      if (!isObject(w)) continue;
      const g = validNode(w.group, 0);
      const r = w.rect;
      if (!g || g.kind !== 'group' || !isObject(r)) continue;
      const nums = [r.x, r.y, r.width, r.height];
      if (!nums.every((n) => typeof n === 'number' && Number.isFinite(n))) continue;
      const z = typeof w.z === 'number' && Number.isFinite(w.z) ? w.z : windows.length + 1;
      windows.push({ group: g, rect: r as unknown as Rect, z });
    }
  }
  const hidden = Array.isArray(value.hidden) ? value.hidden.filter((v): v is string => typeof v === 'string') : [];
  const seq = typeof value.seq === 'number' && Number.isFinite(value.seq) ? Math.max(1, Math.floor(value.seq)) : 1;
  return normalize({ version: LayoutVersion, root, windows, hidden, rail: value.rail !== false, seq });
}

function validNode(v: unknown, depth: number): LayoutNode | undefined {
  if (!isObject(v) || depth > 16) return undefined;
  const size = typeof v.size === 'number' && Number.isFinite(v.size) ? Math.max(MinGroupSize, v.size) : undefined;
  if (v.kind === 'main') return { kind: 'main', id: 'main' };
  if (v.kind === 'group' && typeof v.id === 'string' && Array.isArray(v.views)) {
    const views = v.views.filter((x): x is string => typeof x === 'string');
    return clean({
      kind: 'group' as const,
      id: v.id,
      views,
      active: typeof v.active === 'string' ? v.active : undefined,
      collapsed: v.collapsed === true || undefined,
      size,
    });
  }
  if (
    v.kind === 'split' &&
    typeof v.id === 'string' &&
    (v.dir === 'row' || v.dir === 'column') &&
    Array.isArray(v.children)
  ) {
    const children = v.children.map((c) => validNode(c, depth + 1)).filter((c): c is LayoutNode => !!c);
    return clean({ kind: 'split' as const, id: v.id, dir: v.dir as 'row' | 'column', children, size });
  }
  return undefined;
}

function isObject(v: unknown): v is Record<string, unknown> {
  return typeof v === 'object' && v !== null && !Array.isArray(v);
}

// ---------------------------------------------------------------------------------------------------------------------
// Window placement

/**
 * A window's rectangle kept inside the area it floats over (shell §2): no larger than the area, and placed so that all
 * of it is in when it fits; its header (its top `keep` px) is always in.
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
