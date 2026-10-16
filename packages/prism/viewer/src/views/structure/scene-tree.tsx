/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import {
  combineLatest,
  distinctUntilChanged,
  map,
  merge,
  type Observable,
  of,
  shareReplay,
  skip,
  startWith,
  type Subscription,
} from 'rxjs';
import { Button, colorToHex, keepRows, MenuPopover, VirtualList, type VirtualListHandle } from '@molstar/prism-ui';
import { Icon } from '@molstar/prism-ui/components/icon/icon';
import { uniqueId } from '@molstar/prism-ui/dom/focus';
import { Show } from '@molstar/prism-ui/dom/control-flow';
import { bind, onActive, onCleanup, requireScope } from '@molstar/prism-ui/dom/scope';
import { TreeNavigation } from '@molstar/prism-ui/dom/tree-keys';
import type { ShellModel, ShellSection } from '@molstar/prism-ui/shell/shell';
import type { SceneMenuItem, SceneModel, SceneNode } from '../../model/scene.js';
import { batch, observableState } from '@molstar/prism-ui/reactive/observable-state';
import { classes } from './_generated/scene-tree.module.css.js';

/**
 * Where extensions' rows go (viewer-spec.md §3, `rows`): the extension host registers each as a shell section in this
 * slot, its `when(node)` as the section's `when` on the node; the tree shows them under the nodes they are for.
 */
export { SceneRowSlot } from './scene-row-slot.js';
import { SceneRowSlot } from './scene-row-slot.js';

/**
 * One line of the tree, flattened (viewer-spec.md §3, the scene tree): a node; Show all, for the rest of a long list;
 * Loading…, while a lazy node computes its rows; an extension's row under a node. All are rows of one height.
 */
type Row =
  | { readonly kind: 'node'; readonly key: string; readonly node: SceneNode; readonly depth: number }
  | {
      readonly kind: 'more';
      readonly key: string;
      readonly node: SceneNode;
      readonly depth: number;
      /** The node's children: Show all says how many. */
      readonly total: number;
    }
  | { readonly kind: 'loading'; readonly key: string; readonly node: SceneNode; readonly depth: number }
  | {
      readonly kind: 'row';
      readonly key: string;
      readonly node: SceneNode;
      readonly depth: number;
      readonly section: ShellSection;
    };

/** A row's place among its siblings, and whether it is open (undefined: it has nothing to open). */
interface RowState {
  readonly posinset: number;
  readonly setsize: number;
  readonly expanded?: boolean;
}

/** What the rows are flattened from: the search, its hits, the scene's roots, the extensions' rows. */
interface Search {
  readonly text: string;
  readonly hits: SearchHits;
  readonly roots: readonly SceneNode[];
  readonly sections: readonly ShellSection[];
}

const sameRow = (a: Row, b: Row) =>
  a.kind === b.kind &&
  a.node === b.node &&
  a.depth === b.depth &&
  (a.kind !== 'more' || a.total === (b as typeof a).total) &&
  (a.kind !== 'row' || a.section === (b as typeof a).section);

/**
 * The scene as one tree: structures → components → representations, volumes, measurements. Rows are keyed by node, so
 * a hierarchy update touches only the rows that changed; hover highlights in 3D, a click picks the node for the
 * inspector, visibility toggles inline. Given the shell, extensions' rows (`SceneRowSlot`) follow a node's children.
 *
 * Virtualized (view-models.md §3): the shown rows are flattened (a node's children while it is open) into rows of one
 * height in a `VirtualList`, so only the window is rendered; the flattening follows each shown node's `isExpanded` and
 * `children` (nothing else), and a row follows its node's own states. The keyboard moves over the flattened rows
 * (`TreeNavigation`), scrolling its target into view.
 */
export function SceneTree(props: { scene: SceneModel; shell?: ShellModel; filter?: Observable<string> }) {
  const { scene, shell } = props;
  // A search: the rows that match, with the rows they are in (opened); lazy contents are looked into to be searched.
  const filter = (props.filter ?? of('')).pipe(
    // No search until the first (a debounced filter emits late).
    startWith(''),
    map((t) => t.trim().toLowerCase()),
    distinctUntilChanged(),
  );
  // One walk of the scene per search (and per change of the scene while searching), the rows flattened from its result.
  const sections: Observable<readonly ShellSection[]> = shell ? shell.sections(SceneRowSlot) : of([]);
  const search: Observable<Search> = combineLatest([filter, scene.roots, sections]).pipe(
    map(([text, roots, all]) => ({
      text,
      roots,
      sections: all,
      hits: text
        ? searchHits(roots, text, (node) => !!shell && all.some((r) => rowMatches(shell, r, node, text)))
        : undefined,
    })),
    shareReplay({ bufferSize: 1, refCount: true }),
  );

  const rows = observableState<readonly Row[]>([]);
  const states = observableState<ReadonlyMap<string, RowState>>(new Map());
  const keep = keepRows<Row>((r) => r.key, sameRow);
  // Long lists shown whole (Show all), until closed; the nodes open at the last flattening (a node opening loads).
  const showAll = new Set<SceneNode>();
  let wasOpen = new Set<SceneNode>();
  // The shown nodes' `isExpanded` and `children`: a change flattens again.
  const watched = new Map<SceneNode, Subscription>();
  let current: Search | undefined;
  let building = false;
  let again = false;

  const flatten = (s: Search) => {
    const out: Row[] = [];
    const st = new Map<string, RowState>();
    const open = new Set<SceneNode>();
    const seen = new Set<SceneNode>();
    const { hits, text } = s;
    // The extensions' rows for a node (searching: those whose title matches).
    const extensionRows = (node: SceneNode) =>
      shell
        ? s.sections.filter(
            (r) => (!r.when || r.when({ shell, subject: node })) && (!text || r.title.toLowerCase().includes(text)),
          )
        : [];
    const level = (items: readonly Row[]) => {
      for (let i = 0; i < items.length; i++) {
        const row = items[i];
        out.push(row);
        if (row.kind !== 'node') {
          st.set(row.key, { posinset: i + 1, setsize: items.length });
          continue;
        }
        const { node, depth } = row;
        seen.add(node);
        watch(node);
        const extensions = extensionRows(node);
        let children = node.children.value;
        // A node with an extension's row can be opened before it has children.
        const expandable = children.length + extensions.length > 0 || !!node.lazy;
        // Searching, a row with matches under it is open; otherwise as the person left it.
        const isOpen =
          expandable &&
          (node.isExpanded.value || (!!hits && (children.some((c) => hits.has(c)) || extensions.length > 0)));
        st.set(row.key, { posinset: i + 1, setsize: items.length, expanded: expandable ? isOpen : undefined });
        if (!isOpen) continue;
        open.add(node);
        if (!wasOpen.has(node)) {
          // Children computed on first look (lazy nodes); open only for a search, a look (it stays lazy).
          node.load?.(!node.isExpanded.value);
          children = node.children.value;
        }
        const shown = hits ? children.filter((c) => hits.has(c)) : children;
        // Many children: the first few, then Show all (until the row is closed: opened again, the few again).
        const paged = !hits && Paged.has(node.kind) && !showAll.has(node) && children.length > PageSize + 1;
        const next: Row[] = (paged ? shown.slice(0, PageSize) : shown).map((c) => ({
          kind: 'node',
          key: c.id,
          node: c,
          depth: depth + 1,
        }));
        if (paged) next.push({ kind: 'more', key: `${node.id}::more`, node, depth: depth + 1, total: children.length });
        // Until they arrive, a quiet loading row (not when only extensions' rows are there: nothing is loading).
        if (children.length === 0 && (node.lazy || extensions.length === 0))
          next.push({ kind: 'loading', key: `${node.id}::loading`, node, depth: depth + 1 });
        for (const section of extensions)
          next.push({ kind: 'row', key: `${node.id}::row:${section.id}`, node, depth: depth + 1, section });
        level(next);
      }
    };
    level(
      (hits ? s.roots.filter((r) => hits.has(r)) : s.roots).map((n) => ({
        kind: 'node',
        key: n.id,
        node: n,
        depth: 0,
      })),
    );
    for (const n of wasOpen) if (!open.has(n)) showAll.delete(n);
    wasOpen = open;
    for (const [n, sub] of watched) {
      if (seen.has(n)) continue;
      sub.unsubscribe();
      watched.delete(n);
    }
    batch(() => {
      states.set(st);
      rows.set(keep(out));
    });
  };
  // A change while flattening (a lazy node's rows arriving as it opens) flattens again after.
  const reflow = () => {
    if (!current) return;
    if (building) {
      again = true;
      return;
    }
    building = true;
    try {
      do {
        again = false;
        flatten(current);
      } while (again);
    } finally {
      building = false;
    }
  };
  const watch = (node: SceneNode) => {
    if (watched.has(node)) return;
    watched.set(node, merge(node.isExpanded.pipe(skip(1)), node.children.pipe(skip(1))).subscribe(reflow));
  };
  bind(search, (s) => {
    current = s;
    reflow();
  });
  // Hidden, nothing followed (prism AGENTS.md): the search's binding pauses itself, the nodes' watches go; shown again,
  // the search replays and the flattening watches them again.
  const unwatch = () => {
    for (const sub of watched.values()) sub.unsubscribe();
    watched.clear();
  };
  onActive((active) => {
    if (!active) unwatch();
  });
  onCleanup(unwatch);

  let list: VirtualListHandle | undefined;
  const nav = new TreeNavigation<Row>({
    rows: () => rows.value,
    key: (r) => r.key,
    level: (r) => r.depth + 1,
    scrollTo: (i) => list?.scrollToIndex(i),
    selected: () => scene.picked.value?.id,
  });
  // A row picked from elsewhere (the inspector's list, the canvas) is brought into view once its rows are open: the
  // frame after.
  let frame = 0;
  bind(scene.picked.pipe(distinctUntilChanged()), (node) => {
    if (!node) return;
    cancelAnimationFrame(frame);
    frame = requestAnimationFrame(() => {
      const i = rows.value.findIndex((r) => r.kind === 'node' && r.node === node);
      if (i >= 0) list?.scrollToIndex(i);
    });
  });
  onCleanup(() => cancelAnimationFrame(frame));
  const stateOf = (key: string) =>
    states.asReadonly.pipe(
      map((m) => m.get(key)),
      distinctUntilChanged(
        (a, b) => a?.posinset === b?.posinset && a?.setsize === b?.setsize && a?.expanded === b?.expanded,
      ),
      shareReplay({ bufferSize: 1, refCount: true }),
    );
  return (
    <>
      <Show
        when={scene.roots.pipe(
          map((r) => r.length === 0),
          distinctUntilChanged(),
        )}
      >
        {() => <p class={classes.empty}>Nothing loaded yet.</p>}
      </Show>
      <Show
        when={search.pipe(
          map((s) => !!s.hits && s.roots.length > 0 && !s.roots.some((r) => s.hits!.has(r))),
          distinctUntilChanged(),
        )}
      >
        {() => <p class={classes.empty}>Nothing matches.</p>}
      </Show>
      <VirtualList
        items={rows.asReadonly}
        key={(r) => r.key}
        rowHeight={26}
        class={classes.tree}
        edge
        role="tree"
        aria-label="Scene"
        contentRef={nav.tree}
        ref={(l) => {
          list = l;
        }}
      >
        {(row) => {
          const state = stateOf(row.key);
          switch (row.kind) {
            case 'node':
              return <TreeItem scene={scene} node={row.node} depth={row.depth} state={state} nav={nav} />;
            case 'more':
              return (
                <MoreRow
                  row={row}
                  state={state}
                  nav={nav}
                  onShowAll={() => {
                    showAll.add(row.node);
                    reflow();
                  }}
                />
              );
            case 'loading':
              return (
                <div
                  role="treeitem"
                  aria-level={row.depth + 1}
                  aria-posinset={state.pipe(map((s) => s?.posinset))}
                  aria-setsize={state.pipe(map((s) => s?.setsize))}
                  aria-busy="true"
                  tabIndex={-1}
                  class={classes.item}
                  ref={nav.row(row.key)}
                >
                  <div class={classes.loading} style={{ '--msp-tree-depth': row.depth }}>
                    Loading…
                  </div>
                </div>
              );
            case 'row':
              return (
                <ExtensionRow
                  section={row.section}
                  shell={shell!}
                  node={row.node}
                  depth={row.depth}
                  state={state}
                  register={nav.row(row.key)}
                />
              );
          }
        }}
      </VirtualList>
    </>
  );
}

/** Rows whose children show a few at first (a big structure's entities, an entity's chains): Show all for the rest. */
const Paged: ReadonlySet<SceneNode['kind']> = new Set(['contents', 'entity']);
const PageSize = 5;

/** The rows a search shows: those that match and the rows they are in; undefined: no search. */
type SearchHits = ReadonlySet<SceneNode> | undefined;

/** An extension's row for this node whose title matches a search (it shows under the node). */
const rowMatches = (shell: ShellModel, row: ShellSection, node: SceneNode, text: string) =>
  row.title.toLowerCase().includes(text) && (!row.when || row.when({ shell, subject: node }));

/**
 * The rows that match a search (their label or description, or an extension's row under them by its title) and the
 * rows they are in, in one walk; lazy rows are looked into (`load(true)`: they stay lazy, so later syncs do not
 * compute them).
 */
function searchHits(
  roots: readonly SceneNode[],
  text: string,
  hasRow: (node: SceneNode) => boolean,
): ReadonlySet<SceneNode> {
  const hits = new Set<SceneNode>();
  const hit = (v: string | undefined) => !!v && v.toLowerCase().includes(text);
  const walk = (node: SceneNode): boolean => {
    if (node.lazy) node.load?.(true);
    let any = hit(node.label.value) || hit(node.description.value) || hasRow(node);
    for (const c of node.children.value) if (walk(c)) any = true;
    if (any) hits.add(node);
    return any;
  };
  for (const r of roots) walk(r);
  return hits;
}

type RowStateStream = Observable<RowState | undefined>;

function TreeItem(props: {
  scene: SceneModel;
  node: SceneNode;
  depth: number;
  state: RowStateStream;
  nav: TreeNavigation<Row>;
}) {
  const { scene, node, depth, state } = props;
  const expanded = state.pipe(
    map((s) => s?.expanded),
    distinctUntilChanged(),
  );
  const menuOpen = observableState(false);
  let menuButton: HTMLElement | undefined;
  const picked = scene.picked.pipe(
    map((p) => p === node),
    distinctUntilChanged(),
  );
  // Named by its own label and description.
  const labelId = uniqueId('msp-tree-label');
  return (
    <div
      role="treeitem"
      aria-level={depth + 1}
      aria-posinset={state.pipe(map((s) => s?.posinset))}
      aria-setsize={state.pipe(map((s) => s?.setsize))}
      aria-expanded={expanded}
      aria-selected={picked}
      aria-labelledby={`${labelId} ${labelId}-d`}
      tabIndex={-1}
      class={classes.item}
      ref={props.nav.row(node.id)}
    >
      <div
        data-tip-row=""
        data-tree-row=""
        data-node-id={node.id}
        class={classes.row}
        classList={{ [classes.hidden]: node.isHidden, [classes.picked]: picked }}
        style={{ '--msp-tree-depth': depth }}
        onClick={() => scene.pick(node)}
        onPointerEnter={() => node.highlight(true)}
        onPointerLeave={() => node.highlight(false)}
      >
        <button
          type="button"
          class={classes.expander}
          classList={{ [classes.leaf]: expanded.pipe(map((e) => e === undefined)) }}
          data-tree-toggle=""
          tabIndex={-1}
          aria-label="Expand"
          onClick={(e) => {
            e.stopPropagation();
            node.toggleExpanded();
          }}
        >
          <Icon name="chevron-right" size={14} />
        </button>
        {/* Only top-level rows carry an icon (a structure, a volume, the measurements): below, the indent shows where a
            row belongs and its label and description say what it is. */}
        {depth === 0 ? <Icon name={node.icon} size={15} class={classes.kind} /> : null}
        <span class={classes.label} id={labelId} data-tip-part="">
          {node.label}
        </span>
        {node.swatch ? (
          <span
            class={classes.swatch}
            style={{ background: node.swatch.pipe(map((c) => (c === undefined ? '' : colorToHex(c)))) }}
            aria-hidden="true"
          />
        ) : null}
        <span class={classes.description} id={`${labelId}-d`} data-tip-part="">
          {node.description}
        </span>
        <span class={classes.actions}>
          <button
            type="button"
            class={classes.action}
            classList={{ [classes.keep]: node.isHidden }}
            title={node.isHidden.pipe(map((h) => (h ? 'Show' : 'Hide')))}
            onClick={(e) => {
              e.stopPropagation();
              node.toggleVisibility();
            }}
          >
            <Show when={node.isHidden} fallback={() => <Icon name="eye" size={14} />}>
              {() => <Icon name="eye-off" size={14} />}
            </Show>
          </button>
          {menuOf(node).length ? (
            <button
              type="button"
              class={classes.action}
              classList={{ [classes.keep]: menuOpen }}
              title="More"
              aria-haspopup="menu"
              aria-expanded={menuOpen}
              ref={(el: HTMLElement) => (menuButton = el)}
              onClick={(e) => {
                e.stopPropagation();
                menuOpen.set(!menuOpen.value);
              }}
            >
              <Icon name="ellipsis" size={14} />
            </button>
          ) : (
            // No menu: its slot stays, so every row's eye sits in the same column.
            <span class={classes.actionSlot} aria-hidden="true" />
          )}
        </span>
      </div>
      <Show when={menuOpen}>
        {() => (
          <RowMenu
            items={menuOf(node)}
            anchor={menuButton!}
            onDone={(refocus) => {
              menuOpen.set(false);
              if (refocus) menuButton?.focus();
            }}
          />
        )}
      </Show>
    </div>
  );
}

/** Show all: a quiet row where the rest of a long list would be, at the children's depth. */
function MoreRow(props: {
  row: Extract<Row, { kind: 'more' }>;
  state: RowStateStream;
  nav: TreeNavigation<Row>;
  onShowAll: () => void;
}) {
  const { row, state } = props;
  return (
    <div
      role="treeitem"
      aria-level={row.depth + 1}
      aria-posinset={state.pipe(map((s) => s?.posinset))}
      aria-setsize={state.pipe(map((s) => s?.setsize))}
      tabIndex={-1}
      class={classes.item}
      ref={props.nav.row(row.key)}
    >
      <button
        type="button"
        class={classes.more}
        style={{ '--msp-tree-depth': row.depth }}
        data-tree-action=""
        tabIndex={-1}
        onClick={props.onShowAll}
      >
        {`Show all ${row.total}`}
        <span class={classes.moreCount}>{`${row.total - PageSize} more`}</span>
      </button>
    </div>
  );
}

/**
 * An extension's row under a node: its view, mounted with the node as its subject, at the depth of the node's children
 * (`SceneRow` draws it as a tree row; the tree gives it its place: level, position, set size). In the tree it is one
 * row: its controls go to the host's inspector (`SceneRowPickEvent`), never under it.
 */
function ExtensionRow(props: {
  section: ShellSection;
  shell: ShellModel;
  node: SceneNode;
  depth: number;
  state: RowStateStream;
  register: (el: HTMLElement) => void;
}) {
  const host = (
    <div data-scene-row={props.section.id} style={{ display: 'contents', '--msp-tree-depth': props.depth }} />
  ) as HTMLElement;
  const dispose = props.section.render.mount(
    host,
    { shell: props.shell, subject: props.node },
    requireScope('A scene row'),
  );
  onCleanup(dispose);
  const item = () => host.querySelector<HTMLElement>('[role="treeitem"]');
  // Nothing opens under it here: not expandable.
  item()?.removeAttribute('aria-expanded');
  item()?.setAttribute('aria-level', String(props.depth + 1));
  bind(props.state, (s) => {
    const el = item();
    if (!el || !s) return;
    el.setAttribute('aria-posinset', String(s.posinset));
    el.setAttribute('aria-setsize', String(s.setsize));
  });
  props.register(host);
  return host;
}

/** A row's menu: its own entries, then Remove (set apart, in the danger color) when it can be removed. */
/**
 * A node's actions shown in place (the inspector of a node with nothing to set: its row menu's entries, not hidden
 * behind ⋯); sub-entries open in place, as in the row menu.
 */
export function NodeActions(props: { node: SceneNode }) {
  const items = menuOf(props.node);
  if (!items.length) return null;
  return (
    <div class={classes.menu} data-flat="" role="menu" aria-label={`${props.node.label.value}: actions`}>
      <MenuItems items={items} onDone={() => {}} />
    </div>
  );
}

export function menuOf(node: SceneNode): SceneMenuItem[] {
  const own = node.menu ?? [];
  if (node.removable === false) return [...own];
  return [...own, { id: 'remove', label: 'Remove', icon: 'trash-2', danger: true, run: () => node.remove() }];
}

/** A row's menu, a popover from its ⋯ button (`MenuPopover`): actions run and close it; entries with sub-entries open them in place. */
function RowMenu(props: { items: readonly SceneMenuItem[]; anchor: HTMLElement; onDone(refocus?: boolean): void }) {
  return (
    <MenuPopover anchor={props.anchor} onClose={props.onDone}>
      <div class={classes.menu} data-popover="">
        <MenuItems items={props.items} onDone={() => props.onDone(true)} />
      </div>
    </MenuPopover>
  );
}

function MenuItems(props: { items: readonly SceneMenuItem[]; onDone(): void }) {
  return (
    <>
      {props.items.map((item) =>
        item.items?.length && item.inline ? (
          <InlineItems item={item} onDone={props.onDone} />
        ) : item.items?.length ? (
          <SubMenu item={item} onDone={props.onDone} />
        ) : (
          <MenuItem item={item} onDone={props.onDone} />
        ),
      )}
    </>
  );
}

function MenuItem(props: { item: SceneMenuItem; onDone(): void }) {
  const { item } = props;
  return (
    <button
      type="button"
      role={item.checked !== undefined ? 'menuitemradio' : 'menuitem'}
      aria-checked={item.checked !== undefined ? item.checked : undefined}
      class={classes.menuItem}
      data-danger={item.danger ? '' : undefined}
      onClick={() => {
        item.run?.();
        props.onDone();
      }}
    >
      {/* The column a submenu's chevron takes: a choice's check mark. */}
      <span>{item.checked ? <Icon name="check" size={12} /> : null}</span>
      <span>{item.label}</span>
    </button>
  );
}

/** An entry whose choices are icons on its row (Modify by selection), always shown; each its label as the tooltip. */
function InlineItems(props: { item: SceneMenuItem; onDone(): void }) {
  const { item } = props;
  return (
    <div class={`${classes.menuItem} ${classes.inlineItem}`} role="group" aria-label={item.label}>
      <span />
      <span>{item.label}</span>
      <span class={classes.inlineIcons}>
        {item.items!.map((sub) => (
          <Button
            variant="ghost"
            size="xs"
            icon={<Icon name={sub.icon ?? 'ellipsis'} />}
            title={sub.label}
            {...(sub.disabled ? { disabled: sub.disabled } : {})}
            onClick={() => {
              sub.run?.();
              props.onDone();
            }}
          />
        ))}
      </span>
    </div>
  );
}

function SubMenu(props: { item: SceneMenuItem; onDone(): void }) {
  const { item } = props;
  const open = observableState(false);
  return (
    <div class={classes.subMenu}>
      <button
        type="button"
        role="menuitem"
        class={classes.menuItem}
        aria-haspopup="menu"
        aria-expanded={open}
        onClick={() => open.set(!open.value)}
      >
        <Icon name="chevron-right" size={12} class={classes.subChevron} />
        <span>{item.label}</span>
      </button>
      <Show when={open} reveal>
        {() => (
          <div class={classes.subItems} role="menu" aria-label={item.label}>
            <MenuItems items={item.items!} onDone={props.onDone} />
          </div>
        )}
      </Show>
    </div>
  );
}
