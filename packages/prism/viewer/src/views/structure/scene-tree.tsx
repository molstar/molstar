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
import {
  Button,
  colorToHex,
  keepRows,
  MenuPopover,
  RecycledList,
  VirtualList,
  type VirtualListHandle,
  type VirtualSlot,
} from '@molstar/prism-ui';
import { Icon } from '@molstar/prism-ui/components/icon/icon';
import { uniqueId } from '@molstar/prism-ui/dom/focus';
import { Show } from '@molstar/prism-ui/dom/control-flow';
import { bind, effect, onActive, onCleanup, requireScope, type WatchSource, watch } from '@molstar/prism-ui/dom/scope';
import { workQueue } from '@molstar/prism-ui/dom/scoped';
import { setShown } from '@molstar/prism-ui/dom/writes';
import { setAttr, setClass, setStyle, setText } from '@molstar/prism-ui/dom/writes';
import { TreeNavigation } from '@molstar/prism-ui/dom/tree-keys';
import type { ShellModel, ShellSection } from '@molstar/prism-ui/shell/shell';
import type { SceneMenuItem, SceneModel, SceneNode } from '../../model/scene.js';
import { batch, observableState } from '@molstar/prism-ui/reactive/observable-state';
import { type KeyedSelector, keyedSelector } from '@molstar/prism-ui/reactive/keyed-selector';
import { type KeyedState, keyedState } from '@molstar/prism-ui/reactive/keyed-state';
import { classes } from './_generated/scene-tree.module.css.js';

/**
 * Where extensions' rows go (viewer-spec.md §3, `rows`): the extension host registers each as a shell section in this
 * slot, its `when(node)` as the section's `when` on the node; the tree shows them under the nodes they are for.
 */
export { SceneRowSlot } from './scene-row-slot.js';
import { SceneRowSlot } from './scene-row-slot.js';

import { derive } from '@molstar/prism-ui/reactive/state';
import { cancelFrameCallback, requestFrameCallback } from '@molstar/prism-ui/reactive/frame-state';
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
export function SceneTree(props: {
  scene: SceneModel;
  shell?: ShellModel;
  filter?: Observable<string>;
  /**
   * Node rows as recycled slots (`RecycledList`, the default) or keyed rows (`VirtualList`, false: the measurement's
   * baseline; rendering.md §1).
   */
  recycle?: boolean;
}) {
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
  // The rows' places and open states, by row key: a row listens to its own (view-models.md §3), so a reflow wakes only
  // the rows whose place changed; and the picked node, told to the two rows a pick concerns.
  const states = keyedState<string, RowState>(samePlace);
  const picked = keyedSelector(scene.picked);
  const keep = keepRows<Row>((r) => r.key, sameRow);
  // Long lists shown whole (Show all), until closed; the nodes open at the last flattening (a node opening loads).
  const showAll = new Set<SceneNode>();
  const wasOpen = new Set<SceneNode>();
  // A fixture without the scene's change stream (`SceneModel.changes`): the shown nodes' `isExpanded` and `children`
  // followed one by one, as before it.
  const watched = new Map<SceneNode, Subscription>();
  let current: Search | undefined;

  /** What a flattening makes: rows, their places, the nodes it found open and the nodes it showed. */
  interface Flat {
    readonly out: Row[];
    readonly st: Map<string, RowState>;
    readonly open: Set<SceneNode>;
  }
  // Flattens `row` (at its place among its siblings) and what is open under it, into `f`.
  const emit = (s: Search, f: Flat, row: Row, posinset: number, setsize: number) => {
    f.out.push(row);
    if (row.kind !== 'node') {
      f.st.set(row.key, { posinset, setsize });
      return;
    }
    const { hits, text } = s;
    const { node, depth } = row;
    if (!scene.changes) followNode(node);
    // The extensions' rows for a node (searching: those whose title matches).
    const extensions = shell
      ? s.sections.filter(
          (r) => (!r.when || r.when({ shell, subject: node })) && (!text || r.title.toLowerCase().includes(text)),
        )
      : [];
    let children = node.children.value;
    // A node with an extension's row can be opened before it has children.
    const expandable = children.length + extensions.length > 0 || !!node.lazy;
    // Searching, a row with matches under it is open; otherwise as the person left it.
    const isOpen =
      expandable && (node.isExpanded.value || (!!hits && (children.some((c) => hits.has(c)) || extensions.length > 0)));
    f.st.set(row.key, { posinset, setsize, expanded: expandable ? isOpen : undefined });
    if (!isOpen) return;
    f.open.add(node);
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
    for (let i = 0; i < next.length; i++) emit(s, f, next[i], i + 1, next.length);
  };
  /** The nodes found open now, where they were open before (`was`): Show all is forgotten for a closed one. */
  const opened = (was: Iterable<SceneNode>, open: Set<SceneNode>) => {
    for (const n of was) {
      if (open.has(n)) continue;
      showAll.delete(n);
      wasOpen.delete(n);
    }
    for (const n of open) wasOpen.add(n);
  };

  // Every row (a search, the roots or the extensions' rows changed).
  const flatten = (s: Search) => {
    const f: Flat = { out: [], st: new Map(), open: new Set() };
    const top = s.hits ? s.roots.filter((r) => s.hits!.has(r)) : s.roots;
    for (let i = 0; i < top.length; i++)
      emit(s, f, { kind: 'node', key: top[i].id, node: top[i], depth: 0 }, i + 1, top.length);
    opened([...wasOpen], f.open);
    if (!scene.changes) {
      const seen = new Set<SceneNode>();
      for (const r of f.out) if (r.kind === 'node') seen.add(r.node);
      for (const [n, sub] of watched) {
        if (seen.has(n)) continue;
        sub.unsubscribe();
        watched.delete(n);
      }
    }
    batch(() => {
      states.replace(f.st);
      rows.set(keep(f.out));
    });
  };

  // Only what changed (view-models.md §3, derived views update incrementally): each changed node's subtree flattened
  // again and put in place of the old one, the rest of the rows kept as they are.
  const reflowNodes = (s: Search, ids: ReadonlySet<string>) => {
    const before = rows.value;
    // The changed nodes' rows, in order, without those inside another's subtree (flattened with it).
    const spans: [number, number][] = [];
    for (let i = 0; i < before.length; i++) {
      const r = before[i];
      if (r.kind !== 'node' || !ids.has(r.node.id)) continue;
      let end = i + 1;
      while (end < before.length && before[end].depth > r.depth) end++;
      spans.push([i, end]);
      i = end - 1;
    }
    if (spans.length === 0) return;
    // One pass: the rows between the spans as they are, each span flattened again (a burst of a thousand changes
    // copies the rows once, not once a change).
    const next: Row[] = [];
    const st = new Map<string, RowState>();
    const gone = new Set<string>();
    let at = 0;
    for (const [i, end] of spans) {
      for (let j = at; j < i; j++) next.push(before[j]);
      const row = before[i];
      const place = states.get(row.key) ?? { posinset: 1, setsize: 1 };
      const from = next.length;
      const f: Flat = { out: next, st, open: new Set() };
      emit(s, f, row, place.posinset, place.setsize);
      const was: SceneNode[] = [];
      const old = new Map<string, Row>();
      for (let j = i; j < end; j++) {
        const r = before[j];
        if (r.kind === 'node' && wasOpen.has(r.node)) was.push(r.node);
        gone.add(r.key);
        old.set(r.key, r);
      }
      // A row the same as the span's before is that object again (its view kept), as `keep` does for all rows.
      for (let j = from; j < next.length; j++) {
        const o = old.get(next[j].key);
        if (o && sameRow(o, next[j])) next[j] = o;
      }
      opened(was, f.open);
      at = end;
    }
    for (let j = at; j < before.length; j++) next.push(before[j]);
    batch(() => {
      for (const key of gone) if (!st.has(key)) states.delete(key);
      for (const [key, place] of st) states.set(key, place);
      rows.set(next);
    });
  };

  // Flattened once for a burst of changes (rendering.md §2.7: a task's queue): expanding a thousand nodes at once, or a
  // load's nodes arriving, flattened every change, was quadratic (1,000 nodes: 400 ms; once: a few). Before the rows
  // render (the runtime's flush runs the queue first); a change while flattening (a lazy node's rows arriving as it
  // opens) flattens again in the same flush.
  const work = workQueue({ lane: 'task' });
  let full = true;
  let changed = new Set<string>();
  const reflowNow = () => {
    if (!current) return;
    const ids = changed;
    changed = new Set();
    // Searching, the hits decide what shows under a changed node: all of it again.
    if (full || current.hits) {
      full = false;
      flatten(current);
    } else reflowNodes(current, ids);
  };
  const reflow = () => {
    full = true;
    work.want('flatten', reflowNow);
  };
  const followNode = (node: SceneNode) => {
    if (watched.has(node)) return;
    // Kept in `watched`: unsubscribed when the node leaves the shown rows, and all of them when the tree is hidden.
    const sub = merge(node.isExpanded.pipe(skip(1)), node.children.pipe(skip(1))).subscribe(reflow);
    watched.set(node, sub);
  };
  // The first flattening now, with the mount (the rows are there when it returns, as a binding's first write is); the
  // changes after, once a burst.
  let mounted = false;
  bind(search, (s) => {
    current = s;
    if (mounted) reflow();
    else {
      mounted = true;
      full = false;
      flatten(s);
    }
  });
  // The scene's structure, one listener (`SceneModel.changes`): the changed nodes' subtrees flattened again.
  if (scene.changes) {
    let version = scene.changes.value.version;
    bind(scene.changes, (c) => {
      // Its current value (made, or resumed with nothing changed meanwhile): nothing to do.
      if (c.version === version) return;
      // Resumed after more than one change: what changed meanwhile is not known by node, so all of it.
      if (c.version !== version + 1) full = true;
      version = c.version;
      for (const id of c.nodes) changed.add(id);
      work.want('flatten', reflowNow);
    });
  }
  // Hidden, nothing followed (prism AGENTS.md): the bindings pause themselves, the nodes' watches go; shown again, the
  // search replays and the rows are flattened again.
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
  effect(scene.picked, (node) => {
    if (!node) return;
    cancelFrameCallback(frame);
    frame = requestFrameCallback(() => {
      const i = rows.value.findIndex((r) => r.kind === 'node' && r.node === node);
      if (i >= 0) list?.scrollToIndex(i);
    });
  });
  onCleanup(() => cancelFrameCallback(frame));
  // The rows that are not a node's: Show all, Loading…, an extension's row (keyed, never recycled).
  const otherRow = (row: Row) => {
    switch (row.kind) {
      case 'node':
        return null;
      case 'more':
        return (
          <MoreRow
            row={row}
            states={states}
            nav={nav}
            onShowAll={() => {
              showAll.add(row.node);
              changed.add(row.node.id);
              work.want('flatten', reflowNow);
            }}
          />
        );
      case 'loading':
        return <LoadingRow row={row} states={states} nav={nav} />;
      case 'row':
        return (
          <ExtensionRow
            section={row.section}
            shell={shell!}
            node={row.node}
            depth={row.depth}
            states={states}
            register={nav.row(row.key)}
          />
        );
    }
  };
  return (
    <>
      <Show when={derive(scene.roots, (r) => r.length === 0)}>
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
      {(props.recycle ?? true) ? (
        <RecycledList
          items={rows.asReadonly}
          itemKey={(r) => r.key}
          rowHeight={26}
          class={classes.tree}
          edge
          role="tree"
          aria-label="Scene"
          contentRef={nav.tree}
          ref={(l) => {
            list = l;
          }}
          slot={() => nodeSlot({ scene, nav, states, picked })}
          recycle={(row) => row.kind === 'node'}
        >
          {(row) => otherRow(row)}
        </RecycledList>
      ) : (
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
          {(row) =>
            row.kind === 'node' ? (
              <TreeItem scene={scene} node={row.node} depth={row.depth} states={states} picked={picked} nav={nav} />
            ) : (
              otherRow(row)
            )
          }
        </VirtualList>
      )}
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

/** The rows' places and open states, by key: a row follows its own (`states.key(key)`). */
type RowStates = KeyedState<string, RowState>;

const samePlace = (a: RowState, b: RowState) =>
  a.posinset === b.posinset && a.setsize === b.setsize && a.expanded === b.expanded;

/** A row's place among its siblings (and its open state, if it opens), written by its render. */
function writePlace(el: Element, s: RowState | undefined) {
  setAttr(el, 'aria-posinset', s?.posinset);
  setAttr(el, 'aria-setsize', s?.setsize);
}

/**
 * A node's row as view = f(state) (rendering.md §2): its DOM made once, with refs to the parts that change, and one
 * `watch` render that writes what differs from the node's and the tree's current values: label, description, swatch,
 * hidden, picked, its place and open state, its menu's. What is fixed for the row (its node, depth, icon, whether it has
 * a swatch or a menu) is built in; a row whose node or depth changes is made anew (`For` by item).
 */
function TreeItem(props: {
  scene: SceneModel;
  node: SceneNode;
  depth: number;
  states: RowStates;
  picked: KeyedSelector<SceneNode | undefined>;
  nav: TreeNavigation<Row>;
}) {
  const { scene, node, depth, states } = props;
  const menuOpen = observableState(false);
  const menu = hasMenu(node);
  // Named by its own label and description.
  const labelId = uniqueId('msp-tree-label');
  const label = document.createTextNode('');
  const description = document.createTextNode('');
  let line!: HTMLElement;
  let expander!: HTMLElement;
  let eyeButton!: HTMLElement;
  let swatch: HTMLElement | undefined;
  let menuButton: HTMLElement | undefined;
  // The eye: one icon shown, the other made the first time it is needed (most rows are never hidden).
  const eye = Icon({ name: 'eye', size: 14 });
  let eyeOff: SVGSVGElement | undefined;
  let eyeShown: SVGSVGElement = eye;
  const el = (
    <div
      role="treeitem"
      aria-level={depth + 1}
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
        style={{ '--msp-tree-depth': depth }}
        ref={(e: HTMLElement) => (line = e)}
        onClick={() => scene.pick(node)}
        onPointerEnter={() => node.highlight(true)}
        onPointerLeave={() => node.highlight(false)}
      >
        <button
          type="button"
          class={classes.expander}
          data-tree-toggle=""
          tabIndex={-1}
          aria-label="Expand"
          ref={(e: HTMLElement) => (expander = e)}
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
          {label}
        </span>
        {node.swatch ? <span class={classes.swatch} aria-hidden="true" ref={(e: HTMLElement) => (swatch = e)} /> : null}
        <span class={classes.description} id={`${labelId}-d`} data-tip-part="">
          {description}
        </span>
        <span class={classes.actions}>
          <button
            type="button"
            class={classes.action}
            ref={(e: HTMLElement) => (eyeButton = e)}
            onClick={(e) => {
              e.stopPropagation();
              node.toggleVisibility();
            }}
          >
            {eye}
          </button>
          {menu ? (
            <button
              type="button"
              class={classes.action}
              title="More"
              aria-haspopup="menu"
              ref={(e: HTMLElement) => (menuButton = e)}
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
  ) as HTMLElement;
  // Its own keys of the shared states: a pick tells two rows, a reflow the rows whose place changed.
  const place = states.key(node.id);
  const isPicked = props.picked.key(node);
  const sources: WatchSource[] = [node.label, node.description, node.isHidden, isPicked, place];
  if (node.swatch) sources.push(node.swatch);
  if (menu) sources.push(menuOpen);
  // What the last render wrote, by source: a render that runs for one source's change compares the others' values here,
  // in JS, and touches the DOM only for what changed (a pick renders every row: most write nothing).
  const NotYet = {};
  let lastLabel: unknown = NotYet;
  let lastDescription: unknown = NotYet;
  let lastSwatch: unknown = NotYet;
  let lastHidden: unknown = NotYet;
  let lastPicked: unknown = NotYet;
  let lastPlace: RowState | undefined | typeof NotYet = NotYet;
  let lastOpen: unknown = NotYet;
  watch(sources, () => {
    const text = node.label.value;
    if (text !== lastLabel) setText(label, (lastLabel = text));
    const d = node.description.value;
    if (d !== lastDescription) setText(description, (lastDescription = d) ?? '');
    if (swatch) {
      const c = node.swatch!.value;
      if (c !== lastSwatch) setStyle(swatch, 'background', (lastSwatch = c) === undefined ? '' : colorToHex(c!));
    }
    const hidden = node.isHidden.value;
    if (hidden !== lastHidden) {
      lastHidden = hidden;
      setClass(line, classes.hidden, hidden);
      setClass(eyeButton, classes.keep, hidden);
      setAttr(eyeButton, 'title', hidden ? 'Show' : 'Hide');
      const icon = hidden ? (eyeOff ??= Icon({ name: 'eye-off', size: 14 })) : eye;
      if (icon !== eyeShown) {
        eyeShown.replaceWith(icon);
        eyeShown = icon;
      }
    }
    const picked = isPicked.value;
    if (picked !== lastPicked) {
      lastPicked = picked;
      setAttr(el, 'aria-selected', picked ? 'true' : 'false');
      setClass(line, classes.picked, picked);
    }
    const s = place.value;
    const was = lastPlace as RowState | undefined;
    if (
      lastPlace === NotYet ||
      s?.posinset !== was?.posinset ||
      s?.setsize !== was?.setsize ||
      s?.expanded !== was?.expanded
    ) {
      lastPlace = s;
      writePlace(el, s);
      setAttr(el, 'aria-expanded', s?.expanded === undefined ? undefined : s.expanded ? 'true' : 'false');
      setClass(expander, classes.leaf, s?.expanded === undefined);
    }
    if (menuButton) {
      const open = menuOpen.value;
      if (open !== lastOpen) {
        lastOpen = open;
        setClass(menuButton, classes.keep, open);
        setAttr(menuButton, 'aria-expanded', open ? 'true' : 'false');
      }
    }
  });
  return el;
}

/**
 * A node's row, recycled (`RecycledList`, rendering.md §1): its DOM made once per slot, with every part a node may have
 * (the top-level icon, the swatch, the menu button) shown or hidden per node, and one `watch` re-pointed at the node it
 * shows (`retarget`): its own states and its keys of the tree's (its place, whether it is picked). Re-pointed, it
 * closes its menu, moves the hover highlight, and the keyboard's registration (`TreeNavigation.slot`).
 */
function nodeSlot(props: {
  scene: SceneModel;
  nav: TreeNavigation<Row>;
  states: RowStates;
  picked: KeyedSelector<SceneNode | undefined>;
}): VirtualSlot<Row> {
  const { scene, states } = props;
  let node: SceneNode | undefined;
  let depth = 0;
  let hovered = false;
  const menuOpen = observableState(false);
  const labelId = uniqueId('msp-tree-label');
  const label = document.createTextNode('');
  const description = document.createTextNode('');
  let line!: HTMLElement;
  let expander!: HTMLElement;
  let labelEl!: HTMLElement;
  let swatch!: HTMLElement;
  let eyeButton!: HTMLElement;
  let menuButton!: HTMLElement;
  let actionSlot!: HTMLElement;
  const eye = Icon({ name: 'eye', size: 14 });
  let eyeOff: SVGSVGElement | undefined;
  let eyeShown: SVGSVGElement = eye;
  const el = (
    <div role="treeitem" aria-labelledby={`${labelId} ${labelId}-d`} tabIndex={-1} class={classes.item}>
      <div
        data-tip-row=""
        data-tree-row=""
        class={classes.row}
        ref={(e: HTMLElement) => (line = e)}
        onClick={() => node && scene.pick(node)}
        onPointerEnter={() => {
          hovered = true;
          node?.highlight(true);
        }}
        onPointerLeave={() => {
          hovered = false;
          node?.highlight(false);
        }}
      >
        <button
          type="button"
          class={classes.expander}
          data-tree-toggle=""
          tabIndex={-1}
          aria-label="Expand"
          ref={(e: HTMLElement) => (expander = e)}
          onClick={(e) => {
            e.stopPropagation();
            node?.toggleExpanded();
          }}
        >
          <Icon name="chevron-right" size={14} />
        </button>
        <span class={classes.label} id={labelId} data-tip-part="" ref={(e: HTMLElement) => (labelEl = e)}>
          {label}
        </span>
        <span class={classes.swatch} aria-hidden="true" ref={(e: HTMLElement) => (swatch = e)} />
        <span class={classes.description} id={`${labelId}-d`} data-tip-part="">
          {description}
        </span>
        <span class={classes.actions}>
          <button
            type="button"
            class={classes.action}
            ref={(e: HTMLElement) => (eyeButton = e)}
            onClick={(e) => {
              e.stopPropagation();
              node?.toggleVisibility();
            }}
          >
            {eye}
          </button>
          <button
            type="button"
            class={classes.action}
            title="More"
            aria-haspopup="menu"
            ref={(e: HTMLElement) => (menuButton = e)}
            onClick={(e) => {
              e.stopPropagation();
              menuOpen.set(!menuOpen.value);
            }}
          >
            <Icon name="ellipsis" size={14} />
          </button>
          {/* No menu: its slot stays, so every row's eye sits in the same column. */}
          <span class={classes.actionSlot} aria-hidden="true" ref={(e: HTMLElement) => (actionSlot = e)} />
        </span>
      </div>
      <Show when={menuOpen}>
        {() => (
          <RowMenu
            items={node ? menuOf(node) : []}
            anchor={menuButton}
            onDone={(refocus) => {
              menuOpen.set(false);
              if (refocus) menuButton.focus();
            }}
          />
        )}
      </Show>
    </div>
  ) as HTMLElement;
  const point = props.nav.slot(el);
  // The top-level row's icon: one element per icon the slot has shown, made once, shown or not.
  const kinds = new Map<string, SVGSVGElement>();
  let kind: SVGSVGElement | undefined;
  let place = states.key('');
  let isPicked = props.picked.key(undefined);
  const NotYet = {};
  let lastLabel: unknown = NotYet;
  let lastDescription: unknown = NotYet;
  let lastSwatch: unknown = NotYet;
  let lastHidden: unknown = NotYet;
  let lastPicked: unknown = NotYet;
  let lastPlace: RowState | undefined | typeof NotYet = NotYet;
  let lastOpen: unknown = NotYet;
  let lastNode: SceneNode | undefined;
  const render = () => {
    // Its menu's state first: a slot put away (no node) still shows its menu closed.
    const open = menuOpen.value;
    if (open !== lastOpen) {
      lastOpen = open;
      setClass(menuButton, classes.keep, open);
      setAttr(menuButton, 'aria-expanded', open ? 'true' : 'false');
    }
    const n = node;
    if (!n) return;
    if (n !== lastNode) {
      lastNode = n;
      setAttr(el, 'aria-level', depth + 1);
      setStyle(line, '--msp-tree-depth', String(depth));
      setAttr(line, 'data-node-id', n.id);
      setShown(swatch, !!n.swatch);
      const menu = hasMenu(n);
      setShown(menuButton, menu);
      setShown(actionSlot, !menu);
      const icon = depth === 0 ? n.icon : undefined;
      let next = icon ? kinds.get(icon) : undefined;
      if (icon && !next) {
        next = Icon({ name: icon, size: 15, class: classes.kind });
        kinds.set(icon, next);
        line.insertBefore(next, labelEl);
      }
      if (next !== kind) {
        if (kind) kind.style.display = 'none';
        if (next) next.style.display = '';
        kind = next;
      }
    }
    const text = n.label.value;
    if (text !== lastLabel) setText(label, (lastLabel = text));
    const d = n.description.value;
    if (d !== lastDescription) setText(description, (lastDescription = d) ?? '');
    const c = n.swatch?.value;
    if (c !== lastSwatch) setStyle(swatch, 'background', (lastSwatch = c) === undefined ? '' : colorToHex(c!));
    const hidden = n.isHidden.value;
    if (hidden !== lastHidden) {
      lastHidden = hidden;
      setClass(line, classes.hidden, hidden);
      setClass(eyeButton, classes.keep, hidden);
      setAttr(eyeButton, 'title', hidden ? 'Show' : 'Hide');
      const icon = hidden ? (eyeOff ??= Icon({ name: 'eye-off', size: 14 })) : eye;
      if (icon !== eyeShown) {
        eyeShown.replaceWith(icon);
        eyeShown = icon;
      }
    }
    const picked = isPicked.value;
    if (picked !== lastPicked) {
      lastPicked = picked;
      setAttr(el, 'aria-selected', picked ? 'true' : 'false');
      setClass(line, classes.picked, picked);
    }
    const s = place.value;
    const was = lastPlace as RowState | undefined;
    if (
      lastPlace === NotYet ||
      s?.posinset !== was?.posinset ||
      s?.setsize !== was?.setsize ||
      s?.expanded !== was?.expanded
    ) {
      lastPlace = s;
      writePlace(el, s);
      setAttr(el, 'aria-expanded', s?.expanded === undefined ? undefined : s.expanded ? 'true' : 'false');
      setClass(expander, classes.leaf, s?.expanded === undefined);
    }
  };
  const handle = watch([], render);
  return {
    el,
    show(next) {
      const n = next?.kind === 'node' ? next.node : undefined;
      if (n === node) {
        if (next?.kind === 'node' && next.depth !== depth) {
          depth = next.depth;
          lastNode = undefined;
        }
        render();
        return;
      }
      // Its menu closes on the node it was open for; the hover moves to the node it shows now.
      if (menuOpen.value) menuOpen.set(false);
      if (hovered) node?.highlight(false);
      node = n;
      depth = next?.kind === 'node' ? next.depth : 0;
      point(n?.id);
      if (!n) hovered = false;
      else if (hovered) n.highlight(true);
      if (!n) {
        handle.retarget([]);
        return;
      }
      place = states.key(n.id);
      isPicked = props.picked.key(n);
      const sources: WatchSource[] = [n.label, n.description, n.isHidden, isPicked, place, menuOpen];
      if (n.swatch) sources.push(n.swatch);
      handle.retarget(sources);
    },
  };
}

/** Loading…: a quiet row while a lazy node computes its rows. */
function LoadingRow(props: { row: Extract<Row, { kind: 'loading' }>; states: RowStates; nav: TreeNavigation<Row> }) {
  const { row, states } = props;
  const el = (
    <div
      role="treeitem"
      aria-level={row.depth + 1}
      aria-busy="true"
      tabIndex={-1}
      class={classes.item}
      ref={props.nav.row(row.key)}
    >
      <div class={classes.loading} style={{ '--msp-tree-depth': row.depth }}>
        Loading…
      </div>
    </div>
  ) as HTMLElement;
  const place = states.key(row.key);
  watch([place], () => writePlace(el, place.value));
  return el;
}

/** Show all: a quiet row where the rest of a long list would be, at the children's depth. */
function MoreRow(props: {
  row: Extract<Row, { kind: 'more' }>;
  states: RowStates;
  nav: TreeNavigation<Row>;
  onShowAll: () => void;
}) {
  const { row, states } = props;
  const el = (
    <div role="treeitem" aria-level={row.depth + 1} tabIndex={-1} class={classes.item} ref={props.nav.row(row.key)}>
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
  ) as HTMLElement;
  const place = states.key(row.key);
  watch([place], () => writePlace(el, place.value));
  return el;
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
  states: RowStates;
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
  const key = `${props.node.id}::row:${props.section.id}`;
  const place = props.states.key(key);
  watch([place], () => {
    const el = item();
    const s = place.value;
    if (el && s) writePlace(el, s);
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

/** Whether a row has a "…" menu, without building it (rows are made as the tree scrolls). */
function hasMenu(node: SceneNode): boolean {
  return node.removable !== false || (node.hasMenu ?? (node.menu?.length ?? 0) > 0);
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
