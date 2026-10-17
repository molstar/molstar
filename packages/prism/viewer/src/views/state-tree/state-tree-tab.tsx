/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { flattenTree, keepRows, MenuPopover, VirtualList, type VirtualListHandle } from '@molstar/prism-ui';
import { Icon } from '@molstar/prism-ui/components/icon/icon';
import { uniqueId } from '@molstar/prism-ui/dom/focus';
import { Show } from '@molstar/prism-ui/dom/control-flow';
import { bind, watch } from '@molstar/prism-ui/dom/scope';
import { workQueue } from '@molstar/prism-ui/dom/scoped';
import { setAttr, setClass } from '@molstar/prism-ui/dom/writes';
import { TreeNavigation } from '@molstar/prism-ui/dom/tree-keys';
import { ScrollArea } from '@molstar/prism-ui';
import type { StateTreeInspection, StateTreeModel, StateTreeNode } from '../../model/state-tree.js';
import { batch, observableState } from '@molstar/prism-ui/reactive/observable-state';
import { type KeyedSelector, keyedSelector } from '@molstar/prism-ui/reactive/keyed-selector';
import { type KeyedState, keyedState } from '@molstar/prism-ui/reactive/keyed-state';
import { ActionBlock } from '../home/home-tab.js';
import { classes as tree } from '../structure/_generated/scene-tree.module.css.js';
import { classes as tab } from '../structure/_generated/structure-tab.module.css.js';
import { classes } from './_generated/state-tree-tab.module.css.js';

import { derive } from '@molstar/prism-ui/reactive/state';
/** A shown cell: the node, its depth (`aria-level` − 1). Where it is among its siblings is in `Places`. */
interface CellRow {
  readonly node: StateTreeNode;
  readonly depth: number;
}

/** A shown cell's place among its siblings (`aria-posinset`, `aria-setsize`) and, if it opens, whether it is open. */
interface Place {
  readonly posinset: number;
  readonly setsize: number;
  readonly expanded?: boolean;
}

/** The shown cells' places by ref: a cell follows its own (view-models.md §3), so a reflow wakes only those that moved. */
type Places = KeyedState<string, Place>;

const samePlace = (a: Place, b: Place) =>
  a.posinset === b.posinset && a.setsize === b.setsize && a.expanded === b.expanded;

/**
 * The state tree, in the Structure tab's shape: the cells as rows (expand, show or hide, and a ⋯ with the updates, the
 * applicable actions and Remove), and under them at a fixed height the selected cell: its updates (its own transform's,
 * then the decorators folded into it) and its actions, each a form that folds open (Home's action blocks). The ⋯'s
 * entries select the cell and open theirs.
 *
 * Virtualized (view-models.md §3): the shown cells flattened into rows of one height (`VirtualList`), only the window
 * rendered, each row following its cell's record (a cell that stays the same object keeps its row); the keyboard
 * moves over the flattened rows (`TreeNavigation`), scrolling the target into view.
 */
export function StateTreeTab(props: { model: StateTreeModel }) {
  const { model } = props;
  // Collapsed cells (all open at first: the tree is short for a typical session): the view's own state, by ref
  // (view-models.md §3), changed in place.
  const collapsed = new Set<string>();
  // Which block of the inspector is open: the update (default) or an action.
  const focus = observableState<string>('update');
  const open = (ref: string, block: string) => {
    focus.set(block);
    model.select(ref);
  };
  const shown = derive([model.selected, focus.asReadonly], (ref, f) => (ref ? `${f}\n${ref}` : undefined));
  // The shown cells as rows: a cell that is the same object at the same depth is the same row (its place is apart, so
  // a sibling added does not render its siblings again).
  const keep = keepRows<CellRow>(
    (r) => r.node.ref,
    (a, b) => a.node === b.node && a.depth === b.depth,
  );
  const rows = observableState<readonly CellRow[]>([]);
  const places = keyedState<string, Place>(samePlace);
  // The selected cell, told to the two rows a selection concerns.
  const selected = keyedSelector(model.selected);
  const placeOf = (r: CellRow & { readonly posinset: number; readonly setsize: number }): Place => ({
    posinset: r.posinset,
    setsize: r.setsize,
    ...(r.node.children?.length ? { expanded: !collapsed.has(r.node.ref) } : {}),
  });
  let list: VirtualListHandle | undefined;
  const nav = new TreeNavigation<CellRow>({
    rows: () => rows.value,
    key: (r) => r.node.ref,
    level: (r) => r.depth + 1,
    scrollTo: (i) => list?.scrollToIndex(i),
    selected: () => model.selected.value,
  });
  // Cells opened or closed: each one's rows below put in place of the old (view-models.md §3, incremental), the rest
  // kept. Toggles in one task are one pass (a task's queue, rendering.md §2.7): a thousand cells closed at once copy the
  // rows once, not once a cell.
  const toggled = new Set<string>();
  const work = workQueue({ lane: 'task' });
  const applyToggles = () => {
    const before = rows.value;
    const next: (CellRow & { readonly posinset: number; readonly setsize: number })[] = [];
    const gone: string[] = [];
    const made: (CellRow & { readonly posinset: number; readonly setsize: number })[] = [];
    let at = 0;
    for (let i = 0; i < before.length; i++) {
      const row = before[i];
      if (!toggled.has(row.node.ref)) continue;
      let end = i + 1;
      while (end < before.length && before[end].depth > row.depth) end++;
      for (let j = at; j <= i; j++) next.push(before[j] as (typeof next)[number]);
      for (let j = i + 1; j < end; j++) gone.push(before[j].node.ref);
      const kids = row.node.children;
      if (kids?.length && !collapsed.has(row.node.ref))
        for (const r of flattenTree(
          kids,
          (n) => n.children,
          (n) => !collapsed.has(n.ref),
        )) {
          const shifted = { ...r, depth: r.depth + row.depth + 1 };
          next.push(shifted);
          made.push(shifted);
        }
      at = end;
      i = end - 1;
    }
    for (let j = at; j < before.length; j++) next.push(before[j] as (typeof next)[number]);
    const refs = [...toggled];
    toggled.clear();
    batch(() => {
      for (const ref of gone) places.delete(ref);
      for (const ref of refs) {
        const was = places.get(ref);
        if (was) places.set(ref, { ...was, expanded: !collapsed.has(ref) });
      }
      for (const r of made) places.set(r.node.ref, placeOf(r));
      rows.set(next);
    });
  };
  const toggle = (ref: string) => {
    if (collapsed.has(ref)) collapsed.delete(ref);
    else collapsed.add(ref);
    // Toggled twice in a task: as it was.
    if (toggled.has(ref)) toggled.delete(ref);
    else toggled.add(ref);
    work.want('toggle', applyToggles);
  };
  // The cells changed (the source's flush, coalesced): every row, the rows the same kept.
  bind(model.roots, (roots) => {
    const all = flattenTree(
      roots,
      (n) => n.children,
      (n) => !collapsed.has(n.ref),
    );
    const p = new Map<string, Place>();
    for (const r of all) p.set(r.node.ref, placeOf(r));
    batch(() => {
      places.replace(p);
      rows.set(keep(all));
    });
  });
  return (
    <div class={tab.tab}>
      <div class={tab.tree}>
        <Show when={derive(rows.asReadonly, (r) => r.length === 0)}>
          {() => <p class={tree.empty}>Nothing in the state yet.</p>}
        </Show>
        <VirtualList
          items={rows.asReadonly}
          key={(r) => r.node.ref}
          rowHeight={26}
          class={`${tree.tree} ${classes.tree}`}
          edge
          role="tree"
          aria-label="State tree"
          contentRef={nav.tree}
          ref={(l) => {
            list = l;
          }}
        >
          {(row) => (
            <Cell model={model} row={row} places={places} selected={selected} nav={nav} toggle={toggle} open={open} />
          )}
        </VirtualList>
      </div>
      <ScrollArea edge class={tab.inspector} contentClass={tab.inspectorContent} role="region" aria-label="Cell">
        <Show
          when={shown}
          fallback={() => <p class={tab.inspectorEmpty}>Select a cell to update it or apply an action to it.</p>}
        >
          {(key) => {
            const [block, ref] = key.split('\n');
            const inspection = model.inspect(ref);
            return inspection ? <Inspection inspection={inspection} focus={block} /> : null;
          }}
        </Show>
      </ScrollArea>
    </div>
  );
}

/**
 * A cell's row as view = f(state) (rendering.md §2): its record is fixed for the row (a changed cell is a new record,
 * so a new row), and one `watch` render writes what the tree's states change: its place, open state, selection and
 * its menu's.
 */
function Cell(props: {
  model: StateTreeModel;
  row: CellRow;
  places: Places;
  selected: KeyedSelector<string | undefined>;
  nav: TreeNavigation<CellRow>;
  toggle: (ref: string) => void;
  open: (ref: string, block: string) => void;
}) {
  const { model, row } = props;
  const { node, depth } = row;
  const hasChildren = !!node.children?.length;
  const menuOpen = observableState(false);
  let menuButton!: HTMLElement;
  let line!: HTMLElement;
  const hidden = node.visible === false;
  // Named by its own label and description.
  const labelId = uniqueId('msp-tree-label');
  const el = (
    <div
      role="treeitem"
      aria-level={depth + 1}
      aria-labelledby={`${labelId} ${labelId}-d`}
      tabIndex={-1}
      class={tree.item}
      ref={props.nav.row(node.ref)}
    >
      <div
        data-tip-row=""
        data-tree-row=""
        class={hidden ? `${tree.row} ${tree.hidden}` : tree.row}
        style={{ '--msp-tree-depth': depth }}
        ref={(e: HTMLElement) => (line = e)}
        onClick={() => model.select(node.ref)}
      >
        <button
          type="button"
          class={hasChildren ? tree.expander : `${tree.expander} ${tree.leaf}`}
          data-tree-toggle=""
          tabIndex={-1}
          aria-label="Expand"
          onClick={(e) => {
            e.stopPropagation();
            props.toggle(node.ref);
          }}
        >
          <Icon name="chevron-right" size={14} />
        </button>
        {depth === 0 ? <Icon name={node.icon} size={15} class={tree.kind} /> : null}
        <span class={tree.label} id={labelId} data-tip-part="">
          {node.label}
        </span>
        <span class={tree.description} id={`${labelId}-d`} data-tip-part="">
          {node.status === 'error' ? <span class={classes.error}>{node.error ?? 'Error'}</span> : node.description}
        </span>
        <span class={tree.actions}>
          {node.visible !== undefined ? (
            <button
              type="button"
              class={hidden ? `${tree.action} ${tree.keep}` : tree.action}
              title={hidden ? 'Show' : 'Hide'}
              onClick={(e) => {
                e.stopPropagation();
                model.toggleVisible(node.ref);
              }}
            >
              <Icon name={hidden ? 'eye-off' : 'eye'} size={14} />
            </button>
          ) : null}
          <button
            type="button"
            class={tree.action}
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
        </span>
      </div>
      <Show when={menuOpen}>
        {() => (
          <CellMenu
            model={model}
            node={node}
            anchor={menuButton}
            onClose={(refocus) => {
              menuOpen.set(false);
              if (refocus) menuButton.focus();
            }}
            onPick={(block) => {
              menuOpen.set(false);
              props.open(node.ref, block);
            }}
            onRemove={() => {
              menuOpen.set(false);
              model.remove(node.ref);
            }}
          />
        )}
      </Show>
    </div>
  ) as HTMLElement;
  // Its own keys of the tree's states: a selection tells two cells, a reflow the cells that moved.
  const place = props.places.key(node.ref);
  const isSelected = props.selected.key(node.ref);
  watch([place, isSelected, menuOpen], () => {
    const p = place.value;
    setAttr(el, 'aria-posinset', p?.posinset);
    setAttr(el, 'aria-setsize', p?.setsize);
    // A leaf has no `aria-expanded` (the pattern's end node).
    if (hasChildren) setAttr(el, 'aria-expanded', String(p?.expanded ?? true));
    const selected = isSelected.value;
    setAttr(el, 'aria-selected', String(selected));
    setClass(line, tree.picked, selected);
    const open = menuOpen.value;
    setClass(menuButton, tree.keep, open);
    setAttr(menuButton, 'aria-expanded', String(open));
  });
  return el;
}

/**
 * The updates, the applicable actions, then Remove (set apart, in the danger color): Mol*'s cell menu, as a popover
 * from the row's ⋯ (not inline under the row: rows are of one height in a virtualized tree).
 */
function CellMenu(props: {
  model: StateTreeModel;
  node: StateTreeNode;
  anchor: HTMLElement;
  onClose: (refocus?: boolean) => void;
  onPick: (block: string) => void;
  onRemove: () => void;
}) {
  const inspection = props.model.inspect(props.node.ref);
  const item = (id: string, label: string) => (
    <button type="button" role="menuitem" class={tree.menuItem} onClick={() => props.onPick(id)}>
      <span />
      <span>{label}</span>
    </button>
  );
  return (
    <MenuPopover anchor={props.anchor} onClose={props.onClose} aria-label={`${props.node.label}: actions`}>
      <div class={tree.menu} data-popover="">
        {inspection?.updates.map((u) => item(u.id, `${u.title}…`))}
        {inspection?.actions.map((a) => item(a.id, `${a.title}…`))}
        <button type="button" role="menuitem" class={tree.menuItem} data-danger="" onClick={props.onRemove}>
          <span />
          <span>Remove</span>
        </button>
      </div>
    </MenuPopover>
  );
}

/** The selected cell: what it is, its update form, its actions. */
function Inspection(props: { inspection: StateTreeInspection; focus: string }) {
  const { node, updates, actions } = props.inspection;
  return (
    <div class={tab.inspection}>
      <header class={tab.inspectorHeader}>
        <span class={tab.inspectorTitle}>{node.label}</span>
        <span class={classes.transform}>{node.transform}</span>
      </header>
      <div class={classes.blocks}>
        {updates.map((u, i) => (
          <ActionBlock action={u} open={props.focus === u.id || (i === 0 && props.focus === 'update')} />
        ))}
        {actions.length ? <p class={classes.caption}>Actions</p> : null}
        {actions.map((a) => (
          <ActionBlock action={a} open={props.focus === a.id} />
        ))}
      </div>
    </div>
  );
}
