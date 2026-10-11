/** @jsxImportSource @molstar/plugin-ui-next/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { combineLatest, distinctUntilChanged, map } from 'rxjs';
import { Icon, type IconName } from '../../components/icon/icon.js';
import { For, Show } from '../../dom/control-flow.js';
import type { StateTreeInspection, StateTreeModel, StateTreeNode } from '../../model/state-tree.js';
import { type ObservableState, observableState } from '../../reactive/observable-state.js';
import { ActionBlock } from '../home/home-tab.js';
import { classes as tree } from '../structure/_generated/scene-tree.module.css.js';
import { classes as tab } from '../structure/_generated/structure-tab.module.css.js';
import { classes } from './_generated/state-tree-tab.module.css.js';

/**
 * The state tree, in the Structure tab's shape: the cells as rows (expand, show or hide, and a ⋯ with Update, the
 * applicable actions and Remove), and under them at a fixed height the selected cell: its transform's params to update
 * and its actions, each a form that folds open (Home's action blocks). The ⋯'s entries select the cell and open theirs.
 */
export function StateTreeTab(props: { model: StateTreeModel }) {
  const { model } = props;
  // Collapsed cells (all open at first: the tree is short for a typical session).
  const collapsed = observableState<ReadonlySet<string>>(new Set());
  // Which block of the inspector is open: the update (default) or an action.
  const focus = observableState<string>('update');
  const open = (ref: string, block: string) => {
    focus.set(block);
    model.select(ref);
  };
  const shown = combineLatest([model.selected, focus.asReadonly]).pipe(
    map(([ref, f]) => (ref ? `${f}\n${ref}` : undefined)),
    distinctUntilChanged(),
  );
  return (
    <div class={tab.tab}>
      <div class={tab.tree}>
        <div class={tree.tree} role="tree" aria-label="State tree">
          <For
            each={model.roots}
            key={(n) => n.ref}
            fallback={() => <p class={tree.empty}>Nothing in the state yet.</p>}
          >
            {(node) => <Cell model={model} node={node} depth={0} collapsed={collapsed} open={open} />}
          </For>
        </div>
      </div>
      <section class={tab.inspector} aria-label="Cell">
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
      </section>
    </div>
  );
}

function Cell(props: {
  model: StateTreeModel;
  node: StateTreeNode;
  depth: number;
  collapsed: ObservableState<ReadonlySet<string>>;
  open: (ref: string, block: string) => void;
}) {
  const { model, node, depth, collapsed } = props;
  const hasChildren = !!node.children?.length;
  const expanded = collapsed.asReadonly.pipe(
    map((c) => hasChildren && !c.has(node.ref)),
    distinctUntilChanged(),
  );
  const selected = model.selected.pipe(
    map((s) => s === node.ref),
    distinctUntilChanged(),
  );
  const menuOpen = observableState(false);
  const hidden = node.visible === false;
  const toggle = () => {
    const next = new Set(collapsed.value);
    if (next.has(node.ref)) next.delete(node.ref);
    else next.add(node.ref);
    collapsed.set(next);
  };
  return (
    <div role="treeitem" aria-expanded={expanded} aria-selected={selected}>
      <div
        class={tree.row}
        classList={{ [tree.hidden]: hidden, [tree.picked]: selected }}
        style={{ '--msp-tree-depth': depth }}
        onClick={() => model.select(node.ref)}
      >
        <button
          type="button"
          class={tree.expander}
          classList={{ [tree.leaf]: !hasChildren }}
          aria-label="Expand"
          onClick={(e) => {
            e.stopPropagation();
            toggle();
          }}
        >
          <Icon name="chevron-right" size={14} />
        </button>
        <Icon name={node.icon} size={15} class={tree.kind} />
        <span class={tree.label}>{node.label}</span>
        <span class={tree.description}>
          {node.status === 'error' ? <span class={classes.error}>{node.error ?? 'Error'}</span> : node.description}
        </span>
        <span class={tree.actions}>
          {node.visible !== undefined ? (
            <button
              type="button"
              class={tree.action}
              classList={{ [tree.keep]: hidden }}
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
            classList={{ [tree.keep]: menuOpen }}
            title="More"
            aria-expanded={menuOpen}
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
            depth={depth}
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
      <Show when={expanded}>
        {() => (
          <div role="group">
            <For each={node.children ?? []} key={(n) => n.ref}>
              {(child) => <Cell {...props} node={child} depth={depth + 1} />}
            </For>
          </div>
        )}
      </Show>
    </div>
  );
}

/** Update, the applicable actions, then Remove (set apart, in the danger color): Mol*'s cell menu, inline. */
function CellMenu(props: {
  model: StateTreeModel;
  node: StateTreeNode;
  depth: number;
  onPick: (block: string) => void;
  onRemove: () => void;
}) {
  const inspection = props.model.inspect(props.node.ref);
  const item = (id: string, label: string, icon: IconName) => (
    <button type="button" role="menuitem" class={tree.menuItem} onClick={() => props.onPick(id)}>
      <Icon name={icon} size={14} />
      <span>{label}</span>
    </button>
  );
  return (
    <div class={tree.menu} role="menu" style={{ '--msp-tree-depth': props.depth }}>
      {inspection?.update ? item('update', `${inspection.update.title}…`, 'sliders-horizontal') : null}
      {inspection?.actions.map((a) => item(a.id, `${a.title}…`, a.icon))}
      <button type="button" role="menuitem" class={tree.menuItem} data-danger="" onClick={props.onRemove}>
        <Icon name="trash-2" size={14} />
        <span>Remove</span>
      </button>
    </div>
  );
}

/** The selected cell: what it is, its update form, its actions. */
function Inspection(props: { inspection: StateTreeInspection; focus: string }) {
  const { node, update, actions } = props.inspection;
  return (
    <div class={tab.inspection}>
      <header class={tab.inspectorHeader}>
        <Icon name={node.icon} size={15} />
        <span class={tab.inspectorTitle}>{node.label}</span>
        <span class={classes.transform}>{node.transform}</span>
      </header>
      <div class={classes.blocks}>
        {update ? <ActionBlock action={update} open={props.focus === 'update'} /> : null}
        {actions.length ? <p class={classes.caption}>Actions</p> : null}
        {actions.map((a) => (
          <ActionBlock action={a} open={props.focus === a.id} />
        ))}
      </div>
    </div>
  );
}
