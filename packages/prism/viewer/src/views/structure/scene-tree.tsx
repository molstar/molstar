/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { combineLatest, distinctUntilChanged, map } from 'rxjs';
import { Icon } from '@molstar/prism-ui/components/icon/icon';
import { uniqueId } from '@molstar/prism-ui/dom/focus';
import { For, Show } from '@molstar/prism-ui/dom/control-flow';
import { onCleanup } from '@molstar/prism-ui/dom/scope';
import type { SceneMenuItem, SceneModel, SceneNode } from '../../model/scene.js';
import { observableState } from '@molstar/prism-ui/reactive/observable-state';
import { classes } from './_generated/scene-tree.module.css.js';
import { treeFocus, treeKeys } from '@molstar/prism-ui/dom/tree-keys';

/**
 * The scene as one tree: structures → components → representations, volumes, measurements. Rows are keyed by node, so
 * a hierarchy update touches only the rows that changed; hover highlights in 3D, a click picks the node for the
 * inspector, visibility toggles inline.
 */
export function SceneTree(props: { scene: SceneModel }) {
  const { scene } = props;
  return (
    <div class={classes.tree} role="tree" aria-label="Scene" tabIndex={0} onKeyDown={treeKeys} onFocusIn={treeFocus}>
      <For each={scene.roots} key={(n) => n.id} fallback={() => <p class={classes.empty}>Nothing loaded yet.</p>}>
        {(node) => <TreeNode scene={scene} node={node} depth={0} />}
      </For>
    </div>
  );
}

function TreeNode(props: { scene: SceneModel; node: SceneNode; depth: number }) {
  const { scene, node, depth } = props;
  const hasChildren = node.children.pipe(
    map((c) => c.length > 0 || !!node.lazy),
    distinctUntilChanged(),
  );
  const open = combineLatest([node.isExpanded, hasChildren]).pipe(map(([e, h]) => e && h));
  const menuOpen = observableState(false);
  const picked = scene.picked.pipe(
    map((p) => p === node),
    distinctUntilChanged(),
  );
  // Named by its own label and description (not the whole subtree it holds).
  const labelId = uniqueId('msp-tree-label');
  return (
    <div
      role="treeitem"
      aria-expanded={open}
      aria-selected={picked}
      aria-labelledby={`${labelId} ${labelId}-d`}
      tabIndex={-1}
      class={classes.item}
    >
      <div
        data-tip-row=""
        data-tree-row=""
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
          classList={{ [classes.leaf]: hasChildren.pipe(map((h) => !h)) }}
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
              aria-expanded={menuOpen}
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
        {() => <RowMenu items={menuOf(node)} depth={depth} onDone={() => menuOpen.set(false)} />}
      </Show>
      <Show when={open} reveal>
        {() => {
          // Children computed on first look (lazy nodes); until they arrive, a quiet loading row.
          node.load?.();
          return (
            <div role="group">
              <For
                each={node.children}
                key={(n) => n.id}
                fallback={() => (
                  <div class={classes.loading} style={{ '--msp-tree-depth': depth + 1 }}>
                    Loading…
                  </div>
                )}
              >
                {(child) => <TreeNode scene={scene} node={child} depth={depth + 1} />}
              </For>
            </div>
          );
        }}
      </Show>
    </div>
  );
}

/** A row's menu: its own entries, then Remove (set apart, in the danger color) when it can be removed. */
function menuOf(node: SceneNode): SceneMenuItem[] {
  const own = node.menu ?? [];
  if (node.removable === false) return [...own];
  return [...own, { id: 'remove', label: 'Remove', icon: 'trash-2', danger: true, run: () => node.remove() }];
}

/**
 * A row's menu, inline under the row (nothing to clip in a scrolling tree): actions run and close it; entries with
 * sub-entries open them in place, indented.
 */
function RowMenu(props: { items: readonly SceneMenuItem[]; depth: number; onDone(): void }) {
  // Escape closes the menu first, before the panel it is in (the document's handlers run in capture order).
  const onKey = (e: KeyboardEvent) => {
    if (e.key !== 'Escape' || e.defaultPrevented) return;
    e.preventDefault();
    props.onDone();
  };
  document.addEventListener('keydown', onKey, true);
  onCleanup(() => document.removeEventListener('keydown', onKey, true));
  return (
    <div class={classes.menu} role="menu" style={{ '--msp-tree-depth': props.depth }}>
      <MenuItems items={props.items} onDone={props.onDone} />
    </div>
  );
}

function MenuItems(props: { items: readonly SceneMenuItem[]; onDone(): void }) {
  return (
    <>
      {props.items.map((item) =>
        item.items?.length ? (
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
      role="menuitem"
      class={classes.menuItem}
      data-danger={item.danger ? '' : undefined}
      onClick={() => {
        item.run?.();
        props.onDone();
      }}
    >
      {/* The column a submenu's chevron takes; items carry no icons. */}
      <span />
      <span>{item.label}</span>
    </button>
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
