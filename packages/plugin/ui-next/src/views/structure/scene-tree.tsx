/** @jsxImportSource @molstar/plugin-ui-next/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { combineLatest, distinctUntilChanged, map } from 'rxjs';
import { Icon } from '../../components/icon/icon.js';
import { For, Show } from '../../dom/control-flow.js';
import type { SceneMenuItem, SceneModel, SceneNode } from '../../model/scene.js';
import { observableState } from '../../reactive/observable-state.js';
import { classes } from './_generated/scene-tree.module.css.js';

/**
 * The scene as one tree: structures → components → representations, volumes, measurements. Rows are keyed by node, so
 * a hierarchy update touches only the rows that changed; hover highlights in 3D, a click picks the node for the
 * inspector, visibility toggles inline.
 */
export function SceneTree(props: { scene: SceneModel }) {
  const { scene } = props;
  return (
    <div class={classes.tree} role="tree" aria-label="Scene">
      <For each={scene.roots} key={(n) => n.id} fallback={() => <p class={classes.empty}>Nothing loaded yet.</p>}>
        {(node) => <TreeNode scene={scene} node={node} depth={0} />}
      </For>
    </div>
  );
}

function TreeNode(props: { scene: SceneModel; node: SceneNode; depth: number }) {
  const { scene, node, depth } = props;
  const hasChildren = node.children.pipe(
    map((c) => c.length > 0),
    distinctUntilChanged(),
  );
  const open = combineLatest([node.isExpanded, hasChildren]).pipe(map(([e, h]) => e && h));
  const menuOpen = observableState(false);
  const picked = scene.picked.pipe(
    map((p) => p === node),
    distinctUntilChanged(),
  );
  return (
    <div role="treeitem" aria-expanded={open} aria-selected={picked}>
      <div
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
          aria-label="Expand"
          onClick={(e) => {
            e.stopPropagation();
            node.toggleExpanded();
          }}
        >
          <Icon name="chevron-right" size={14} />
        </button>
        <Icon name={node.icon} size={15} class={classes.kind} />
        <span class={classes.label}>{node.label}</span>
        <span class={classes.description}>{node.description}</span>
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
          {node.menu?.length ? (
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
          ) : null}
          <button
            type="button"
            class={classes.action}
            title="Remove"
            onClick={(e) => {
              e.stopPropagation();
              node.remove();
            }}
          >
            <Icon name="trash-2" size={14} />
          </button>
        </span>
      </div>
      {node.menu?.length ? (
        <Show when={menuOpen}>
          {() => <RowMenu items={node.menu!} depth={depth} onDone={() => menuOpen.set(false)} />}
        </Show>
      ) : null}
      <Show when={open}>
        {() => (
          <div role="group">
            <For each={node.children} key={(n) => n.id}>
              {(child) => <TreeNode scene={scene} node={child} depth={depth + 1} />}
            </For>
          </div>
        )}
      </Show>
    </div>
  );
}

/**
 * A row's menu, inline under the row (nothing to clip in a scrolling tree): actions run and close it; entries with
 * sub-entries open them in place, indented.
 */
function RowMenu(props: { items: readonly SceneMenuItem[]; depth: number; onDone(): void }) {
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
      onClick={() => {
        item.run?.();
        props.onDone();
      }}
    >
      {item.icon ? <Icon name={item.icon} size={14} /> : <span />}
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
      <Show when={open}>
        {() => (
          <div class={classes.subItems} role="menu" aria-label={item.label}>
            <MenuItems items={item.items!} onDone={props.onDone} />
          </div>
        )}
      </Show>
    </div>
  );
}
