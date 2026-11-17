/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { ContextMenuItem } from '../components/context-menu/context-menu.js';
import { allGroups, columnOf, findView, inRegion, mainNode, type Rect, RightGroup, type ViewRef } from './layout.js';
import { keyOf, type ShellModel, viewIdOf } from './shell-model.js';

/** What a window's keyboard adjusting does (its menu's Move, Resize). */
export type WindowAdjust = 'move' | 'resize';

/**
 * A view's layout menu (shell §2, §3.3: pop out and dock back are moves): a docked view's Float, Split below, Move up,
 * Move down (in a column), Dock right (away from the right group); a window's Dock (to its home), Move and Resize (from
 * the keyboard), Reset position; then into the main area's region (§2.1), Move to each other window, Hide (a single
 * view docked: `close`), Reset layout.
 */
export function viewMenu(
  shell: ShellModel,
  ref: ViewRef,
  options: { anchorRect?: () => Rect | undefined; adjust?(how: WindowAdjust): void; resetPosition?(): void } = {},
): ContextMenuItem[] {
  const layout = shell.layout.value;
  const place = findView(layout, ref);
  const def = shell.view(ref);
  if (!place || !def) return [];
  const t = shell.labels;
  const items: ContextMenuItem[] = [];
  if (place.window) {
    items.push({ id: 'dock', label: t.dock, icon: 'panel-right', run: () => shell.dock(ref) });
    if (options.adjust) {
      items.push({ id: 'move', label: t.moveWithKeyboard, icon: 'move-3d', run: () => options.adjust?.('move') });
      items.push({
        id: 'resize',
        label: t.resizeWithKeyboard,
        icon: 'maximize-2',
        run: () => options.adjust?.('resize'),
      });
    }
    if (options.resetPosition) {
      items.push({ id: 'reset-position', label: t.resetPosition, icon: 'rotate-ccw', run: options.resetPosition });
    }
  } else {
    items.push({
      id: 'float',
      label: t.float,
      icon: 'app-window',
      run: () => shell.float(ref, options.anchorRect?.()),
    });
    if (place.group.views.length > 1) {
      items.push({ id: 'split', label: t.splitBelow, icon: 'panel-bottom', run: () => shell.splitBelow(ref) });
    }
    const column = columnOf(layout, place.group.id);
    if (column) {
      const groups = column.children.filter((c) => c.kind === 'group');
      const i = groups.findIndex((g) => g.id === place.group.id);
      if (i > 0) items.push({ id: 'up', label: t.moveUp, icon: 'arrow-up', run: () => shell.moveInColumn(ref, -1) });
      if (i >= 0 && i < groups.length - 1) {
        items.push({ id: 'down', label: t.moveDown, icon: 'arrow-down', run: () => shell.moveInColumn(ref, 1) });
      }
    }
    if (place.group.id !== RightGroup) {
      items.push({ id: 'dock', label: t.dockRight, icon: 'panel-right', run: () => shell.dock(ref, RightGroup) });
    }
  }
  items.push(...stageItems(shell, ref));
  // Into a window that is open (a window holds several views as tabs).
  for (const g of allGroups(layout)) {
    if (g.id === place.group.id || !layout.windows.some((w) => w.group.id === g.id)) continue;
    const first = shell.view(g.active ?? g.views[0]);
    if (!first) continue;
    items.push({
      id: `into:${g.id}`,
      label: t.moveToWindow(first.title),
      icon: 'app-window',
      run: () => shell.move(ref, { kind: 'group', id: g.id }),
    });
  }
  if (!place.window && (def.instances ?? 'single') === 'single') {
    items.push({ id: 'hide', label: t.hide, icon: 'eye-off', run: () => shell.close(ref) });
  }
  items.push({ id: 'reset', label: t.resetLayout, icon: 'list-restart', run: () => shell.resetLayout() });
  return items;
}

/**
 * Into the main area's region (shell §2.1): beside the app's content, above it, on its top (as the shell's `top`
 * option allows). Not for a view already there (its pane menu moves it).
 */
function stageItems(shell: ShellModel, ref: ViewRef): ContextMenuItem[] {
  const l = shell.layout.value;
  const place = findView(l, ref);
  if (place && !place.window && inRegion(l, place.group.id)) return [];
  const t = shell.labels;
  const items: ContextMenuItem[] = [
    { id: 'stage-right', label: t.besideMain, icon: 'panel-right', run: () => shell.split(ref, 'right') },
    { id: 'stage-below', label: t.belowMain, icon: 'panel-bottom', run: () => shell.split(ref, 'bottom') },
  ];
  if (shell.top !== 'overlay') {
    items.push({
      id: 'stage-above',
      label: t.aboveMain,
      icon: 'panel-top',
      run: () => shell.split(ref, 'top'),
    });
  }
  if (shell.top !== 'split') {
    items.push({ id: 'stage-top', label: t.overMainTop, icon: 'panel-top', run: () => shell.toTop(ref) });
  }
  return items;
}

/**
 * A pane's menu (shell §2.1: a view in the main area's region, or on its top): Maximize or Show every pane, Split right
 * and below (a group of several), onto the top, Arrange as a grid (or, on the top, Above the main view); Float; Dock
 * right; Close.
 */
export function paneMenu(
  shell: ShellModel,
  ref: ViewRef,
  options: { anchorRect?: () => Rect | undefined } = {},
): ContextMenuItem[] {
  const l = shell.layout.value;
  const place = findView(l, ref);
  const def = shell.view(ref);
  if (!place || !def) return [];
  const t = shell.labels;
  const items: ContextMenuItem[] = [];
  const onTop = mainNode(l).top?.id === place.group.id;
  if (!onTop) {
    items.push(
      shell.maximized.value === ref
        ? { id: 'restore', label: t.showEveryPane, icon: 'minimize-2', run: () => shell.maximize(undefined) }
        : { id: 'maximize', label: t.maximize, icon: 'maximize-2', run: () => shell.maximize(ref) },
    );
    if (place.group.views.length > 1) {
      items.push(
        { id: 'split-right', label: t.splitRight, icon: 'panel-right', run: () => shell.split(ref, 'right', ref) },
        { id: 'split-below', label: t.splitBelow, icon: 'panel-bottom', run: () => shell.split(ref, 'bottom', ref) },
      );
    }
    if (shell.top !== 'split') {
      items.push({
        id: 'stage-top',
        label: t.overMainTop,
        icon: 'panel-top',
        run: () => shell.toTop(ref),
      });
    }
    items.push({ id: 'grid', label: t.arrangeGrid, icon: 'layout-grid', run: () => shell.arrangeGrid() });
  } else if (shell.top !== 'overlay') {
    items.push({
      id: 'stage-above',
      label: t.aboveMain,
      icon: 'panel-top',
      run: () => shell.split(ref, 'top'),
    });
  }
  items.push({ id: 'float', label: t.float, icon: 'app-window', run: () => shell.float(ref, options.anchorRect?.()) });
  items.push({ id: 'dock', label: t.dockRight, icon: 'panel-right', run: () => shell.dock(ref, RightGroup) });
  items.push({ id: 'close', label: t.closePane, icon: 'x', run: () => shell.close(ref) });
  return items;
}

/**
 * A rail button's menu (a context menu, Shift+F10): Hide (`close`), Show each hidden view, Hide the rail (the rail
 * toggle in the main area's corner brings it back), Reset layout.
 */
export function railMenu(shell: ShellModel, ref: ViewRef | undefined): ContextMenuItem[] {
  const t = shell.labels;
  const items: ContextMenuItem[] = [];
  if (ref) items.push({ id: 'hide', label: t.hide, icon: 'eye-off', run: () => shell.close(ref) });
  for (const h of shell.layout.value.hidden) {
    const def = shell.view(h);
    // An instance shows as itself (its key kept).
    if (def) {
      items.push({
        id: `show:${h}`,
        label: t.showView(def.title),
        icon: def.icon,
        run: () => shell.show(viewIdOf(h), { key: keyOf(h) }),
      });
    }
  }
  items.push({ id: 'hide-rail', label: t.hideRail, icon: 'eye-off', run: () => shell.toggleRail(false) });
  items.push({ id: 'reset', label: t.resetLayout, icon: 'list-restart', run: () => shell.resetLayout() });
  return items;
}
