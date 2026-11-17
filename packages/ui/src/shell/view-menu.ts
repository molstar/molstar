/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { ContextMenuItem } from '../components/context-menu/context-menu.js';
import { allGroups, columnOf, findView, inRegion, mainNode, type Rect, RightGroup, type ViewRef } from './layout.js';
import type { ShellModel } from './shell-model.js';

/** What a window's keyboard adjusting does (its menu's Move, Resize). */
export type WindowAdjust = 'move' | 'resize';

/**
 * A view's layout menu (shell §2: pop out and dock back are moves): Float or Dock; Split below, Move up, Move down in a
 * docked column; Move to another group; Hide; Reset layout. A window's adds Move, Resize (from the keyboard) and Reset
 * position.
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
  const items: ContextMenuItem[] = [];
  if (place.window) {
    items.push({ id: 'dock', label: 'Dock', icon: 'panel-right', run: () => shell.dock(ref) });
    if (options.adjust) {
      items.push({ id: 'move', label: 'Move with the keyboard', icon: 'move-3d', run: () => options.adjust?.('move') });
      items.push({
        id: 'resize',
        label: 'Resize with the keyboard',
        icon: 'maximize-2',
        run: () => options.adjust?.('resize'),
      });
    }
    if (options.resetPosition) {
      items.push({ id: 'reset-position', label: 'Reset position', icon: 'rotate-ccw', run: options.resetPosition });
    }
  } else {
    items.push({
      id: 'float',
      label: 'Float',
      icon: 'app-window',
      run: () => shell.float(ref, options.anchorRect?.()),
    });
    if (place.group.views.length > 1) {
      items.push({ id: 'split', label: 'Split below', icon: 'panel-bottom', run: () => shell.splitBelow(ref) });
    }
    const column = columnOf(layout, place.group.id);
    if (column) {
      const groups = column.children.filter((c) => c.kind === 'group');
      const i = groups.findIndex((g) => g.id === place.group.id);
      if (i > 0) items.push({ id: 'up', label: 'Move up', icon: 'arrow-up', run: () => shell.moveInColumn(ref, -1) });
      if (i >= 0 && i < groups.length - 1) {
        items.push({ id: 'down', label: 'Move down', icon: 'arrow-down', run: () => shell.moveInColumn(ref, 1) });
      }
    }
    if (place.group.id !== RightGroup && (def.home ?? RightGroup) === RightGroup) {
      items.push({ id: 'dock', label: 'Dock right', icon: 'panel-right', run: () => shell.dock(ref) });
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
      label: `Move to the ${first.title} window`,
      icon: 'app-window',
      run: () => shell.move(ref, { kind: 'group', id: g.id }),
    });
  }
  if (!place.window && (def.instances ?? 'single') === 'single') {
    items.push({ id: 'hide', label: 'Hide', icon: 'eye-off', run: () => shell.hide(ref) });
  }
  items.push({ id: 'reset', label: 'Reset layout', icon: 'list-restart', run: () => shell.resetLayout() });
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
  const items: ContextMenuItem[] = [
    { id: 'stage-right', label: 'Beside the main view', icon: 'panel-right', run: () => shell.split(ref, 'right') },
    { id: 'stage-below', label: 'Below the main view', icon: 'panel-bottom', run: () => shell.split(ref, 'bottom') },
  ];
  if (shell.top !== 'overlay') {
    items.push({
      id: 'stage-above',
      label: 'Above the main view',
      icon: 'panel-top',
      run: () => shell.split(ref, 'top'),
    });
  }
  if (shell.top !== 'split') {
    items.push({ id: 'stage-top', label: 'Over the main view’s top', icon: 'panel-top', run: () => shell.toTop(ref) });
  }
  return items;
}

/**
 * A pane's menu (shell §2.1: a view in the main area's region, or on its top): Maximize or Show every pane; split it
 * beside another pane's; Float; back to its home (Dock right); onto or off the top; Close.
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
  const items: ContextMenuItem[] = [];
  const onTop = mainNode(l).top?.id === place.group.id;
  if (!onTop) {
    const maximized = shell.maximized.value === place.group.id;
    items.push(
      maximized
        ? { id: 'restore', label: 'Show every pane', icon: 'minimize-2', run: () => shell.maximize(undefined) }
        : { id: 'maximize', label: 'Maximize', icon: 'maximize-2', run: () => shell.maximize(place.group.id) },
    );
    if (place.group.views.length > 1) {
      items.push(
        { id: 'split-right', label: 'Split right', icon: 'panel-right', run: () => shell.split(ref, 'right', ref) },
        { id: 'split-below', label: 'Split below', icon: 'panel-bottom', run: () => shell.split(ref, 'bottom', ref) },
      );
    }
    if (shell.top !== 'split') {
      items.push({
        id: 'stage-top',
        label: 'Over the main view’s top',
        icon: 'panel-top',
        run: () => shell.toTop(ref),
      });
    }
    items.push({ id: 'grid', label: 'Arrange as a grid', icon: 'layout-grid', run: () => shell.arrangeGrid() });
  } else if (shell.top !== 'overlay') {
    items.push({
      id: 'stage-above',
      label: 'Above the main view',
      icon: 'panel-top',
      run: () => shell.split(ref, 'top'),
    });
  }
  items.push({ id: 'float', label: 'Float', icon: 'app-window', run: () => shell.float(ref, options.anchorRect?.()) });
  if ((def.home ?? RightGroup) === RightGroup) {
    items.push({ id: 'dock', label: 'Dock right', icon: 'panel-right', run: () => shell.dock(ref) });
  }
  items.push({ id: 'close', label: 'Close', icon: 'x', run: () => shell.close(ref) });
  return items;
}

/**
 * A rail button's menu (a context menu, Shift+F10): Hide, Show each hidden view, Hide the rail (the Mol* star in the main
 * area's corner brings it back), Reset layout.
 */
export function railMenu(shell: ShellModel, ref: ViewRef | undefined): ContextMenuItem[] {
  const items: ContextMenuItem[] = [];
  if (ref) items.push({ id: 'hide', label: 'Hide', icon: 'eye-off', run: () => shell.hide(ref) });
  for (const h of shell.layout.value.hidden) {
    const def = shell.view(h);
    if (def) items.push({ id: `show:${h}`, label: `Show ${def.title}`, icon: def.icon, run: () => shell.show(h) });
  }
  items.push({ id: 'hide-rail', label: 'Hide the rail', icon: 'eye-off', run: () => shell.toggleRail(false) });
  items.push({ id: 'reset', label: 'Reset layout', icon: 'list-restart', run: () => shell.resetLayout() });
  return items;
}
