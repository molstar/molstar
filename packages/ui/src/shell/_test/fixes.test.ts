/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { State } from '@molstar/core/reactive';
import { describe, expect, it } from 'vitest';
import {
  findView,
  type LayoutNode,
  mainNode,
  MinWindow,
  moveView,
  normalize,
  RightGroup,
  setWindowRect,
  type SplitNode,
} from '../layout.js';
import { ShellModel, type ShellViewDef } from '../shell-model.js';
import { paneMenu, viewMenu } from '../view-menu.js';
import { rightLayout } from './layouts.js';

// Regressions of the shell's review: what a pane's tab does, what the active and the maximized pane follow, which
// panes count, what a layout write keeps, what a registration replaces, what a reset or a restore forgets.

const def = (id: string, extra: Partial<ShellViewDef> = {}): ShellViewDef => ({
  id,
  title: id,
  icon: 'file',
  render: { mount: () => () => {} },
  ...extra,
});

/** Two stage views in one region group (tabs), a stage-homed `c`, and a right-column view. */
const tabbed = () => {
  const s = new ShellModel();
  s.registerViews([def('a', { home: 'stage' }), def('b', { home: 'stage' }), def('c', { home: 'stage' }), def('home')]);
  s.show('a');
  s.move('b', { kind: 'group', id: findView(s.layout.value, 'a')!.group.id });
  s.select('a');
  return s;
};

describe('a pane’s tabs', () => {
  it('pressing a region pane’s open tab never folds it; a side group’s does', () => {
    const s = tabbed();
    const group = () => findView(s.layout.value, 'a')!.group;
    s.select('a');
    expect(group().collapsed).toBeUndefined();
    expect(group().active).toBe('a');
    s.show('home');
    s.select('home');
    expect(findView(s.layout.value, 'home')!.group.collapsed).toBe(true);
  });

  it('the active pane follows a tab switch in its pane (and in a window)', () => {
    const s = tabbed();
    s.setActivePane('a');
    s.select('b');
    expect(s.activePane.value).toBe('b');
    expect(s.paneStatus('b').value.active).toBe(true);
    expect(s.paneStatus('a').value.active).toBe(false);
    s.float('b');
    s.move('c', { kind: 'group', id: findView(s.layout.value, 'b')!.group.id });
    s.setActivePane('c');
    s.select('b');
    expect(s.activePane.value).toBe('b');
  });

  it('a maximize is held by a view’s ref and follows its group’s tabs; maximizing a hidden tab shows it', () => {
    const s = tabbed();
    s.split('c', 'right');
    s.maximize('b');
    expect(findView(s.layout.value, 'b')!.group.active).toBe('b');
    expect(s.maximized.value).toBe('b');
    expect(s.paneStatus('b').value.maximized).toBe(true);
    s.select('a');
    expect(s.maximized.value).toBe('a');
    expect(paneMenu(s, 'a').find((i) => i.id === 'restore')).toBeDefined();
  });
});

describe('the panes that count', () => {
  it('a stage view closed from a window leaves the layout (not hidden in the region); one pane is not several', () => {
    const s = new ShellModel();
    s.registerViews([def('a', { home: 'stage' }), def('t', { home: 'top' })]);
    s.show('a');
    s.float('a');
    s.close('a');
    expect(findView(s.layout.value, 'a')).toBeUndefined();
    expect(s.layout.value.hidden).toEqual([]);
    s.show('t');
    s.float('t');
    s.close('t');
    expect(findView(s.layout.value, 't')).toBeUndefined();
    expect(mainNode(s.layout.value).root).toBeUndefined();
    expect(s.paneStatus('content').value.several).toBe(false);
  });

  it('a region pane that does not show (its view flagged off) is no pane: several, the corners and the grid skip it', () => {
    const visible = new State(true);
    const s = new ShellModel();
    s.registerViews([def('a', { home: 'stage', visible }), def('b', { home: 'stage', region: { side: 'bottom' } })]);
    s.show('a');
    s.show('b');
    s.setOverChrome({ bottomLeft: 100, bottomRight: 50 });
    expect(s.paneStatus('content').value.several).toBe(true);
    visible.set(false);
    expect(s.paneStatus('content').value.several).toBe(true);
    s.close('b');
    expect(s.paneStatus('content').value.several).toBe(false);
    // The content alone in what shows: both corners are its.
    expect(s.paneInsets('content').value).toEqual({ bottomLeft: 100, bottomRight: 50 });
    s.show('b');
    s.arrangeGrid();
    const root = mainNode(s.layout.value).root as SplitNode;
    // The grid's lines: the content and `b`; `a`, not showing, after them.
    expect(root.children.at(-1)).toMatchObject({ kind: 'group', views: ['a'] });
  });
});

describe('layout writes keep what did not change', () => {
  it('normalize keeps a normal layout and every node a write does not touch', () => {
    let l = rightLayout(['home', 'layers']);
    l = moveView(l, 'x', { kind: 'edge', of: 'content', side: 'right' });
    l = moveView(l, 'w', { kind: 'float', rect: { x: 10, y: 10, width: 300, height: 300 } });
    expect(normalize(l)).toBe(l);
    const moved = setWindowRect(l, l.windows[0].group.id, { x: 50, y: 50, width: 300, height: 300 });
    expect(moved.root).toBe(l.root);
    expect(moved.windows[0]).not.toBe(l.windows[0]);
  });

  it('a window’s move through the model rebuilds no tree node; a pane’s resize keeps the right column', () => {
    const s = new ShellModel();
    s.registerViews([def('home'), def('a', { home: 'stage' }), def('w', { home: 'float' })]);
    s.show('a');
    s.show('w');
    const l1 = s.layout.value;
    s.setWindowRect(l1.windows[0].group.id, { x: 50, y: 50, width: 300, height: 300 });
    expect(s.layout.value.root).toBe(l1.root);
    const right = (s.layout.value.root as SplitNode).children[1];
    s.resize(findView(s.layout.value, 'a')!.group.id, 200);
    expect((s.layout.value.root as SplitNode).children[1]).toBe(right);
    expect(mainNode(s.layout.value).root).not.toBe(mainNode(l1).root);
  });
});

describe('registrations', () => {
  it('a view registered again follows only the new `visible` flag; the old undo leaves the new one', () => {
    const s = new ShellModel();
    const old = new State(true);
    const undoOld = s.registerViews([def('a', { visible: old })]);
    const fresh = new State(true);
    const replacement = def('a', { visible: fresh });
    s.registerViews([replacement]);
    old.set(false);
    expect(s.absent.value.has('a')).toBe(false);
    fresh.set(false);
    expect(s.absent.value.has('a')).toBe(true);
    undoOld();
    expect(s.view('a')).toBe(replacement);
    expect(s.absent.value.has('a')).toBe(true);
  });
});

describe('restore and reset forget', () => {
  it('an instance’s subject and focus when it is not in the restored layout, and on reset', () => {
    const s = new ShellModel();
    s.registerViews([def('p', { instances: 'many' })]);
    const saved = s.save();
    s.show('p', { key: 'k', subject: 1, section: 'x' });
    s.restore(saved);
    expect(s.subject('p#k').value).toBeUndefined();
    expect(s.focus('p#k').value).toBeUndefined();
    s.show('p', { key: 'k', subject: 2 });
    s.resetLayout();
    expect(s.subject('p#k').value).toBeUndefined();
  });

  it('a grid asked for before the main area was measured, at a reset or a restore', () => {
    for (const end of ['reset', 'restore'] as const) {
      const s = new ShellModel();
      s.registerViews([def('a', { home: 'stage' }), def('b', { home: 'stage' })]);
      const saved = s.save();
      s.show('a');
      s.show('b');
      s.arrangeGrid();
      if (end === 'reset') s.resetLayout();
      else s.restore(saved);
      s.show('a');
      const before = mainNode(s.layout.value).root;
      s.setArea(800, 600);
      expect((mainNode(s.layout.value).root as LayoutNode & { grid?: true }).grid, end).toBeUndefined();
      expect(mainNode(s.layout.value).root?.kind).toBe(before?.kind);
    }
  });
});

describe('sizes and docking', () => {
  it('a view’s least size applies to every window it is in: floated, opened, resized', () => {
    const s = new ShellModel();
    const minSize = { width: 500, height: 400 };
    s.registerViews([def('a', { minSize }), def('f', { home: 'float', size: { width: 200, height: 100 }, minSize })]);
    s.float('a');
    expect(findView(s.layout.value, 'a')!.window!.rect).toMatchObject(minSize);
    s.show('f');
    const w = findView(s.layout.value, 'f')!.window!;
    expect(w.rect).toMatchObject(minSize);
    s.setWindowRect(w.group.id, { ...w.rect, width: MinWindow.width, height: MinWindow.height });
    expect(findView(s.layout.value, 'f')!.window!.rect).toMatchObject(minSize);
    s.move('a', { kind: 'float', rect: { x: 0, y: 0, width: 10, height: 10 } });
    expect(findView(s.layout.value, 'a')!.window!.rect).toMatchObject(minSize);
  });

  it('Dock right works for any docked view and any pane, whatever its home', () => {
    const s = new ShellModel();
    s.registerViews([def('home'), def('a', { home: 'stage' }), def('side', { home: 'left' })]);
    s.show('a');
    paneMenu(s, 'a')
      .find((i) => i.id === 'dock')!
      .run?.();
    expect(findView(s.layout.value, 'a')!.group.id).toBe(RightGroup);
    viewMenu(s, 'side')
      .find((i) => i.id === 'dock')!
      .run?.();
    expect(findView(s.layout.value, 'side')!.group.id).toBe(RightGroup);
  });
});
