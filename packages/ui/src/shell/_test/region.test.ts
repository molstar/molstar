/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import {
  arrangeGrid,
  TopGroup,
  defaultLayout,
  findView,
  inRegion,
  type LayoutNode,
  mainNode,
  moveView,
  normalize,
  removeView,
  resizeNode,
  type ShellLayout,
  validateLayout,
} from '../layout.js';
import { PaneMap } from '../pane-map.js';
import { shellCommands } from '../shell-commands.js';
import { ShellModel, type ShellViewDef } from '../shell-model.js';

// The main area's region (shell §2.1): views beside the app's content, the content always there and taking the rest,
// no tree once only the content is left; the top overlay; stored layouts (an old bare main stays as it was).

const base = () => defaultLayout(['home', 'layers']);
const leaves = (n: LayoutNode | undefined): string[] =>
  !n
    ? []
    : n.kind === 'split'
      ? n.children.flatMap(leaves)
      : n.kind === 'group'
        ? [`g:${n.views.join('+')}`]
        : [n.kind];

describe('the main area’s region', () => {
  it('a view docked beside the content starts a tree: the content takes the rest, the view has a size', () => {
    const l = moveView(base(), 'seq', { kind: 'edge', of: 'content', side: 'top', size: 96 });
    const region = mainNode(l).root!;
    expect(region.kind).toBe('split');
    expect(leaves(region)).toEqual(['g:seq', 'content']);
    if (region.kind !== 'split') return;
    expect(region.dir).toBe('column');
    // A pane's least size is a strip's, not a panel's (`MinPaneSize`).
    expect((region.children[0] as { size?: number }).size).toBe(96);
    expect((region.children[1] as { size?: number }).size).toBeUndefined();
    expect(findView(l, 'seq')).toBeDefined();
    expect(inRegion(l, findView(l, 'seq')!.group.id)).toBe(true);
    expect(inRegion(l, 'right')).toBe(false);
  });

  it('splits again beside a view in the region (a stage beside a stage, a plot beside it)', () => {
    let l = moveView(base(), '3d#b', { kind: 'edge', of: 'content', side: 'right' });
    const b = findView(l, '3d#b')!.group.id;
    l = moveView(l, 'plot', { kind: 'edge', of: b, side: 'bottom' });
    expect(leaves(mainNode(l).root)).toEqual(['content', 'g:3d#b', 'g:plot']);
  });

  it('removing the last docked view leaves no tree: the content alone, as before', () => {
    const l = moveView(base(), 'seq', { kind: 'edge', of: 'content', side: 'top' });
    const back = removeView(l, 'seq');
    expect(mainNode(back).root).toBeUndefined();
    expect(back.root).toEqual(base().root);
  });

  it('keeps exactly one content in the region and none outside it; one main at the top', () => {
    const l = moveView(base(), 'seq', { kind: 'edge', of: 'content', side: 'top' });
    const main = mainNode(l);
    // A second content and a stray main in the region go; a content outside the region goes.
    const broken: ShellLayout = {
      ...l,
      root: {
        kind: 'split',
        id: 'root',
        dir: 'row',
        children: [
          {
            ...main,
            root: {
              kind: 'split',
              id: 'r',
              dir: 'row',
              children: [main.root!, { kind: 'content', id: 'content' }, { kind: 'main', id: 'main' }],
            },
          },
          { kind: 'content', id: 'content' },
          { kind: 'group', id: 'right', views: ['home'] },
        ],
      },
    };
    const n = normalize(broken);
    expect(leaves(mainNode(n).root)).toEqual(['g:seq', 'content']);
    expect(leaves(n.root).filter((x) => x === 'main' || x === 'content')).toEqual(['main']);
  });

  it('a region without its content gets it back at the left', () => {
    const l = normalize({
      ...base(),
      root: {
        kind: 'split',
        id: 'root',
        dir: 'row',
        children: [
          { kind: 'main', id: 'main', root: { kind: 'group', id: 'p', views: ['plot'] } },
          { kind: 'group', id: 'right', views: ['home'], size: 340 },
        ],
      },
    });
    expect(leaves(mainNode(l).root)).toEqual(['content', 'g:plot']);
  });

  it('a region group resizes; the content has no size to set', () => {
    let l = moveView(base(), 'seq', { kind: 'edge', of: 'content', side: 'top', size: 200 });
    const g = findView(l, 'seq')!.group.id;
    l = resizeNode(l, g, 260);
    expect(findView(l, 'seq')!.group.size).toBe(260);
    expect(resizeNode(l, 'content', 400)).toEqual(l);
    expect(findView(resizeNode(l, g, 10), 'seq')!.group.size).toBe(64);
  });

  it('the top overlay: a group of its own, found, emptied away', () => {
    let l = moveView(base(), 'seq', { kind: 'top' });
    expect(mainNode(l).top?.views).toEqual(['seq']);
    expect(findView(l, 'seq')?.group.id).toBe(TopGroup);
    expect(inRegion(l, TopGroup)).toBe(true);
    l = moveView(l, 'seq', { kind: 'edge', of: 'content', side: 'top' });
    expect(mainNode(l).top).toBeUndefined();
    expect(inRegion(l, findView(l, 'seq')!.group.id)).toBe(true);
  });

  it('a stored layout keeps its region and top overlay; an old one (a bare main) reads as before', () => {
    let l = moveView(base(), '3d#b', { kind: 'edge', of: 'content', side: 'right' });
    l = moveView(l, 'seq', { kind: 'top' });
    const back = validateLayout(JSON.parse(JSON.stringify(l)))!;
    expect(back).toEqual(l);
    const old = validateLayout(JSON.parse(JSON.stringify(base())))!;
    expect(mainNode(old).root).toBeUndefined();
    expect(old).toEqual(base());
  });
});

describe('the region as a grid', () => {
  /** A region of the content and three stage views. */
  const four = () => {
    let l = moveView(base(), '3d#b', { kind: 'edge', of: 'content', side: 'right' });
    l = moveView(l, '3d#c', { kind: 'edge', of: 'content', side: 'bottom' });
    return moveView(l, '3d#d', { kind: 'edge', of: findView(l, '3d#b')!.group.id, side: 'bottom' });
  };
  const area = { width: 1000, height: 800 };
  const grid = (l: ShellLayout) => mainNode(l).root as LayoutNode & { kind: 'split' };
  const column = (l: ShellLayout, i: number) => grid(l).children[i] as LayoutNode & { kind: 'split' };
  const sizeAt = (l: ShellLayout, c: number, r: number) => (column(l, c).children[r] as { size?: number }).size;

  it('arranges the panes in columns, the content first, equal sizes, lined up', () => {
    const l = arrangeGrid(four(), area);
    const g = grid(l);
    expect(g.kind).toBe('split');
    expect(g.dir).toBe('row');
    expect(g.grid).toBe(true);
    // In the panes' order as they were (top to bottom, left to right), filling columns.
    expect(g.children.map(leaves)).toEqual([
      ['content', 'g:3d#c'],
      ['g:3d#b', 'g:3d#d'],
    ]);
    expect((g.children[1] as { size?: number }).size).toBe(500);
    expect(sizeAt(l, 1, 1)).toBe(400);
  });

  it('resizing a cell moves the line across the grid: the same row in every column', () => {
    let l = arrangeGrid(four(), area);
    const cell = column(l, 0).children[1].id;
    l = resizeNode(l, cell, 300);
    expect(sizeAt(l, 0, 1)).toBe(300);
    expect(sizeAt(l, 1, 1)).toBe(300);
  });

  it('a grid whose shape breaks is no grid (a pane closed): cells resize alone again', () => {
    let l = arrangeGrid(four(), area);
    l = removeView(l, '3d#d');
    expect(grid(l).grid).toBeUndefined();
    // Three panes: two columns, a short one; not lined up.
    const three = arrangeGrid(l, area);
    expect(grid(three).grid).toBeUndefined();
    expect(grid(three).children.map(leaves)).toEqual([['content', 'g:3d#c'], ['g:3d#b']]);
  });

  it('a stored grid stays one', () => {
    const l = arrangeGrid(four(), area);
    expect(validateLayout(JSON.parse(JSON.stringify(l)))).toEqual(l);
  });
});

describe('ShellModel: the main area’s region, headless', () => {
  const view = (id: string, extra: Partial<ShellViewDef> = {}): ShellViewDef => ({
    id,
    title: id,
    icon: 'box',
    render: { mount: () => () => {} },
    ...extra,
  });
  const shell = (top?: 'overlay' | 'split' | 'both') => {
    const s = new ShellModel({ top });
    s.registerViews([
      view('home'),
      view('3d', { stage: true, home: 'stage', instances: 'many' }),
      view('seq', { home: 'top' }),
      view('plot', { home: 'stage', stageSide: 'bottom', stageSize: 200 }),
    ]);
    return s;
  };

  it('a stage-home view opens beside the content when shown, not when registered', () => {
    const s = shell();
    expect(mainNode(s.layout.value).root).toBeUndefined();
    s.show('3d', { key: 'b' });
    expect(leaves(mainNode(s.layout.value).root)).toEqual(['content', 'g:3d#b']);
    s.show('plot');
    expect(findView(s.layout.value, 'plot')!.group.size).toBe(200);
  });

  it('a top view goes where the top option says: split above (default), or overlaid on the top', () => {
    const docked = shell();
    docked.show('seq');
    expect(leaves(mainNode(docked.layout.value).root)).toEqual(['g:seq', 'content']);
    for (const mode of ['overlay', 'both'] as const) {
      const s = shell(mode);
      s.show('seq');
      expect(mainNode(s.layout.value).top?.views).toEqual(['seq']);
    }
  });

  it('split moves a placed view beside a pane; close takes it out; shown again it goes home', () => {
    const s = shell();
    s.show('home');
    s.show('3d', { key: 'b' });
    s.split('home', 'bottom', '3d#b');
    expect(leaves(mainNode(s.layout.value).root)).toEqual(['content', 'g:3d#b', 'g:home']);
    s.close('home');
    expect(findView(s.layout.value, 'home')).toBeUndefined();
    s.show('home');
    expect(findView(s.layout.value, 'home')?.group.id).toBe('right');
  });

  it('maximize: one pane alone; the commands follow the active pane; closing a maximized pane ends it', () => {
    const s = shell();
    const commands = shellCommands(s);
    const maximize = commands.find((c) => c.id === 'view.maximize-pane')!;
    const closePane = commands.find((c) => c.id === 'view.close-pane')!;
    expect(maximize.when!(undefined)).toBe('There is one pane');
    s.show('3d', { key: 'b' });
    expect(maximize.when!(undefined)).toBe(true);
    s.setActivePane('3d#b');
    maximize.run({}, undefined, new AbortController().signal);
    const node = s.paneNode('3d#b')!;
    expect(s.maximized.value).toBe(node);
    expect(closePane.when!(undefined)).toBe(true);
    closePane.run({}, undefined, new AbortController().signal);
    expect(s.maximized.value).toBeUndefined();
    expect(s.activePane.value).toBe('content');
    expect(closePane.when!(undefined)).toBe('Nothing to close');
    // Not a pane of the region: nothing.
    s.maximize('right');
    expect(s.maximized.value).toBeUndefined();
  });

  it('arranges a grid: at once with a size, or at the main area’s first size; the command needs panes', () => {
    const s = shell();
    const grid = shellCommands(s).find((c) => c.id === 'view.grid-panes')!;
    expect(grid.when!(undefined)).toBe('There is one pane');
    s.show('3d', { key: 'b' });
    s.show('3d', { key: 'c' });
    s.show('3d', { key: 'd' });
    grid.run({}, undefined, new AbortController().signal);
    const root = () => mainNode(s.layout.value).root as LayoutNode & { kind: 'split' };
    expect(root().grid).toBe(true);
    // Not measured yet: default sizes, then equal ones at the first size.
    s.setArea(1000, 800);
    expect((root().children[1] as { size?: number }).size).toBe(500);
  });

  it('a pane’s status: active, hovered, focused, maximized, several; told only when it changes; none outside', () => {
    const s = shell();
    const content = s.paneStatus('content');
    expect(content.value).toEqual({ active: true, hovered: false, focused: false, maximized: false, several: false });
    s.show('3d', { key: 'b' });
    const b = s.paneStatus('3d#b');
    let told = 0;
    const off = b.subscribe(() => told++);
    const before = told;
    s.setActivePane('3d#b');
    s.setPaneHovered('3d#b');
    s.maximize(s.paneNode('3d#b'));
    // Maximized, it shows alone: not one of several.
    expect(b.value).toEqual({ active: true, hovered: true, focused: false, maximized: true, several: false });
    expect(content.value.active).toBe(false);
    s.setPaneFocused('content');
    expect(told - before).toBe(3);
    off.unsubscribe();
    // A view outside the region has no pane status.
    s.show('home');
    expect(s.paneStatus('home').value).toMatchObject({ active: false, hovered: false });
  });

  it('a pane’s port: the layout tells its status and insets, the pane tells only its title; one port a pane', () => {
    const s = shell();
    s.show('3d', { key: 'b' });
    const port = s.panePort('3d#b');
    expect(s.panePort('3d#b')).toBe(port);
    expect(s.paneTitle('3d#b').value).toBe('3d');
    const undo = port.setTitle('Stage B');
    expect(s.paneTitle('3d#b').value).toBe('Stage B');
    undo();
    expect(s.paneTitle('3d#b').value).toBe('3d');
    port.toggleMaximized();
    expect(s.maximized.value).toBe(s.paneNode('3d#b'));
    port.toggleMaximized();
    expect(s.maximized.value).toBeUndefined();
    // The content is no view: it neither floats nor closes.
    s.panePort('content').close();
    expect(mainNode(s.layout.value).root).toBeDefined();
  });

  it('insets: the shared chrome’s widths go to the panes in the main area’s bottom corners (alone: both)', () => {
    const s = shell();
    s.show('3d', { key: 'b' });
    const content = s.panePort('content').insets;
    const b = s.panePort('3d#b').insets;
    expect(b.value).toEqual({ bottomLeft: 0, bottomRight: 0 });
    s.setOverChrome({ bottomLeft: 120, bottomRight: 96 });
    // Side by side: the content holds the bottom left, B the bottom right.
    expect(content.value).toEqual({ bottomLeft: 120, bottomRight: 0 });
    expect(b.value).toEqual({ bottomLeft: 0, bottomRight: 96 });
    s.maximize(s.paneNode('3d#b'));
    expect(b.value).toEqual({ bottomLeft: 120, bottomRight: 96 });
    expect(content.value).toEqual({ bottomLeft: 0, bottomRight: 0 });
  });

  it('the app’s pane map: what the app keeps a pane, the active one’s entry following the shell’s active pane', () => {
    const s = shell();
    s.show('3d', { key: 'b' });
    const panes = new PaneMap<{ name: string }>(s.activePane);
    const a = { name: 'A' };
    const b = { name: 'B' };
    panes.set('content', a);
    const undo = panes.set('3d#b', b);
    expect(panes.active.value).toEqual({ pane: 'content', value: a });
    s.setActivePane('3d#b');
    expect(panes.active.value).toEqual({ pane: '3d#b', value: b });
    undo();
    expect(panes.active.value).toBeUndefined();
    expect(panes.get('content')).toBe(a);
  });

  it('saves and restores the region (a stage beside the content, the top overlay)', () => {
    const s = shell('overlay');
    s.show('3d', { key: 'b' });
    s.show('seq');
    const saved = JSON.parse(JSON.stringify(s.save()));
    const t = shell('overlay');
    expect(t.restore(saved)).toBe(true);
    expect(leaves(mainNode(t.layout.value).root)).toEqual(['content', 'g:3d#b']);
    expect(mainNode(t.layout.value).top?.views).toEqual(['seq']);
  });
});
