/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import {
  columnOf,
  dockedGroups,
  findGroup,
  findView,
  inRegion,
  inRightColumn,
  type LayoutNode,
  layoutIndex,
  mainNode,
  MinGroupSize,
  railGroup,
  regionHas,
  type SplitNode,
  validateLayout,
} from '../layout.js';
import { shellCommands } from '../shell-commands.js';
import { DefaultShellLabels, keyOf, ShellModel, type ShellViewDef } from '../shell-model.js';
import { paneMenu, railMenu, viewMenu } from '../view-menu.js';

// What follows every layout write (shell §2.1): a maximize, the active pane and what is kept per pane follow the layout
// they refer to; the region's sizes keep the content a least size; the layout's index; the shell's words (shell §6).

const def = (id: string, extra: Partial<ShellViewDef> = {}): ShellViewDef => ({
  id,
  title: id[0].toUpperCase() + id.slice(1),
  icon: 'file',
  render: { mount: () => () => {} },
  ...extra,
});

const shellWith = (options = {}) => {
  const s = new ShellModel(options);
  s.registerViews([
    def('home', { place: 'head' }),
    def('stage', { title: 'Stage B', canvas: true, home: 'stage', instances: 'many' }),
    def('plot', { home: 'stage', region: { side: 'bottom' } }),
    def('panel'),
  ]);
  return s;
};

const groupOf = (s: ShellModel, ref: string) => findView(s.layout.value, ref)!.group.id;

describe('a maximize follows the layout', () => {
  it('ends when its pane leaves the region (the plot group maximized, the plot floated)', () => {
    const s = shellWith();
    s.show('stage', { key: 'b' });
    s.show('plot');
    s.maximize('plot');
    expect(s.maximized.value).toBe('plot');
    s.float('plot');
    expect(s.maximized.value).toBeUndefined();
  });

  it('ends when a view is shown or moved into the region outside it (it would stay hidden behind it)', () => {
    const s = shellWith();
    s.show('stage', { key: 'b' });
    const g1 = groupOf(s, 'stage#b');
    s.maximize('stage#b');
    s.show('plot');
    expect(findView(s.layout.value, 'plot')).toBeDefined();
    expect(s.maximized.value).toBeUndefined();
    // Into the maximized pane's group itself: it stays, on what the group shows now.
    s.maximize('stage#b');
    s.move('panel', { kind: 'group', id: g1 });
    expect(s.maximized.value).toBe('panel');
  });

  it('a window’s view maximized ends when it docks; a restore without the pane ends it', () => {
    const s = shellWith();
    s.show('stage', { key: 'b' });
    s.float('stage#b');
    s.maximize('stage#b');
    expect(s.maximized.value).toBe('stage#b');
    s.dock('stage#b');
    expect(s.maximized.value).toBeUndefined();
    s.show('plot');
    s.maximize('plot');
    const stored = s.save();
    s.close('plot');
    expect(s.maximized.value).toBeUndefined();
    s.restore(stored);
    s.maximize('plot');
    expect(s.maximized.value).toBe('plot');
    s.restore({ ...stored, root: { kind: 'main', id: 'main' } });
    expect(s.maximized.value).toBeUndefined();
  });
});

describe('the active pane follows the layout', () => {
  it('falls back to the content when its pane closes (✕, the port, the command)', () => {
    const s = shellWith();
    s.show('stage', { key: 'b' });
    s.show('plot');
    s.setActivePane('plot');
    s.close('plot');
    expect(s.activePane.value).toBe('content');
    s.setActivePane('stage#b');
    s.panePort('stage#b').close();
    expect(s.activePane.value).toBe('content');
    s.show('plot');
    s.setActivePane('plot');
    const close = shellCommands(s).find((c) => c.id === 'view.close-pane')!;
    close.run({}, undefined, new AbortController().signal);
    expect(findView(s.layout.value, 'plot')).toBeUndefined();
    expect(s.activePane.value).toBe('content');
  });

  it('stays on a pane that moves to a window, falls back when it docks out of the main area', () => {
    const s = shellWith();
    s.show('panel');
    s.split('panel', 'right');
    s.setActivePane('panel');
    s.float('panel');
    expect(s.activePane.value).toBe('panel');
    s.dock('panel');
    expect(inRegion(s.layout.value, groupOf(s, 'panel'))).toBe(false);
    expect(s.activePane.value).toBe('content');
  });

  it('maximize-pane maximizes the active pane in a window (its ref), not the content', () => {
    const s = shellWith();
    s.show('stage', { key: 'b' });
    s.float('stage#b');
    s.setActivePane('stage#b');
    const max = shellCommands(s).find((c) => c.id === 'view.maximize-pane')!;
    expect(max.when?.(undefined)).toBe(true);
    max.run({}, undefined, new AbortController().signal);
    expect(s.maximized.value).toBe('stage#b');
  });
});

describe('what is kept for a pane', () => {
  it('a closed instance in the region loses its subject and section focus', () => {
    const s = shellWith();
    s.show('stage', { key: 'b', subject: 42, section: 'x' });
    expect(s.subject('stage#b').value).toBe(42);
    s.close('stage#b');
    expect(s.subject('stage#b').value).toBeUndefined();
    expect(s.focus('stage#b').value).toBeUndefined();
  });

  it('a ref that leaves loses its title and states; a port the app holds keeps working when it comes back', () => {
    const s = shellWith();
    s.show('stage', { key: 'b' });
    const port = s.panePort('stage#b');
    port.setTitle('Stage D');
    expect(s.paneTitle('stage#b').value).toBe('Stage D');
    s.close('stage#b');
    expect(s.paneTitle('stage#b').value).toBe('Stage B');
    expect(s.panePort('stage#b')).not.toBe(port);
    expect(port.status.value.active).toBe(false);
    s.show('stage', { key: 'b' });
    s.setActivePane('stage#b');
    expect(port.status.value).toMatchObject({ active: true, several: true });
    port.toggleMaximized();
    expect(s.maximized.value).toBe('stage#b');
    // A title set before the view is placed stays (only a ref that leaves loses it).
    s.setPaneTitle('stage#c', 'Stage C');
    s.show('plot');
    s.show('stage', { key: 'c' });
    expect(s.paneTitle('stage#c').value).toBe('Stage C');
  });
});

describe('select and Show keep an instance’s key', () => {
  it('select of an unplaced instance opens that instance', () => {
    const s = shellWith();
    s.select('stage#x');
    expect(findView(s.layout.value, 'stage#x')).toBeDefined();
    expect(findView(s.layout.value, 'stage#default')).toBeUndefined();
    expect(keyOf('stage#x')).toBe('x');
    expect(keyOf('stage')).toBeUndefined();
  });

  it('the rail menu’s Show brings back that instance (hidden in a stored layout)', () => {
    const s = shellWith();
    s.show('stage', { key: 'x' });
    s.restore({ ...s.save(), hidden: ['stage#x'] });
    expect(s.layout.value.hidden).toEqual(['stage#x']);
    const show = railMenu(s, undefined).find((i) => i.id === 'show:stage#x')!;
    expect(show.label).toBe('Show Stage B');
    show.run?.();
    expect(s.layout.value.hidden).not.toContain('stage#x');
    expect(findView(s.layout.value, 'stage#default')).toBeUndefined();
  });
});

describe('sizes', () => {
  it('a stored pane keeps its size under a docked group’s least (the level’s own least applies)', () => {
    const s = shellWith();
    s.show('plot');
    const stored = JSON.parse(JSON.stringify(s.save()));
    const region = mainNode(stored).root as SplitNode;
    const plot = region.children.find((c) => c.kind === 'group') as { size?: number };
    plot.size = 96;
    const l = validateLayout(stored)!;
    expect((findView(l, 'plot')!.group as { size?: number }).size).toBe(96);
    // A docked group's least is still `MinGroupSize`.
    const small = JSON.parse(JSON.stringify(stored));
    (small.root.children[1] as { size: number }).size = 20;
    expect((validateLayout(small)!.root as SplitNode).children[1]).toMatchObject({ size: MinGroupSize });
  });

  it('a smaller main area shrinks the region’s panes so the content keeps its least size; the grid stays', () => {
    const s = shellWith();
    s.show('stage', { key: 'b' });
    s.show('stage', { key: 'c' });
    s.setArea(1200, 800);
    s.arrangeGrid(3);
    const columns = () => (mainNode(s.layout.value).root as SplitNode).children;
    const sized = () => columns().map((c) => (c as { size?: number }).size ?? 0);
    expect(sized().reduce((a, b) => a + b)).toBe(800);
    s.setArea(500, 400);
    const content = 500 - sized().reduce((a, b) => a + b);
    expect(content).toBeGreaterThanOrEqual(160);
    expect(sized().filter((x) => x > 0)).toEqual([166, 166]);
    // A lined-up grid shrinks line by line: it stays one.
    s.show('stage', { key: 'd' });
    s.setArea(1200, 800);
    s.arrangeGrid(2);
    expect((mainNode(s.layout.value).root as SplitNode).grid).toBe(true);
    s.setArea(400, 300);
    const grid = mainNode(s.layout.value).root as SplitNode;
    expect(grid.grid).toBe(true);
    expect(grid.children.map((c) => (c as { size?: number }).size)).toEqual([undefined, 200]);
    const cells = grid.children.map((c) => (c as SplitNode).children.map((x) => (x as { size?: number }).size));
    expect(cells).toEqual([
      [undefined, 150],
      [undefined, 150],
    ]);
  });

  it('resizing a pane stops where the content keeps its least size', () => {
    const s = shellWith();
    s.show('stage', { key: 'b' });
    s.setArea(1000, 600);
    const g = groupOf(s, 'stage#b');
    s.resize(g, 5000);
    expect((findGroup(s.layout.value, g)!.group as { size?: number }).size).toBe(1000 - 240);
    s.resize(g, 300);
    expect((findGroup(s.layout.value, g)!.group as { size?: number }).size).toBe(300);
  });
});

describe('the layout’s index', () => {
  it('is built once a layout value and answers as the walks did', () => {
    const s = shellWith();
    s.show('stage', { key: 'b' });
    s.show('plot');
    s.show('panel');
    s.splitBelow('home');
    s.float('stage#b');
    const l = s.layout.value;
    expect(layoutIndex(l)).toBe(layoutIndex(l));
    const walk = (n: LayoutNode, out: LayoutNode[] = []): LayoutNode[] => {
      out.push(n);
      if (n.kind === 'split') for (const c of n.children) walk(c, out);
      if (n.kind === 'main') {
        if (n.root) walk(n.root, out);
        if (n.top) walk(n.top, out);
      }
      return out;
    };
    const nodes = walk(l.root);
    for (const g of dockedGroups(l)) {
      expect(findGroup(l, g.id)!.group).toBe(g);
      for (const v of g.views) expect(findView(l, v)!.group).toBe(g);
    }
    expect(findView(l, 'stage#b')!.window).toBeDefined();
    expect(inRightColumn(l, railGroup(l)!.id)).toBe(true);
    const region = mainNode(l).root!;
    for (const n of nodes) expect(regionHas(l, n.id)).toBe(walk(region).includes(n));
    const plot = groupOf(s, 'plot');
    expect(inRegion(l, plot)).toBe(true);
    expect(inRightColumn(l, plot)).toBe(false);
    expect(columnOf(l, groupOf(s, 'home'))?.dir).toBe('column');
  });
});

describe('the shell’s words (shell §6)', () => {
  it('the app’s replace the defaults once, for the menus and the commands (headless too)', () => {
    const s = shellWith({
      labels: { float: 'Flotter', maximize: 'Agrandir', showView: (t: string) => `Montrer ${t}`, gridTitle: 'Grille' },
    });
    expect(Object.isFrozen(s.labels)).toBe(true);
    expect(s.labels.dock).toBe(DefaultShellLabels.dock);
    s.show('panel');
    expect(viewMenu(s, 'panel').find((i) => i.id === 'float')!.label).toBe('Flotter');
    s.show('plot');
    expect(paneMenu(s, 'plot').find((i) => i.id === 'maximize')!.label).toBe('Agrandir');
    s.close('panel');
    expect(railMenu(s, undefined).find((i) => i.id === 'show:panel')!.label).toBe('Montrer Panel');
    const commands = shellCommands(s);
    expect(commands.find((c) => c.id === 'view.grid-panes')!.title).toBe('Grille');
    expect(commands[0].category).toBe('View');
    expect(new ShellModel({ labels: { fullscreenUnavailable: 'Non' } }).fullscreenAvailable.value).toBe('Non');
    expect(new ShellModel({ labels: { mainTitle: 'Scène' } }).paneTitle('content').value).toBe('Scène');
  });
});
