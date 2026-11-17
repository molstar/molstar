/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import {
  allGroups,
  clampRect,
  columnOf,
  defaultLayout,
  dockedGroups,
  findView,
  type LayoutNode,
  type LayoutTarget,
  MinGroupSize,
  MinWindow,
  mergeGroup,
  moveInColumn,
  moveView,
  placedViews,
  placeWindow,
  RightGroup,
  railGroup,
  raiseWindow,
  removeView,
  resizeNode,
  rightColumn,
  selectView,
  setCollapsed,
  setHidden,
  type ShellLayout,
  stacked,
  validateLayout,
} from '../layout.js';

/** Every invariant `normalize` promises (layout.ts). */
function checkInvariants(l: ShellLayout) {
  let mains = 0;
  const ids = new Set<string>();
  const views = new Set<string>();
  const walk = (n: LayoutNode, parentDir?: 'row' | 'column') => {
    expect(ids.has(n.id), `duplicate id ${n.id}`).toBe(false);
    ids.add(n.id);
    if (n.kind === 'main') mains++;
    else if (n.kind === 'group') {
      expect(n.views.length).toBeGreaterThan(0);
      for (const v of n.views) {
        expect(views.has(v), `view ${v} twice`).toBe(false);
        views.add(v);
      }
      const shown = n.views.filter((v) => !l.hidden.includes(v));
      if (shown.length) expect(shown).toContain(n.active);
    } else {
      expect(n.children.length).toBeGreaterThanOrEqual(2);
      expect(n.dir).not.toBe(parentDir);
      // One child takes the rest: the one holding the main area, else the first.
      const unsized = n.children.filter((c) => c.kind === 'main' || c.size === undefined);
      expect(unsized.length).toBeGreaterThanOrEqual(1);
      for (const c of n.children) walk(c, n.dir);
    }
  };
  walk(l.root);
  expect(mains).toBe(1);
  for (const w of l.windows) {
    expect(ids.has(w.group.id)).toBe(false);
    ids.add(w.group.id);
    expect(w.group.views.length).toBeGreaterThan(0);
    for (const v of w.group.views) {
      expect(views.has(v)).toBe(false);
      views.add(v);
    }
    expect(w.rect.width).toBeGreaterThanOrEqual(MinWindow.width);
    expect(w.rect.height).toBeGreaterThanOrEqual(MinWindow.height);
  }
  expect(stacked(l).map((w) => w.z)).toEqual(l.windows.map((_, i) => i + 1));
  for (const h of l.hidden) expect(views.has(h)).toBe(true);
}

const base = () => defaultLayout(['home', 'layers', 'settings'], { active: 'home' });

describe('the shell layout', () => {
  it('starts as the main area and the right group beside it, 340 wide, its first view active', () => {
    const l = defaultLayout(['home', 'layers']);
    checkInvariants(l);
    expect(l.root).toEqual({
      kind: 'split',
      id: 'root',
      dir: 'row',
      children: [
        { kind: 'main', id: 'main' },
        { kind: 'group', id: RightGroup, views: ['home', 'layers'], active: 'home', size: 340 },
      ],
    });
    expect(railGroup(l)?.id).toBe(RightGroup);
  });

  it('floats a view out and docks it back; a group left empty goes', () => {
    let l = moveView(base(), 'layers', { kind: 'float', rect: { x: 10, y: 10, width: 300, height: 200 } });
    checkInvariants(l);
    expect(findView(l, 'layers')?.window).toBeDefined();
    expect(dockedGroups(l)[0].views).toEqual(['home', 'settings']);
    l = moveView(l, 'layers', { kind: 'group', id: RightGroup, index: 1 });
    checkInvariants(l);
    expect(l.windows).toEqual([]);
    expect(dockedGroups(l)[0]).toMatchObject({ views: ['home', 'layers', 'settings'], active: 'layers' });
  });

  it('splits the right column below a group: the upper one keeps the rail and the rest, the lower one a size', () => {
    let l = moveView(base(), 'settings', { kind: 'edge', of: RightGroup, side: 'bottom', size: 200 });
    checkInvariants(l);
    const column = rightColumn(l)!;
    expect(column).toMatchObject({ kind: 'split', dir: 'column', size: 340 });
    const [upper, lower] = (column as any).children;
    expect(upper).toMatchObject({ id: RightGroup, views: ['home', 'layers'] });
    expect(upper.size).toBeUndefined();
    expect(lower).toMatchObject({ views: ['settings'], size: 200 });
    expect(railGroup(l)?.id).toBe(RightGroup);
    expect(columnOf(l, lower.id)).toBe(column);
    // Move up: back into the upper group; the split collapses and the right group gets its width back.
    l = moveInColumn(l, 'settings', -1);
    checkInvariants(l);
    expect(rightColumn(l)).toMatchObject({ kind: 'group', id: RightGroup, size: 340 });
  });

  it('a window holds several views as tabs, and merges into another group', () => {
    let l = moveView(base(), 'layers', { kind: 'float', rect: { x: 0, y: 0, width: 300, height: 200 }, id: 'w' });
    l = moveView(l, 'settings', { kind: 'group', id: 'w' });
    checkInvariants(l);
    expect(l.windows[0].group).toMatchObject({ views: ['layers', 'settings'], active: 'settings' });
    l = mergeGroup(l, 'w', RightGroup);
    checkInvariants(l);
    expect(l.windows).toEqual([]);
    expect(dockedGroups(l)[0]).toMatchObject({ views: ['home', 'layers', 'settings'], active: 'settings' });
  });

  it('raises a window to the top of the stack', () => {
    let l = moveView(base(), 'layers', { kind: 'float', rect: { x: 0, y: 0, width: 300, height: 200 }, id: 'a' });
    l = moveView(l, 'settings', { kind: 'float', rect: { x: 30, y: 30, width: 300, height: 200 }, id: 'b' });
    expect(stacked(l).map((w) => w.group.id)).toEqual(['a', 'b']);
    l = raiseWindow(l, 'a');
    checkInvariants(l);
    expect(stacked(l).map((w) => w.group.id)).toEqual(['b', 'a']);
  });

  it('a hidden view keeps its place, another becomes active; selecting it shows it again', () => {
    let l = setHidden(base(), 'home', true);
    checkInvariants(l);
    expect(dockedGroups(l)[0]).toMatchObject({ views: ['home', 'layers', 'settings'], active: 'layers' });
    l = selectView(l, 'home');
    expect(l.hidden).toEqual([]);
    expect(dockedGroups(l)[0].active).toBe('home');
  });

  it('removes a view; the last one out takes its group, and the main area stays', () => {
    let l = removeView(base(), 'home');
    l = removeView(l, 'layers');
    l = removeView(l, 'settings');
    checkInvariants(l);
    expect(l.root).toEqual({ kind: 'main', id: 'main' });
    expect(railGroup(l)).toBeUndefined();
  });

  it('resizes a sized node (at least the minimum); collapses a group', () => {
    let l = resizeNode(base(), RightGroup, 20);
    expect((rightColumn(l) as any).size).toBe(MinGroupSize);
    l = resizeNode(l, RightGroup, 412.6);
    expect((rightColumn(l) as any).size).toBe(413);
    l = setCollapsed(l, RightGroup, true);
    expect(dockedGroups(l)[0].collapsed).toBe(true);
    l = selectView(l, 'layers');
    expect(dockedGroups(l)[0].collapsed).toBeUndefined();
  });

  it('keeps its invariants over random sequences of operations', () => {
    let seed = 7;
    const rand = (n: number) => {
      seed = (seed * 1103515245 + 12345) & 0x7fffffff;
      return seed % n;
    };
    const refs = ['a', 'b', 'c', 'd', 'e', 'f#1', 'f#2'];
    for (let run = 0; run < 40; run++) {
      let l = defaultLayout(refs.slice(0, 3));
      for (let step = 0; step < 60; step++) {
        const ref = refs[rand(refs.length)];
        const groups = allGroups(l);
        const g = groups[rand(Math.max(1, groups.length))];
        const op = rand(9);
        const targets: LayoutTarget[] = [
          { kind: 'group', id: g?.id ?? RightGroup, index: rand(4) },
          {
            kind: 'edge',
            of: rand(2) ? 'main' : (g?.id ?? 'main'),
            side: (['left', 'right', 'top', 'bottom'] as const)[rand(4)],
            size: 150 + rand(300),
          },
          { kind: 'float', rect: { x: rand(500), y: rand(400), width: rand(600), height: rand(500) } },
        ];
        if (op < 4) l = moveView(l, ref, targets[rand(3)]);
        else if (op === 4) l = removeView(l, ref);
        else if (op === 5) l = setHidden(l, ref, !!rand(2));
        else if (op === 6) l = selectView(l, ref);
        else if (op === 7 && g) l = mergeGroup(l, g.id, groups[rand(groups.length)].id);
        else if (g) l = raiseWindow(l, g.id);
        checkInvariants(l);
        // A stored layout comes back as itself.
        expect(validateLayout(JSON.parse(JSON.stringify(l)))).toEqual(l);
      }
      expect(placedViews(l).length).toBeLessThanOrEqual(refs.length);
    }
  });
});

describe('stored layouts', () => {
  it('drops another version or a wrong shape', () => {
    expect(validateLayout(undefined)).toBeUndefined();
    expect(validateLayout({ ...base(), version: 2 })).toBeUndefined();
    expect(validateLayout({ ...base(), root: { kind: 'nope' } })).toBeUndefined();
  });

  it('keeps unknown view ids, puts the main area back, clamps sizes, drops broken windows', () => {
    const stored = {
      version: 1,
      root: { kind: 'group', id: 'right', views: ['later', 3], size: 5 },
      windows: [{ group: { kind: 'group', id: 'w', views: ['x'] }, rect: { x: 0, y: 0, width: 'wide', height: 1 } }],
      hidden: ['later', 'gone'],
      rail: false,
      seq: 4,
    };
    const l = validateLayout(stored)!;
    checkInvariants(l);
    expect(l.root).toMatchObject({ kind: 'split', dir: 'row' });
    expect(dockedGroups(l)[0]).toMatchObject({ id: 'right', views: ['later'], size: MinGroupSize });
    expect(l.windows).toEqual([]);
    expect(l.hidden).toEqual(['later']);
    expect(l.rail).toBe(false);
  });
});

describe('window placement', () => {
  const area = { width: 800, height: 600 };
  it('clamps a window into the area, no larger than it', () => {
    expect(clampRect({ x: 700, y: -20, width: 300, height: 200 }, area)).toEqual({
      x: 500,
      y: 0,
      width: 300,
      height: 200,
    });
    expect(clampRect({ x: 10, y: 10, width: 2000, height: 50 }, area)).toEqual({
      x: 0,
      y: 10,
      width: 800,
      height: MinWindow.height,
    });
  });

  it('opens beside its anchor, else cascades past the windows already there', () => {
    const size = { width: 300, height: 200 };
    expect(placeWindow(size, area, [], { x: 700, y: 100, width: 32, height: 32 })).toEqual({
      x: 392,
      y: 100,
      width: 300,
      height: 200,
    });
    const first = placeWindow(size, area, []);
    const second = placeWindow(size, area, [first]);
    expect(second.x - first.x).toBe(24);
    expect(second.y - first.y).toBe(24);
  });
});
