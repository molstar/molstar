/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { dockedGroups, findView, RightGroup, railGroup, rightColumn, stacked } from '../layout.js';
import { ShellModel, type ShellViewDef } from '../shell-model.js';

const view = (id: string, extra: Partial<ShellViewDef> = {}): ShellViewDef => ({
  id,
  title: id,
  icon: 'box',
  render: { mount: () => () => {} },
  ...extra,
});

function shell() {
  const s = new ShellModel();
  s.setArea(1000, 700);
  s.registerViews([
    view('home', { order: 0, place: 'head' }),
    view('settings', { order: 10, place: 'foot' }),
    view('layers', { order: 1 }),
    view('export', { home: 'card' }),
    view('inspector', { home: 'page' }),
    view('confirm', { home: 'modal' }),
    view('plot', { home: 'float', instances: 'many', size: { width: 400, height: 300 } }),
  ]);
  return s;
}

describe('ShellModel', () => {
  it('places docked views in their home group by order, folded until one is shown', () => {
    const s = shell();
    const right = railGroup(s.layout.value)!;
    expect(right.id).toBe(RightGroup);
    expect(right.views).toEqual(['home', 'layers', 'settings']);
    expect(right.collapsed).toBe(true);
    s.show('layers');
    expect(railGroup(s.layout.value)).toMatchObject({ active: 'layers' });
    expect(railGroup(s.layout.value)!.collapsed).toBeUndefined();
  });

  it('a press on the open view folds the group; another press opens it', () => {
    const s = shell();
    s.show('home');
    s.select('home');
    expect(railGroup(s.layout.value)!.collapsed).toBe(true);
    s.select('home');
    expect(railGroup(s.layout.value)!.collapsed).toBeUndefined();
    s.select('layers');
    expect(railGroup(s.layout.value)!.active).toBe('layers');
  });

  it('cards, pages and modals open in their layers; a page closes a card, a docked view leaves a card', () => {
    const s = shell();
    s.show('export', { section: 'image' });
    expect(s.card.value?.ref).toBe('export');
    expect(s.focus('export').value).toBe('image');
    s.show('inspector');
    expect(s.card.value).toBeUndefined();
    expect(s.page.value?.ref).toBe('inspector');
    s.show('export');
    s.show('layers');
    expect(s.card.value).toBeUndefined();
    s.show('confirm');
    expect(s.modal.value?.ref).toBe('confirm');
    s.close('confirm');
    expect(s.modal.value).toBeUndefined();
  });

  it('instances open as windows by key, each with its own subject, cascaded; closing one leaves the other', () => {
    const s = shell();
    s.show('plot', { key: 'a', subject: 'structure A' });
    s.show('plot', { key: 'b', subject: 'structure B' });
    const windows = stacked(s.layout.value);
    expect(windows.map((w) => w.group.views[0])).toEqual(['plot#a', 'plot#b']);
    expect(windows[0].rect).toMatchObject({ width: 400, height: 300 });
    expect(windows[1].rect.x - windows[0].rect.x).toBe(24);
    expect(s.subject('plot#b').value).toBe('structure B');
    s.show('plot', { key: 'a' });
    expect(stacked(s.layout.value).map((w) => w.group.views[0])).toEqual(['plot#b', 'plot#a']);
    s.close('plot#a');
    expect(findView(s.layout.value, 'plot#a')).toBeUndefined();
    expect(findView(s.layout.value, 'plot#b')).toBeDefined();
  });

  it('pops a docked view out and docks it back; closing a popped-out view hides it at home', () => {
    const s = shell();
    s.show('layers');
    s.float('layers');
    expect(findView(s.layout.value, 'layers')?.window).toBeDefined();
    s.dock('layers');
    expect(findView(s.layout.value, 'layers')?.group.id).toBe(RightGroup);
    s.float('layers');
    s.close('layers');
    expect(findView(s.layout.value, 'layers')?.group.id).toBe(RightGroup);
    expect(s.layout.value.hidden).toEqual(['layers']);
    s.show('layers');
    expect(s.layout.value.hidden).toEqual([]);
  });

  it('splits below, moves up and down in the column', () => {
    const s = shell();
    s.show('layers');
    s.splitBelow('settings');
    const column = rightColumn(s.layout.value)!;
    expect(column.kind).toBe('split');
    const groups = dockedGroups(s.layout.value);
    expect(groups.map((g) => g.views)).toEqual([['home', 'layers'], ['settings']]);
    s.moveInColumn('home', 1);
    expect(dockedGroups(s.layout.value).map((g) => g.views)).toEqual([['layers'], ['settings', 'home']]);
    expect(railGroup(s.layout.value)!.views).toEqual(['layers']);
  });

  it('keeps windows inside the main area as it shrinks', () => {
    const s = shell();
    s.show('plot', { key: 'a' });
    const id = stacked(s.layout.value)[0].group.id;
    s.setWindowRect(id, { x: 700, y: 500, width: 400, height: 300 });
    expect(stacked(s.layout.value)[0].rect).toMatchObject({ x: 600, y: 400 });
    s.setArea(500, 400);
    expect(stacked(s.layout.value)[0].rect).toEqual({ x: 100, y: 100, width: 400, height: 300 });
  });

  it('saves and restores the layout; a stored layout keeps where views were, new views go home', () => {
    const s = shell();
    s.show('layers');
    s.float('layers');
    const saved = JSON.parse(JSON.stringify(s.save()));
    const t = new ShellModel();
    expect(t.restore(saved)).toBe(true);
    t.registerViews([view('home'), view('layers'), view('settings'), view('help')]);
    expect(findView(t.layout.value, 'layers')?.window).toBeDefined();
    expect(findView(t.layout.value, 'help')?.group.id).toBe(RightGroup);
    expect(t.restore({ version: 99 })).toBe(false);
  });

  it('a view that is not remembered is left out of the saved layout', () => {
    const s = new ShellModel();
    s.registerViews([view('home'), view('scratch', { home: 'float', remember: false })]);
    s.show('scratch');
    expect(findView(s.layout.value, 'scratch')).toBeDefined();
    expect(findView(s.save(), 'scratch')).toBeUndefined();
  });

  it('resets the layout: every view home, nothing hidden, no windows', () => {
    const s = shell();
    s.close('home');
    s.show('plot', { key: 'a' });
    s.splitBelow('settings');
    s.resetLayout();
    expect(s.layout.value.windows).toEqual([]);
    expect(s.layout.value.hidden).toEqual([]);
    expect(dockedGroups(s.layout.value).map((g) => g.views)).toEqual([['home', 'layers', 'settings']]);
  });

  it('modes: auto follows the width; embedded starts with the rail hidden', () => {
    const s = shell();
    s.setWidth(1200);
    expect(s.mode.value).toBe('desktop');
    s.setWidth(600);
    expect(s.mode.value).toBe('compact');
    const e = new ShellModel({ mode: 'embedded' });
    expect(e.mode.value).toBe('embedded');
    expect(e.railShown.value).toBe(false);
  });

  it('an unregistered view keeps its place; its layers close', () => {
    const s = shell();
    const undo = s.registerViews([view('extra', { order: 5 }), view('extra-card', { home: 'card' })]);
    s.show('extra-card');
    undo();
    expect(s.view('extra')).toBeUndefined();
    expect(s.card.value).toBeUndefined();
    expect(findView(s.layout.value, 'extra')).toBeDefined();
  });

  it('sections by slot and order; rail items', () => {
    const s = shell();
    const r = { mount: () => () => {} };
    const undo = s.registerSections([
      { id: 'b', slot: 'export', title: 'B', order: 2, render: r },
      { id: 'a', slot: 'export', title: 'A', order: 1, render: r },
      { id: 'c', slot: 'other', title: 'C', render: r },
    ]);
    expect(s.sections('export').value.map((x) => x.id)).toEqual(['a', 'b']);
    undo();
    expect(s.sections('export').value).toEqual([]);
  });
});

describe('ShellModel: Expand and Fullscreen, headless', () => {
  it('expands and comes back; showing a docked view leaves it; fullscreen is unavailable without a view', () => {
    const s = shell();
    s.toggleExpanded(true);
    expect(s.expanded.value).toBe(true);
    s.show('layers');
    expect(s.expanded.value).toBe(false);
    s.toggleExpanded();
    expect(s.expanded.value).toBe(true);
    s.toggleExpanded(false);
    expect(s.expanded.value).toBe(false);
    expect(typeof s.fullscreenAvailable.value).toBe('string');
    s.requestFullscreen(true);
    expect(s.fullscreenRequest.value).toBeUndefined();
  });
});
