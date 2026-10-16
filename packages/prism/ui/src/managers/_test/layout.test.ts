/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { domView } from '../../dom/view.js';
import { ShellModel } from '../../shell/shell.js';
import { NavigationManager } from '../navigation.js';
import { ViewStateManager } from '../view-state.js';

const tab = (id: string, order = 0) => ({ id, title: id, icon: 'square' as const, order, render: domView(() => null) });

describe('viewState and navigation', () => {
  it('hold the shell state, and the shell follows them', () => {
    const viewState = new ViewStateManager();
    const navigation = new NavigationManager();
    const shell = new ShellModel({ viewState, navigation });
    shell.registerTabs([tab('home'), tab('structure', 1)]);
    expect(viewState.layout.value.docks.right.tabs).toEqual(['home', 'structure']);
    expect(navigation.tab.value).toBe('home');

    shell.show('structure', 'components');
    expect(navigation.tab.value).toBe('structure');
    expect(navigation.section.value).toBe('components');

    // Writing a manager moves the shell.
    navigation.setDock('right', { collapsed: true });
    expect(shell.dock('right').active.value).toBeUndefined();
    viewState.setDock('right', { size: 500 });
    expect(shell.dock('right').state.value.size).toBe(500);
  });

  it('save and restore, dropping other versions', () => {
    const a = new ShellModel();
    a.registerTabs([tab('home'), tab('structure', 1)]);
    a.show('structure');
    a.setSize('right', 420);
    const layout = a.viewState.save();
    const nav = a.navigation.save();

    const viewState = new ViewStateManager();
    const navigation = new NavigationManager();
    expect(viewState.restore(JSON.parse(JSON.stringify(layout)))).toBe(true);
    expect(navigation.restore(JSON.parse(JSON.stringify(nav)))).toBe(true);
    const b = new ShellModel({ viewState, navigation });
    b.registerTabs([tab('home'), tab('structure', 1)]);
    expect(b.dock('right').state.value.size).toBe(420);
    expect(b.dock('right').active.value?.id).toBe('structure');

    expect(viewState.restore({ version: 2, layout: {} })).toBe(false);
    expect(navigation.restore('nonsense')).toBe(false);
  });

  it('keeps the minimum dock size', () => {
    const v = new ViewStateManager();
    v.setDock('right', { size: 3 });
    expect(v.layout.value.docks.right.size).toBe(180);
  });

  it('restores only a layout with known modes and finite sizes, and keeps the size floor', () => {
    const v = new ViewStateManager();
    const stored = (layout: object) => ({ version: 1, layout });
    expect(v.restore(stored({ mode: 'foo', docks: {} }))).toBe(false);
    expect(v.restore(stored({ docks: { right: { tabs: [], size: null } } }))).toBe(false);
    expect(v.restore(stored({ docks: { right: { tabs: [], size: Number.NaN } } }))).toBe(false);
    expect(v.restore(stored({ docks: { right: { tabs: [], size: 300, pinned: 'yes' } } }))).toBe(false);
    expect(v.layout.value.mode).toBe('auto');
    expect(v.restore(stored({ mode: 'compact', docks: { right: { tabs: ['a'], size: 5 } } }))).toBe(true);
    expect(v.layout.value.mode).toBe('compact');
    expect(v.layout.value.docks.right.size).toBe(180);
  });

  it('keeps hidden tabs in a saved layout, and drops a malformed list', () => {
    const v = new ViewStateManager();
    v.setTabHidden('history', true);
    v.setTabHidden('history', true);
    const saved = JSON.parse(JSON.stringify(v.save()));
    const w = new ViewStateManager();
    expect(w.restore(saved)).toBe(true);
    expect(w.layout.value.hiddenTabs).toEqual(['history']);
    expect(w.restore({ version: 1, layout: { docks: {}, hiddenTabs: [1] } })).toBe(false);
    // A layout kept before hidden tabs existed restores with none.
    const old = new ViewStateManager();
    expect(old.restore({ version: 1, layout: { docks: { right: { tabs: ['x'] } } } })).toBe(true);
    expect(old.layout.value.hiddenTabs).toEqual([]);
    expect(old.layout.value.docks.right.tabs).toEqual(['x']);
  });

  it('resets to the layout it started with', () => {
    const v = new ViewStateManager({ compactRail: true, docks: { right: { size: 400 } }, buttons: { xr: false } });
    v.setDock('right', { size: 500, pinned: true });
    v.setTabHidden('a', true);
    v.setButtons({ xr: true });
    v.reset();
    expect(v.layout.value).toEqual(v.initial);
    expect(v.layout.value.docks.right).toMatchObject({ size: 400, pinned: false });
    expect(v.layout.value.buttons).toEqual({ xr: false });
    expect(v.layout.value.hiddenTabs).toEqual([]);
    expect(v.layout.value.compactRail).toBe(true);
  });
});
