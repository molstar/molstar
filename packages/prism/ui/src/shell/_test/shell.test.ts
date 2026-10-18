/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { ViewStateManager } from '../../managers/view-state.js';
import type { ViewRenderer } from '../../dom/view.js';
import { ShellModel, type ShellPage, type ShellTab, type ShellViewContext } from '../shell.js';

const view: ViewRenderer<ShellViewContext> = { mount: () => () => {} };
const tab = (id: string, extra: Partial<ShellTab> = {}): ShellTab => ({
  id,
  title: id,
  icon: 'house',
  render: view,
  ...extra,
});

describe('ShellModel', () => {
  it('places tabs in their docks in order; the first becomes active', () => {
    const shell = new ShellModel();
    shell.registerTabs([tab('b', { order: 2 }), tab('a', { order: 1 }), tab('seq', { dock: 'top' })]);
    expect(shell.dock('right').state.value).toMatchObject({ tabs: ['a', 'b'], active: 'a' });
    expect(shell.dock('right').active.value?.id).toBe('a');
    expect(shell.dock('top').tabs.value.map((t) => t.id)).toEqual(['seq']);
  });

  it('remembers sections left open or closed; a section request is forgotten when the tab is shown plainly', () => {
    const shell = new ShellModel();
    shell.registerTabs([tab('selection'), tab('structure')]);
    expect(shell.sectionOpen('selection', 'measure')).toBeUndefined();
    shell.setSectionOpen('selection', 'measure', true);
    expect(shell.sectionOpen('selection', 'measure')).toBe(true);

    shell.show('selection', 'measure');
    expect(shell.tabFocus.value).toEqual({ tab: 'selection', section: 'measure' });
    shell.show('structure');
    expect(shell.tabFocus.value).toEqual({ tab: 'selection', section: 'measure' });
    shell.select('selection');
    expect(shell.tabFocus.value).toBeUndefined();
  });

  it('asking again for the same tab section or card section emits again', () => {
    const shell = new ShellModel();
    shell.registerTabs([tab('selection')]);
    shell.registerOverlays([{ id: 'export', title: 'Export', render: view }]);
    const focus: unknown[] = [];
    const overlay: unknown[] = [];
    shell.tabFocus.subscribe((f) => focus.push(f));
    shell.overlay.subscribe((o) => overlay.push(o));
    shell.show('selection', 'measure');
    shell.show('selection', 'measure');
    expect(focus.filter(Boolean)).toHaveLength(2);
    shell.openOverlay('export', undefined, 'mp4');
    shell.openOverlay('export', undefined, 'mp4');
    expect(overlay.filter(Boolean)).toHaveLength(2);
  });

  it('places a tab registered later by its order among the tabs already there', () => {
    const shell = new ShellModel();
    shell.registerTabs([tab('home', { order: 0 }), tab('structure', { order: 1 }), tab('tree', { order: 2 })]);
    shell.registerTabs([tab('selection', { order: 1.2 })]);
    shell.registerTabs([tab('extra', { order: 9 })]);
    expect(shell.dock('right').state.value.tabs).toEqual(['home', 'structure', 'selection', 'tree', 'extra']);
  });

  it('selecting the active tab collapses its dock; selecting again (or another tab) expands it', () => {
    const shell = new ShellModel();
    shell.registerTabs([tab('a'), tab('b')]);
    const right = shell.dock('right');
    shell.select('a');
    expect(right.state.value.collapsed).toBe(true);
    expect(right.active.value).toBeUndefined();
    shell.select('b');
    expect(right.state.value).toMatchObject({ collapsed: false, active: 'b' });
    shell.show('a');
    expect(right.active.value?.id).toBe('a');
  });

  it('moves a tab between docks and removes tabs on unregister', () => {
    const shell = new ShellModel();
    const remove = shell.registerTabs([tab('a'), tab('b')]);
    shell.move('a', 'left');
    expect(shell.dock('right').state.value).toMatchObject({ tabs: ['b'], active: 'b' });
    expect(shell.dock('left').state.value).toMatchObject({ tabs: ['a'], active: 'a' });
    remove();
    expect(shell.dock('left').tabs.value).toEqual([]);
    expect(shell.dock('right').state.value).toMatchObject({ tabs: [], active: undefined });
  });

  it('toggles overlays and resolves the mode from the width', () => {
    const shell = new ShellModel({ compactWidth: 700 });
    const close = shell.registerOverlays([{ id: 'save', title: 'Save', render: view }]);
    shell.toggleOverlay('save');
    expect(shell.overlay.value?.overlay.id).toBe('save');
    shell.toggleOverlay('save');
    expect(shell.overlay.value).toBeUndefined();
    shell.openOverlay('save');
    close();
    expect(shell.overlay.value).toBeUndefined();

    expect(shell.mode.value).toBe('desktop');
    shell.setWidth(500);
    expect(shell.mode.value).toBe('compact');
    // Embedded is compact too where it is narrow (shell.md §3).
    shell.setMode('embedded');
    expect(shell.mode.value).toBe('compact');
    shell.setWidth(900);
    expect(shell.mode.value).toBe('embedded');
  });

  it('gives embedded a layout of its own: the compact rail, a floating right panel, nothing pinned', () => {
    const viewState = new ViewStateManager({ docks: { right: { floating: false, pinned: true } } });
    const shell = new ShellModel({ viewState });
    shell.registerTabs([tab('a'), tab('b')]);
    shell.show('a');
    expect(shell.compactRail.value).toBe(false);
    expect(shell.dock('right').state.value).toMatchObject({ floating: false, pinned: true });

    shell.setMode('embedded');
    expect(shell.mode.value).toBe('embedded');
    expect(shell.compactRail.value).toBe(true);
    expect(shell.railShown.value).toBe(false);
    // The rail starts hidden, and what opened from it closed.
    expect(shell.dock('right').active.value).toBeUndefined();
    expect(shell.dock('right').state.value).toMatchObject({ floating: true, pinned: false });
    shell.toggleRail();
    expect(shell.railShown.value).toBe(true);
    // The layout keeps its own values for another mode.
    expect(viewState.layout.value.docks.right).toMatchObject({ floating: false, pinned: true });
    expect(viewState.layout.value.compactRail).toBe(false);

    shell.setMode('desktop');
    expect(shell.compactRail.value).toBe(false);
    expect(shell.dock('right').state.value).toMatchObject({ floating: false, pinned: true });
  });

  it('pins the right panel by default on the desktop only: compact reads it as unpinned and keeps the pin', () => {
    const viewState = new ViewStateManager();
    const shell = new ShellModel({ viewState });
    shell.registerTabs([tab('a')]);
    shell.show('a');
    expect(shell.dock('right').state.value.pinned).toBe(true);
    // A phone: the panel is a bottom sheet, never pinned (Escape and a click outside close it).
    shell.setMode('compact');
    expect(shell.dock('right').state.value.pinned).toBe(false);
    expect(viewState.layout.value.docks.right.pinned).toBe(true);
    shell.setMode('desktop');
    expect(shell.dock('right').state.value.pinned).toBe(true);
  });

  it('pins nothing beside a compact rail, and keeps the pin for the full rail', () => {
    const viewState = new ViewStateManager();
    const shell = new ShellModel({ viewState });
    shell.registerTabs([tab('a')]);
    shell.show('a');
    shell.setCompactRail(true);
    expect(shell.dock('right').state.value.pinned).toBe(false);
    expect(viewState.layout.value.docks.right.pinned).toBe(true);
    shell.setCompactRail(false);
    expect(shell.dock('right').state.value.pinned).toBe(true);
  });

  it('hides a tab without unregistering it, and shows it again where it was', () => {
    const shell = new ShellModel();
    shell.registerTabs([tab('a', { order: 1 }), tab('b', { order: 2 }), tab('c', { order: 3 })]);
    shell.show('b');
    shell.hideTab('b');
    expect(shell.dock('right').tabs.value.map((t) => t.id)).toEqual(['a', 'c']);
    // Its open panel closed; it keeps its place in the layout, which keeps it hidden.
    expect(shell.dock('right').active.value).toBeUndefined();
    expect(shell.dock('right').state.value.tabs).toEqual(['a', 'b', 'c']);
    expect(shell.viewState.layout.value.hiddenTabs).toEqual(['b']);
    expect(shell.hiddenTabs.value.map((t) => t.id)).toEqual(['b']);

    shell.showTab('b');
    expect(shell.dock('right').tabs.value.map((t) => t.id)).toEqual(['a', 'b', 'c']);
    expect(shell.hiddenTabs.value).toEqual([]);

    // Shown by name (a command), a hidden tab comes back open.
    shell.hideTab('c');
    shell.show('c');
    expect(shell.dock('right').active.value?.id).toBe('c');
    expect(shell.hiddenTabs.value).toEqual([]);
  });

  it('resets the layout to its seed, putting registered tabs back and showing hidden ones', () => {
    const viewState = new ViewStateManager({ docks: { right: { size: 400 } } });
    const shell = new ShellModel({ viewState });
    shell.registerTabs([tab('a', { order: 1 }), tab('b', { order: 2 }), tab('seq', { dock: 'top' })]);
    shell.show('b');
    shell.move('a', 'left');
    shell.setSize('right', 600);
    shell.setPinned('right', false);
    shell.hideTab('seq');
    shell.resetLayout();
    expect(viewState.layout.value.docks.right).toMatchObject({ tabs: ['a', 'b'], size: 400, pinned: true });
    expect(viewState.layout.value.docks.left.tabs).toEqual([]);
    expect(viewState.layout.value.hiddenTabs).toEqual([]);
    expect(shell.dock('top').tabs.value.map((t) => t.id)).toEqual(['seq']);
    // What is open stays open where it still is; a dock left without its open tab has none.
    expect(shell.dock('right').active.value?.id).toBe('b');
    expect(shell.dock('left').state.value.active).toBeUndefined();
  });

  it('floats the right and top docks by default and can dock them', () => {
    const shell = new ShellModel();
    expect(shell.dock('right').state.value.floating).toBe(true);
    expect(shell.dock('top').state.value.floating).toBe(true);
    expect(shell.dock('left').state.value.floating).toBe(false);
    shell.setFloating('right', false);
    expect(shell.dock('right').state.value.floating).toBe(false);
  });

  it('clamps dock sizes and keeps the layout in the view state', () => {
    const viewState = new ViewStateManager();
    const shell = new ShellModel({ viewState });
    shell.setSize('right', 20);
    expect(shell.dock('right').state.value.size).toBe(180);
    expect(viewState.layout.value.docks.right.size).toBe(180);
    expect(viewState.layout.value.mode).toBe('auto');
  });

  it('shows toolbar buttons as the spec says, and changes them at runtime', () => {
    const shell = new ShellModel({ buttons: { xr: false } });
    shell.registerToolbar([
      { id: 'reset', icon: 'focus', title: 'Reset camera' },
      { id: 'xr', icon: 'glasses', title: 'Enter XR' },
    ]);
    expect(shell.toolbar.value.map((i) => i.id)).toEqual(['reset']);
    shell.setButtons({ xr: true, reset: false });
    expect(shell.toolbar.value.map((i) => i.id)).toEqual(['xr']);
  });
});

describe('pages', () => {
  const page = (id: string): ShellPage => ({ id, title: id, render: view });

  it('opens a registered page at a focus, closes an open card, and closes', () => {
    const shell = new ShellModel();
    shell.registerOverlays([{ id: 'export', title: 'Export', render: view }]);
    shell.registerPages([page('inspector')]);
    shell.openOverlay('export');
    shell.openPage('inspector', 'tools');
    expect(shell.page.value?.page.id).toBe('inspector');
    expect(shell.page.value?.focus).toBe('tools');
    expect(shell.overlay.value).toBeUndefined();
    shell.closePage();
    expect(shell.page.value).toBeUndefined();
  });

  it('ignores an unknown page, and closes a page when it is unregistered', () => {
    const shell = new ShellModel();
    shell.openPage('nothing');
    expect(shell.page.value).toBeUndefined();
    const dispose = shell.registerPages([page('inspector')]);
    shell.openPage('inspector');
    dispose();
    expect(shell.page.value).toBeUndefined();
  });

  it('hides the rail of a compact rail until it is opened, and closes the panel and the card with it', () => {
    const shell = new ShellModel({ compactRail: true });
    shell.registerTabs([tab('a')]);
    expect(shell.railShown.value).toBe(false);
    shell.toggleRail();
    expect(shell.railShown.value).toBe(true);
    shell.show('a');
    shell.setPinned('right', true);
    shell.toggleRail();
    expect(shell.railShown.value).toBe(false);
    expect(shell.dock('right').active.value).toBeUndefined();
    // Not compact: always shown, the toggle does nothing.
    const plain = new ShellModel();
    plain.toggleRail();
    expect(plain.railShown.value).toBe(true);
  });

  it('keeps a tab open when it is registered again before the old one goes (a lazy stub replaced by the real tab)', () => {
    const shell = new ShellModel();
    shell.registerTabs([tab('home')]);
    const removeStub = shell.registerTabs([tab('tree')]);
    shell.show('tree');
    const seen: (string | undefined)[] = [];
    const sub = shell.dock('right').active.subscribe((t) => seen.push(t?.id));
    // The real tab comes first, then the stub goes: the panel never closes or falls back to another tab.
    shell.registerTabs([tab('tree', { title: 'State tree' })]);
    removeStub();
    sub.unsubscribe();
    expect(seen.every((id) => id === 'tree')).toBe(true);
    expect(shell.dock('right').active.value?.title).toBe('State tree');
  });

  it('closes the panel of a compact rail before it hides the rail (hidden UI is paused: no flash on showing it again)', () => {
    const shell = new ShellModel({ compactRail: true });
    shell.registerTabs([tab('home'), tab('structure')]);
    shell.setRailOpen(true);
    shell.show('structure');
    const seen: string[] = [];
    const a = shell.railShown.subscribe((shown) => seen.push(`rail:${shown}`));
    const b = shell.dock('right').active.subscribe((t) => seen.push(`panel:${t ? t.id : 'none'}`));
    seen.length = 0;
    shell.setRailOpen(false);
    a.unsubscribe();
    b.unsubscribe();
    expect(seen.indexOf('panel:none')).toBeGreaterThanOrEqual(0);
    expect(seen.indexOf('panel:none')).toBeLessThan(seen.indexOf('rail:false'));
  });
});
