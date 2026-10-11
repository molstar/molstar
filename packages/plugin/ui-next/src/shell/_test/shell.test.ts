/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { ViewRenderer } from '../../dom/view.js';
import { ShellModel, type ShellTab, type ShellViewContext } from '../shell.js';

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
    shell.setMode('embedded');
    expect(shell.mode.value).toBe('embedded');
  });

  it('floats the right dock by default and can dock it', () => {
    const shell = new ShellModel();
    expect(shell.dock('right').state.value.floating).toBe(true);
    expect(shell.dock('left').state.value.floating).toBe(false);
    shell.setFloating('right', false);
    expect(shell.layout.docks.right.floating).toBe(false);
  });

  it('clamps dock sizes and reports the layout', () => {
    const shell = new ShellModel();
    shell.setSize('right', 20);
    expect(shell.layout.docks.right.size).toBe(180);
    expect(shell.layout.mode).toBe('auto');
  });
});
