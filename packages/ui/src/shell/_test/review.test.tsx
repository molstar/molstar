/** @vitest-environment jsdom */
/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { State } from '@molstar/core/reactive';
import { afterEach, describe, expect, it } from 'vitest';
import { flush } from '../../dom/index.js';
import { findView, MinGroupSize, MinMainSize, RailWidth, railGroup, rightColumn, validateLayout } from '../layout.js';
import { key, microtasks, mountShell, probeView } from './fixture.js';

// Regressions from an independent review of the shell (each `it` names what went wrong).

afterEach(() => {
  document.body.innerHTML = '';
});

const views = () => [
  probeView('home', { place: 'head', order: 0, compactBar: true }),
  probeView('layers', { order: 1, compactBar: true }),
  probeView('inspect', { order: 2 }),
  probeView('settings', { place: 'foot', order: 10 }),
];

describe('shell review regressions', () => {
  it('raising a window does not move its node (a move drops the focus in it)', () => {
    const m = mountShell(views());
    m.shell.float('home');
    m.shell.float('layers');
    m.shell.float('inspect');
    flush();
    const nodes = () => m.qa('[data-window]');
    const before = nodes();
    for (const w of before) m.shell.raise(w.dataset.window!);
    flush();
    expect(nodes()).toEqual(before);
    m.dispose();
  });

  it('a docked group whose views all go away is hidden (visibility and unregistering are followed)', () => {
    const shown = new State(true);
    const m = mountShell([probeView('a', { visible: shown }), probeView('b', { visible: shown })]);
    m.shell.show('a');
    flush();
    expect(m.groupNode('right')?.getAttribute('data-shown')).toBeNull();
    shown.set(false);
    flush();
    expect(m.groupNode('right')?.getAttribute('data-shown')).toBe('false');
    m.dispose();
  });

  it('a view registered anew under its id gets a new body (not the old definition’s)', () => {
    const m = mountShell([probeView('home')]);
    const extra = probeView('extra');
    const undo = m.shell.registerViews([extra.def]);
    m.shell.show('extra');
    flush();
    m.shell.show('home');
    flush();
    undo();
    flush();
    // Unregistered while paused: its body goes.
    expect(extra.probe.cleanups).toBe(1);
    const again = probeView('extra');
    m.shell.registerViews([again.def]);
    m.shell.show('extra');
    flush();
    expect(again.probe.mounts).toBe(1);
    m.dispose();
  });

  it('folding the rail’s group of a split column does not rebuild the tree (the main area stays put)', () => {
    const m = mountShell(views());
    m.shell.show('layers');
    m.shell.splitBelow('inspect');
    flush();
    const main = m.main();
    const parent = main.parentElement;
    m.shell.collapse('right');
    flush();
    m.shell.toggleRail(false);
    flush();
    m.shell.toggleRail(true);
    flush();
    expect(main.parentElement).toBe(parent);
    expect(main.isConnected).toBe(true);
    m.dispose();
  });

  it('a window none of whose views can show is hidden (no empty box without a header)', () => {
    const shown = new State(true);
    const m = mountShell([probeView('home'), probeView('tool', { home: 'float', visible: shown })]);
    m.shell.show('tool');
    flush();
    const w = m.q('[data-window]')!;
    expect(w.hidden).toBe(false);
    shown.set(false);
    flush();
    expect(w.hidden).toBe(true);
    m.dispose();
  });

  it('a rail item follows its own visible flag', () => {
    const m = mountShell(views());
    const visible = new State(true);
    m.shell.registerRailItems([{ id: 'pin', icon: 'pin', title: 'Pin', group: 'app', visible, run: () => {} }]);
    flush();
    expect(m.q('nav button[aria-label="Pin"]')).not.toBeNull();
    visible.set(false);
    flush();
    expect(m.q('nav button[aria-label="Pin"]')).toBeNull();
    m.dispose();
  });

  it('compact: with a window open, a docked view picked from the bar shows', () => {
    const m = mountShell(views(), { mode: 'compact' });
    m.shell.float('inspect');
    flush();
    expect(m.frame('inspect')?.getAttribute('data-placement')).toBe('sheet');
    const bar = m.q('[role="tablist"][aria-label="Views"]')!;
    [...bar.querySelectorAll<HTMLElement>('button')].find((b) => b.textContent?.includes('Layers'))!.click();
    flush();
    expect(m.frame('layers')?.getAttribute('data-placement')).toBe('sheet');
    m.dispose();
  });

  it('show() of a view under the rail’s group brings the hidden rail back', () => {
    const m = mountShell(views());
    m.shell.show('layers');
    m.shell.splitBelow('inspect');
    m.shell.toggleRail(false);
    flush();
    m.shell.show('inspect');
    flush();
    expect(m.shell.railShown.value).toBe(true);
    m.dispose();
  });

  it('a card covering the right column makes its groups inert and paused; a page covers the windows', () => {
    const m = mountShell([...views(), probeView('export', { home: 'card' }), probeView('report', { home: 'page' })]);
    m.shell.show('layers');
    flush();
    m.shell.show('export');
    flush();
    expect(m.groupNode('right')?.hasAttribute('inert')).toBe(true);
    expect(m.probes.layers.scope?.active).toBe(false);
    m.shell.close('export');
    flush();
    expect(m.groupNode('right')?.hasAttribute('inert')).toBe(false);
    m.shell.float('inspect');
    m.shell.show('report');
    flush();
    expect(m.q('[data-window]')!.closest('[inert]')).not.toBeNull();
    m.dispose();
  });

  it('closing a split keeps what the person hid hidden', () => {
    const m = mountShell(views());
    m.shell.show('layers');
    m.shell.splitBelow('inspect');
    const lower = findView(m.shell.layout.value, 'inspect')!.group.id;
    m.shell.move('settings', { kind: 'group', id: lower });
    m.shell.close('settings');
    m.shell.closeGroup(lower);
    expect(m.shell.layout.value.hidden).toEqual(['settings']);
    expect(railGroup(m.shell.layout.value)?.views).toContain('settings');
    m.dispose();
  });

  it('docked sizes keep the main area its least width (a long drag, a stored size)', () => {
    const m = mountShell(views());
    m.shell.setWidth(900);
    m.shell.show('home');
    m.shell.resize('right', 5000);
    const size = (rightColumn(m.shell.layout.value) as { size?: number }).size!;
    expect(size).toBe(900 - RailWidth - MinMainSize);
    // Narrower: it gives way, down to its own least size.
    m.shell.setWidth(500);
    expect((rightColumn(m.shell.layout.value) as { size?: number }).size).toBe(500 - RailWidth - MinMainSize);
    m.shell.setWidth(300);
    expect((rightColumn(m.shell.layout.value) as { size?: number }).size).toBe(MinGroupSize);
    m.dispose();
  });

  it('the busy overlay holds the focus and makes the rest inert', async () => {
    const busy = new State<{ label: string; cancel?: () => void } | undefined>(undefined);
    const m = mountShell(views(), {}, { busy });
    busy.set({ label: 'Loading', cancel: () => busy.set(undefined) });
    flush();
    await microtasks();
    const dialog = m.q('[role="alertdialog"]')!;
    expect(dialog.contains(document.activeElement)).toBe(true);
    expect(m.main().closest('[inert]')).not.toBeNull();
    const tab = key(document.activeElement!, 'Tab');
    expect(tab.defaultPrevented).toBe(true);
    m.dispose();
  });

  it('a stored layout with a stacking order out of sequence is read renumbered; one with no number is refused', () => {
    const rect = { x: 0, y: 0, width: 300, height: 200 };
    const stored = (za: number, zb: number) => ({
      version: 1,
      root: { kind: 'main', id: 'main' },
      windows: [
        { group: { kind: 'group', id: 'a', views: ['x'] }, rect, z: za },
        { group: { kind: 'group', id: 'b', views: ['y'] }, rect, z: zb },
      ],
      hidden: [],
      rail: true,
      seq: 1,
    });
    const l = validateLayout(stored(7, 3))!;
    expect(l.windows.map((w) => [w.group.id, w.z])).toEqual([
      ['b', 1],
      ['a', 2],
    ]);
    expect(validateLayout(stored(Number.NaN, 1))).toBeUndefined();
  });

  it('a section shows while its `when` holds, asked again as its `whenChanged` changes', async () => {
    const { SectionList } = await import('../section-list.js');
    const { domView } = await import('../../dom/index.js');
    const ready = new State(false);
    const m = mountShell([
      probeView('export', {
        render: domView(({ ctx }) => <SectionList shell={ctx.shell} slot="export" ctx={ctx} />),
      }),
    ]);
    m.shell.registerSections([
      {
        id: 'late',
        slot: 'export',
        title: 'Late',
        when: () => ready.value,
        whenChanged: ready,
        render: domView(() => <p>late</p>),
      },
    ]);
    m.shell.show('export');
    flush();
    expect(m.q('[data-section="late"]')).toBeNull();
    ready.set(true);
    flush();
    expect(m.q('[data-section="late"]')).not.toBeNull();
    m.dispose();
  });
});
