/** @vitest-environment jsdom */
/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ShellModel } from '@molstar/prism-ui/shell/shell';
import { DefaultPrismViewerSpec } from '@molstar/prism-viewer/default-spec';
import { createPrismViewerHeadless } from '@molstar/prism-viewer/headless';
import type { PrismViewer } from '@molstar/prism-viewer/viewer';
import { afterEach, describe, expect, it, vi } from 'vitest';
import type { KinemageInfo } from '../operations.js';
import { KinemageTab, kinemageUI } from '../prism-panel.js';
import {
  animateMenu,
  framesDescription,
  groupDescription,
  kinemagePanelLayout,
  kinemageDescription,
  kinemageShown,
  mastersDescription,
  sameLayout,
  sameShown,
  shownKey,
} from '../prism-panel-model.js';
import { KinemageExtension } from '../prism-viewer.js';
import { Kin } from './kin-fixture.js';

// The Kinemage panel (prism-panel.tsx): its model from a hand-written listing, then the panel in jsdom on a headless viewer
// with the test kinemage: its rows' eyes and its header's buttons call the `kinemage` operations as the person, and follow the state.

const info: KinemageInfo = {
  index: 1,
  title: 'Test',
  ref: 'k1',
  views: [{ key: '1', name: 'Overview' }],
  groups: [
    { name: 'Contacts', visible: true, nobutton: false, subgroups: [] },
    {
      name: 'Structure',
      visible: false,
      nobutton: false,
      animate: 1,
      subgroups: [
        { name: 'Backbone', visible: true, nobutton: false },
        { name: 'Hidden', visible: true, nobutton: true },
      ],
    },
    { name: 'Alternate', visible: true, nobutton: false, animate: 1, subgroups: [] },
    { name: 'Internal', visible: true, nobutton: true, subgroups: [] },
  ],
  masters: [{ name: 'clashes', visible: false }],
  animate: ['Structure', 'Alternate'],
  animate2: [],
};

describe('the Kinemage panel model', () => {
  it('lays out the groups with their subgroups (as Group:Subgroup) and the masters, without the nobutton ones', () => {
    const [k] = kinemagePanelLayout([info]);
    expect(k.groups.map((g) => g.key)).toEqual(['Contacts', 'Structure', 'Alternate']);
    expect(k.groups[1]).toMatchObject({
      animate: 1,
      subgroups: [{ kind: 'subgroup', key: 'Structure:Backbone', label: 'Backbone' }],
    });
    expect(k.masters).toEqual([{ kind: 'master', key: 'clashes', label: 'clashes' }]);
    expect([k.animate, k.animate2]).toEqual([2, 0]);
  });

  it('keeps the layout when only what is shown changes, and tells the shown apart', () => {
    const hidden: KinemageInfo = { ...info, masters: [{ name: 'clashes', visible: true }] };
    expect(sameLayout(kinemagePanelLayout([info]), kinemagePanelLayout([hidden]))).toBe(true);
    expect(sameLayout(kinemagePanelLayout([info]), kinemagePanelLayout([info, { ...info, ref: 'k2' }]))).toBe(false);
    const a = kinemageShown([info]);
    expect(a.get(shownKey('k1', { kind: 'master', key: 'clashes' }))).toBe(false);
    expect(a.get(shownKey('k1', { kind: 'subgroup', key: 'Structure:Backbone' }))).toBe(true);
    expect(sameShown(a, kinemageShown([info]))).toBe(true);
    expect(sameShown(a, kinemageShown([hidden]))).toBe(false);
  });

  it("describes the rows: a kinemage's contents, the masters on, a group's subgroups", () => {
    const [k] = kinemagePanelLayout([info]);
    expect(kinemageDescription(k)).toBe('3 groups · 1 master · 1 view');
    expect(mastersDescription(info)).toBe('0 of 1 on');
    expect(groupDescription(k.groups[1])).toBe('1 subgroup');
    expect(groupDescription(k.groups[0])).toBe('');
  });

  it('lists the animate sequences for Next frame, named with their kinemage when there are several', () => {
    const one = kinemagePanelLayout([info]);
    expect(animateMenu(one)).toEqual([{ id: 'k1:animate', label: 'Frames', sequence: 'animate', kinemage: 'k1' }]);
    const two = kinemagePanelLayout([info, { ...info, ref: 'k2', index: 2, title: 'Other', animate2: ['Contacts'] }]);
    expect(animateMenu(two).map((f) => f.label)).toEqual([
      'Test › Frames',
      'Other › Frames',
      'Other › Frames (2animate)',
    ]);
  });

  it('says which animate frame is shown, or none', () => {
    expect(framesDescription(info, 'animate')).toBe('Alternate · 2 of 2');
    const none = { ...info, groups: info.groups.map((g) => ({ ...g, visible: false })) };
    expect(framesDescription(none, 'animate')).toBe('none of 2 shown');
  });
});

const viewers: PrismViewer[] = [];
afterEach(async () => {
  for (const v of viewers.splice(0)) await v.dispose();
  document.body.replaceChildren();
});

async function viewer() {
  const base = DefaultPrismViewerSpec();
  const v = await createPrismViewerHeadless({
    spec: { ...base, extensions: [...(base.extensions ?? []), KinemageExtension] },
  });
  viewers.push(v);
  return v;
}

/** A state's value, read through its stream (the toolbar's `visible`). */
function current<T>(o: { subscribe(fn: (v: T) => void): { unsubscribe(): void } }): T | undefined {
  let value: T | undefined;
  o.subscribe((v) => {
    value = v;
  }).unsubscribe();
  return value;
}

async function until(check: () => boolean, ms = 5000): Promise<void> {
  const t0 = Date.now();
  while (!check()) {
    if (Date.now() - t0 > ms) throw new Error(`timed out after ${ms} ms`);
    await new Promise((r) => setTimeout(r, 5));
  }
}

const load = async (v: PrismViewer) => {
  const r = await v.api.load.files({ files: [new File([Kin], 'test.kin')] });
  if (!r.ok) throw new Error(r.message);
};

describe('the Kinemage panel', () => {
  it('is a tab beside Structure while a kinemage is loaded, opened by itself when one loads', async () => {
    const v = await viewer();
    // The view area stood in for: a headless viewer has no rail to open the tab on.
    const person = v.api.as({ kind: 'person' });
    const open = vi.fn(async (_: { tab: string }) => ({ ok: true }));
    const api = new Proxy(person, { get: (t, k) => (k === 'view' ? { open } : Reflect.get(t, k)) });
    const ui = kinemageUI({ api, sources: v.sources, managers: v.managers });
    const [tab] = ui.tabs!;
    expect(tab).toMatchObject({ id: KinemageTab, icon: 'kinemage', order: 1.1 });
    expect(current(tab.visible!)).toBe(false);
    await load(v);
    expect(current(tab.visible!)).toBe(true);
    await until(() => open.mock.calls.length === 1);
    expect(open).toHaveBeenCalledWith({ tab: KinemageTab });
  });

  it('draws a tree of the views, frames, masters and groups; a view row applies it, an eye is a step that Undo takes back', async () => {
    const v = await viewer();
    await load(v);
    const person = v.api.as({ kind: 'person' });
    // `kinemage.view` watched (no canvas here: the operation itself would say so); the rest is the area as it is.
    const area = person.kinemage!;
    const applied = vi.fn(async (_: { view: string; kinemage: string }) => ({ ok: true }));
    const kinemage = new Proxy(area, {
      get: (t, k) => (k === 'view' ? applied : (Reflect.get(t, k) as unknown)),
    });
    const api = new Proxy(person, { get: (t, k) => (k === 'kinemage' ? kinemage : Reflect.get(t, k)) });
    const ui = kinemageUI({ api, sources: v.sources, managers: v.managers });
    const host = document.createElement('div');
    document.body.append(host);
    const unmount = ui.tabs![0].render.mount(host, { shell: new ShellModel() } as never, undefined as never);

    const label = (item: Element) => item.querySelector('[data-tree-row] [data-tip-part]')?.textContent;
    const item = (name: string) => [...host.querySelectorAll('[role="treeitem"]')].find((i) => label(i) === name)!;
    const row = (name: string) => item(name).querySelector<HTMLElement>('[data-tree-row]')!;
    const eye = (name: string) =>
      row(name).querySelector<HTMLButtonElement>('button[data-tip="Show"], button[data-tip="Hide"]')!;
    // Shown or hidden as its eye offers the other.
    const shown = (name: string) => eye(name)?.getAttribute('data-tip') === 'Hide';
    const rows = () =>
      [...host.querySelectorAll('[role="treeitem"]')].map((i) => {
        const name = label(i)!;
        return eye(name) ? [name, shown(name)] : [name];
      });
    expect(rows()).toEqual([
      ['Test kinemage for the Prism extension'],
      ['Views'],
      ['Overview'],
      // The animate groups under their row, which says the frame shown.
      ['Frames'],
      ['Structure', true],
      ['Alternate', false],
      ['Masters'],
      ['bad overlap', true],
      ['clashes', true],
      ['Contacts', true],
    ]);
    expect(row('Frames').textContent).toContain('Structure · 1 of 2');
    expect(row('Structure').textContent).toContain('2 subgroups');

    // A view's row applies it.
    row('Overview').click();
    expect(applied).toHaveBeenCalledWith({ view: '1', kinemage: expect.any(String) });

    const steps = v.managers.history.steps.value.length;
    eye('Contacts').click();
    await until(() => !shown('Contacts') && v.managers.history.steps.value.length === steps + 1);
    const listing = await v.api.kinemage!.get({});
    expect(listing.ok && (listing.value as { kinemages: KinemageInfo[] }).kinemages[0].groups[0].visible).toBe(false);

    // A group's row opens its subgroups.
    row('Structure').click();
    await until(() => !!item('Sidechains'));
    eye('Sidechains').click();
    // Its step recorded: the operation has ended (the eye follows the state before the shapes are built again).
    await until(() => !shown('Sidechains') && v.managers.history.steps.value.length === steps + 2);

    // Undo is the state's: the eyes follow it.
    await v.api.history.undo();
    await until(() => shown('Sidechains'));
    await v.api.history.undo();
    await until(() => shown('Contacts'));

    // Next frame, the header's button: one sequence, so it steps it.
    expect(
      [...host.querySelectorAll('button')].filter((b) => !b.closest('[data-tree-row]')).map((b) => b.textContent),
    ).toEqual(['Next frame']);
    [...host.querySelectorAll('button')].find((b) => b.textContent === 'Next frame')!.click();
    await until(() => !shown('Structure') && shown('Alternate'));
    unmount();
  });

  it('says so with no kinemage, and draws each of several', async () => {
    const v = await viewer();
    const ui = kinemageUI({ api: v.api.as({ kind: 'person' }), sources: v.sources, managers: v.managers });
    const host = document.createElement('div');
    document.body.append(host);
    const unmount = ui.tabs![0].render.mount(host, { shell: new ShellModel() } as never, undefined as never);
    expect(host.textContent).toContain('No kinemage is loaded.');
    await load(v);
    await load(v);
    const roots = () => [...host.querySelectorAll('[role="tree"] > [role="treeitem"]')];
    await until(() => roots().length === 2);
    unmount();
  });
});
