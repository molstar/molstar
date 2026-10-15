/** @vitest-environment jsdom */
/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { ShellTab } from '@molstar/prism-ui';
import { ShellModel } from '@molstar/prism-ui/shell/shell';
import { afterEach, describe, expect, it } from 'vitest';
import { historyTab } from '../ui/history.js';
import { settingsTab } from '../ui/settings.js';
import { expectOk, headlessViewer, until } from './harness.js';

// The History and Settings tabs over a headless viewer's managers, sources and person API (the page's shell mounts the
// same tabs): what they show, and that their buttons and forms act through the API.

afterEach(() => {
  document.body.innerHTML = '';
});

function mount(tab: ShellTab) {
  const host = document.createElement('div');
  document.body.append(host);
  const shell = new ShellModel();
  tab.render.mount(host, { shell } as never, undefined as never);
  return host;
}

describe('the History tab', () => {
  it('lists the steps, and Undo undoes the latest', async () => {
    const v = await headlessViewer();
    const person = v.api.as({ kind: 'person' });
    const host = mount(historyTab(v.managers, person));
    expect(host.textContent).toContain('Nothing to undo yet.');
    const undo = () => host.querySelector('button') as HTMLButtonElement;
    expect(undo().disabled).toBe(true);
    expectOk(await v.api.settings.set({ group: 'viewer', values: { checkpointing: false } }));
    await until(() => host.querySelectorAll('[data-origin]').length === 1);
    expect(host.textContent).toContain('Changed the Viewer settings');
    expect(undo().disabled).toBe(false);
    undo().click();
    await until(() => host.querySelectorAll('[data-origin]').length === 0);
    expect(host.textContent).toContain('Nothing to undo yet.');
  });

  it('shows who made a step, and the reason it cannot be undone', async () => {
    const v = await headlessViewer();
    const host = mount(historyTab(v.managers, v.api.as({ kind: 'person' })));
    expectOk(
      await v.api.as({ kind: 'agent', turn: 't1' }).settings.set({ group: 'viewer', values: { checkpointing: false } }),
    );
    v.managers.history.invalidate('Changed outside the viewer.');
    await until(() => host.querySelector('[data-unsafe]') !== null);
    const row = host.querySelector('[data-origin]') as HTMLElement;
    expect(row.dataset.origin).toBe('agent');
    expect(row.textContent).toContain('Changed outside the viewer.');
    expect(row.textContent).toContain('Undo turn');
  });
});

describe('the Settings tab', () => {
  it('shows the groups the host has, and an edit goes through settings.set', async () => {
    const v = await headlessViewer();
    const host = mount(settingsTab(v.sources, v.managers, v.api.as({ kind: 'person' })));
    const titles = [...host.querySelectorAll('section[aria-label]')].map((s) => s.getAttribute('aria-label'));
    expect(titles).toContain('Viewer');
    expect(titles).toContain('Behaviors');
    expect(titles).not.toContain('Appearance');
    // The first group (Viewer) is open: its switch is Undo for each step.
    const box = host.querySelector(
      'section[aria-label="Viewer"] input[type=checkbox], section[aria-label="Viewer"] [role=switch]',
    ) as HTMLElement;
    expect(box).toBeTruthy();
    box.click();
    await until(() => v.managers.history.checkpointing.value === false);
  });
});
