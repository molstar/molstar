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
import { cbs, expectOk, headlessViewer, until } from './harness.js';

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
  it('has the log as its second half; compact, Steps and Log are tabs', async () => {
    const v = await headlessViewer();
    v.managers.notifications.post('Said once', 'info');
    const host = mount(historyTab(v, v.api.as({ kind: 'person' }), v.managers.notifications));
    expect([...host.querySelectorAll('h3')].map((h) => h.textContent)).toEqual(['Steps', 'Log']);
    expect(host.textContent).toContain('Said once');

    const compact = document.createElement('div');
    document.body.append(compact);
    const shell = new ShellModel({ mode: 'compact' } as never);
    historyTab(v, v.api.as({ kind: 'person' }), v.managers.notifications).render.mount(
      compact,
      { shell, focus: 'log' } as never,
      undefined as never,
    );
    expect(compact.querySelector('h3')).toBeNull();
    expect(compact.textContent).toContain('Said once');
    (compact.querySelector('button[aria-pressed="false"]') as HTMLButtonElement).click();
    expect(compact.textContent).toContain('Nothing has happened yet.');
  });

  it('lists the runs of the record, and Undo undoes the latest', async () => {
    const v = await headlessViewer();
    const person = v.api.as({ kind: 'person' });
    const host = mount(historyTab(v, person, v.managers.notifications));
    expect(host.textContent).toContain('Nothing has happened yet.');
    const undo = () => host.querySelector('button') as HTMLButtonElement;
    expect(undo().disabled).toBe(true);
    expectOk(await v.api.settings.set({ group: 'viewer', values: { checkpointing: false } }));
    await until(() => host.querySelectorAll('[data-origin]').length === 1);
    expect(host.textContent).toContain('Changed the Viewer settings');
    expect(undo().disabled).toBe(false);
    undo().click();
    // The record keeps the step, marked undone.
    await until(() => host.querySelectorAll('[data-undone]').length === 1);
    expect(host.querySelector('[data-undone]')?.textContent).toContain('Undone');
    expect(undo().disabled).toBe(true);
  });

  it("shows an agent's turn as a group with its run, and Undo turn for the latest turn", async () => {
    const v = await headlessViewer();
    const host = mount(historyTab(v, v.api.as({ kind: 'person' }), v.managers.notifications));
    const turn = v.record.beginTurn('t1', 'Turn the checkpoints off');
    expectOk(
      await v.api.as({ kind: 'agent', turn: 't1' }).settings.set({ group: 'viewer', values: { checkpointing: false } }),
    );
    await turn.settle({ reply: 'Done.', outcome: 'done' });
    await until(() => host.textContent?.includes('Undo turn') === true);
    const row = host.querySelector('[data-origin]') as HTMLElement;
    expect(row.dataset.origin).toBe('agent');
    expect(row.textContent).toContain('Turn the checkpoints off');
    expect(row.textContent).toContain('Done.');
    expect(row.textContent).toContain('1 change');
  });

  it("shows a run's targets and changes, a turn's model and streaming reply, and who undid a turn", async () => {
    const v = await headlessViewer({ files: [cbs()] });
    const host = mount(historyTab(v, v.api.as({ kind: 'person' }), v.managers.notifications));
    expectOk(await v.api.as({ kind: 'person' }).select.set({ target: 'ligand' }));
    await until(() => host.textContent?.includes('ligand · selection') === true);
    const turn = v.record.beginTurn('t1', 'Make it a surface', undefined, 'anthropic:claude-sonnet');
    const end = v.managers.history.beginGroup('Assistant', { kind: 'agent', turn: 't1' });
    expectOk(await v.api.as({ kind: 'agent', turn: 't1' }).styles.apply({ style: 'spacefill' }));
    end();
    turn.reply('Working on it');
    await until(() => host.textContent?.includes('Under way') === true);
    expect(host.textContent).toContain('claude-sonnet');
    expect(host.querySelector('[title="anthropic:claude-sonnet"]')).not.toBeNull();
    // The reply streams (coalesced to ~300 ms), then the final text replaces it.
    await until(() => host.textContent?.includes('Working on it') === true, 3000);
    await turn.settle({ reply: 'Done.', outcome: 'done' });
    await until(() => host.textContent?.includes('Done.') === true);
    expect(host.textContent).not.toContain('Working on it');
    expectOk(await v.api.as({ kind: 'agent', turn: 't2' }).history.undoTurn({ turn: 't1' }));
    await until(() => host.textContent?.includes('Undone by the assistant') === true);
  });

  it('keeps the record on a button, and says so', async () => {
    const v = await headlessViewer();
    const host = mount(historyTab(v, v.api.as({ kind: 'person' }), v.managers.notifications));
    // An icon in the bar: its name says what it does, and that nothing is kept yet.
    expect(host.querySelector('button[aria-label^="Keep this session"]')).not.toBeNull();
    expect(host.querySelector('button[aria-label^="Kept in this browser"]')).toBeNull();
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
