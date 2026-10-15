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
import { afterEach, describe, expect, it, vi } from 'vitest';
import { headlessViewer, until } from '../../_test/harness.js';
import type { StorageEntry, StorageSource } from '../../storage/model.js';
import { settingsTab } from '../settings.js';

// The Storage view in the Settings tab (journal.md §7): folded, it reads nothing; open, each source's entries with their
// sizes and Forget, and Forget everything asked first.

afterEach(() => {
  document.body.innerHTML = '';
});

function mount(tab: ShellTab) {
  const host = document.createElement('div');
  document.body.append(host);
  tab.render.mount(host, { shell: new ShellModel() } as never, undefined as never);
  return host;
}

function fakeSource(entries: StorageEntry[]): StorageSource & { list: ReturnType<typeof vi.fn> } {
  const kept = [...entries];
  return {
    id: 'fake',
    title: 'Fake things',
    where: 'memory',
    list: vi.fn(async () => [...kept]),
    forget: async (id) => {
      kept.splice(
        kept.findIndex((e) => e.id === id),
        1,
      );
    },
  };
}

const buttonNamed = (host: HTMLElement, text: string) =>
  [...host.querySelectorAll('button')].find((b) => b.textContent?.trim() === text || b.title === text) as
    | HTMLButtonElement
    | undefined;

describe('the Storage view', () => {
  it('lists while open; Forget removes an entry; Forget everything asks first', async () => {
    const v = await headlessViewer();
    const source = fakeSource([
      { id: 'a', label: 'Thing A', bytes: 2048, detail: 'Kept a while ago' },
      { id: 'b', label: 'Thing B', stale: true },
    ]);
    v.storage.register(source);
    const host = mount(settingsTab(v.sources, v.managers, v.api.as({ kind: 'person' }), v.storage));
    const section = host.querySelector('section[aria-label="Storage"]') as HTMLElement;
    expect(section).toBeTruthy();
    // Folded: nothing is read.
    expect(source.list).not.toHaveBeenCalled();
    (section.querySelector('button[aria-expanded]') as HTMLButtonElement).click();
    await until(() => section.textContent?.includes('Thing A') === true);
    expect(section.textContent).toContain('Fake things');
    expect(section.textContent).toContain('2.0 KB');
    expect(section.textContent).toContain('Kept a while ago');
    expect(section.textContent).toContain('Kept records');
    expect(section.querySelector('[data-stale]')?.textContent).toContain('Thing B');

    buttonNamed(section, 'Forget Thing A')!.click();
    await until(() => !section.textContent?.includes('Thing A'));
    expect(section.textContent).toContain('Thing B');

    buttonNamed(section, 'Forget everything')!.click();
    expect(section.textContent).toContain('cannot be undone');
    // Cancel keeps everything.
    buttonNamed(section, 'Cancel')!.click();
    expect(section.textContent).toContain('Thing B');
    buttonNamed(section, 'Forget everything')!.click();
    const confirm = section.querySelector('[role=group] button') as HTMLButtonElement;
    confirm.click();
    await until(() => !section.textContent?.includes('Thing B'));
    await until(() => section.textContent?.includes('Nothing is kept.') === true);
  });

  it('is not shown without a storage model', async () => {
    const v = await headlessViewer();
    const host = mount(settingsTab(v.sources, v.managers, v.api.as({ kind: 'person' })));
    expect(host.querySelector('section[aria-label="Storage"]')).toBeNull();
  });
});
