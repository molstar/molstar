/** @vitest-environment jsdom */
/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { afterEach, describe, expect, it } from 'vitest';
import { render } from '../../../dom/render.js';
import { NotificationsManager } from '../../../managers/notifications.js';
import { LogView } from '../log-view.js';

afterEach(() => {
  document.body.innerHTML = '';
});

describe('the Log tab', () => {
  it('lists the newest message first', () => {
    const status = new NotificationsManager();
    const host = document.createElement('div');
    document.body.append(host);
    render(() => <LogView status={status} />, host);
    status.post('first', 'info');
    status.post('second', 'error');
    const texts = [...host.querySelectorAll('[data-kind]')].map((e) => e.textContent);
    expect(texts[0]).toContain('second');
    expect(texts[1]).toContain('first');
  });
});

describe('the Log tab filters', () => {
  const texts = (host: HTMLElement) => [...host.querySelectorAll('[data-kind]')].map((e) => e.textContent ?? '');
  const pick = (host: HTMLElement, label: string) =>
    ([...host.querySelectorAll('button')].find((b) => b.textContent?.trim() === label) as HTMLButtonElement).click();

  it('shows the details, a filter per source written, the problems, and filters by text', async () => {
    const status = new NotificationsManager();
    const host = document.createElement('div');
    document.body.append(host);
    render(() => <LogView status={status} />, host);
    status.post('Colored the polymer', 'info', { detail: '12 ms' });
    status.post('No such chain', 'warning');
    status.log('Created Model', 'info', 'Transforms', { detail: '3 ms' });
    status.log('Created Structure', 'info', 'Transforms', { detail: '5 ms' });
    await Promise.resolve();
    await Promise.resolve();
    expect(texts(host)).toHaveLength(4);
    expect(texts(host)[0]).toContain('5 ms');
    expect(texts(host).find((t) => t.includes('Colored'))).toContain('12 ms');

    pick(host, 'Transforms');
    await Promise.resolve();
    expect(texts(host)).toHaveLength(2);
    expect(texts(host)[0]).toMatch(/^Created Structure5 ms/);
    expect(texts(host)[1]).toMatch(/^Created Model3 ms/);
    pick(host, 'Messages');
    await Promise.resolve();
    expect(texts(host)).toHaveLength(2);
    pick(host, 'Problems');
    await Promise.resolve();
    expect(texts(host)).toHaveLength(1);
    expect(texts(host)[0]).toContain('No such chain');

    pick(host, 'All');
    const input = host.querySelector('input') as HTMLInputElement;
    input.value = 'model';
    input.dispatchEvent(new Event('input'));
    await Promise.resolve();
    expect(texts(host)).toHaveLength(1);
    expect(texts(host)[0]).toContain('Created Model');
    input.value = 'nothing like it';
    input.dispatchEvent(new Event('input'));
    await Promise.resolve();
    expect(host.textContent).toContain('Nothing matches.');
  });
});
