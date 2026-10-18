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
