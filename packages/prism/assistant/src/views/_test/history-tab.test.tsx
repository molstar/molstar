/** @vitest-environment jsdom */
/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { NotificationsManager, observableState } from '@molstar/prism-ui/core';
import { render } from '@molstar/prism-ui/dom/render';
import { afterEach, describe, expect, it, vi } from 'vitest';
import type { AgentModel, CommandRun, HistoryItem } from '../../model.js';
import { HistoryTab } from '../history-tab.js';

afterEach(() => {
  vi.useRealTimers();
  document.body.innerHTML = '';
});

const run = (id: string, time: Date): HistoryItem => ({
  kind: 'run',
  run: { id, command: 'color', title: `Color ${id}`, source: { kind: 'person' }, status: 'done', time } as CommandRun,
});

describe('the History tab', () => {
  it('says when an item ran as a time of day, which does not go stale (perf audit 11/B4)', () => {
    vi.useFakeTimers();
    const start = new Date(2026, 9, 22, 9, 5);
    vi.setSystemTime(start);
    const history = observableState<readonly HistoryItem[]>([run('a', start)]);
    const agent = { history: history.asReadonly } as unknown as AgentModel;
    const host = document.createElement('div');
    document.body.append(host);
    render(() => <HistoryTab agent={agent} status={new NotificationsManager()} />, host);
    const when = () => [...host.querySelectorAll('time')].map((t) => t.textContent);
    const first = when()[0];
    expect(first).not.toMatch(/now|ago/);
    expect(host.querySelector('time')?.getAttribute('datetime')).toBe(start.toISOString());
    // Twenty minutes on, a newer item: the first still says its time, not "just now".
    vi.setSystemTime(new Date(2026, 9, 22, 9, 25));
    history.set([run('a', start), run('b', new Date())]);
    expect(when()).toHaveLength(2);
    expect(when()).toContain(first);
    expect(new Set(when()).size).toBe(2);
  });
});
