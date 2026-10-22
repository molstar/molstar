/** @vitest-environment jsdom */
/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { of } from 'rxjs';
import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { render } from '../../../dom/render.js';
import { flush } from '../../../dom/scope.js';
import { NotificationsManager } from '../../../managers/notifications.js';
import { StatusLine, StatusPill } from '../canvas-chrome.js';

// The status column's lines: a finished one goes on its own, held while hovered (a ring around its ✕ counts it down);
// an error stays until its ✕ (shell.md, the status column).

beforeEach(() => vi.useFakeTimers());
afterEach(() => {
  vi.useRealTimers();
  document.body.innerHTML = '';
});

function mount(child: () => Node) {
  const host = document.createElement('div');
  document.body.append(host);
  render(child, host);
  return host;
}

describe('the status column', () => {
  it('dismisses a line on its own after its time, held while the pointer is over it', () => {
    const dismissed = vi.fn();
    const host = mount(() => <StatusPill icon={null} text="load 1tqn" onDismiss={dismissed} autoDismiss={of(1000)} />);
    expect(host.querySelector('svg circle')).not.toBeNull();
    vi.advanceTimersByTime(600);
    const line = host.firstElementChild as HTMLElement;
    line.dispatchEvent(new Event('pointerenter'));
    vi.advanceTimersByTime(2000);
    expect(dismissed).not.toHaveBeenCalled();
    line.dispatchEvent(new Event('pointerleave'));
    vi.advanceTimersByTime(300);
    expect(dismissed).not.toHaveBeenCalled();
    vi.advanceTimersByTime(200);
    expect(dismissed).toHaveBeenCalledTimes(1);
  });

  it('keeps an error until its ✕', () => {
    const status = new NotificationsManager();
    const host = mount(() => <StatusLine status={status} />);
    status.post('state-tree could not load', 'error');
    vi.advanceTimersByTime(60_000);
    const line = () => host.querySelector('[data-kind="error"]');
    expect(line()?.textContent).toContain('state-tree could not load');
    (host.querySelector('button[aria-label="Dismiss"]') as HTMLButtonElement).click();
    expect(line()).toBeNull();
  });

  it('opens the log at the problems from an error, at everything from another message', () => {
    const status = new NotificationsManager();
    const opened = vi.fn();
    const host = mount(() => <StatusLine status={status} onOpenLog={opened} />);
    status.post('volume streaming failed', 'error');
    (host.querySelector('button[title="Show the problems in the log"]') as HTMLButtonElement).click();
    expect(opened).toHaveBeenCalledTimes(1);
    expect(status.logFilter.value).toBe('problems');
    status.post('loaded 1iep', 'info');
    (host.querySelector('button[title="Show the log"]') as HTMLButtonElement).click();
    expect(opened).toHaveBeenCalledTimes(2);
    expect(status.logFilter.value).toBe('all');
  });

  it('keeps each task pill while its task progresses: the spinner, the bar and Cancel stay (perf audit 11/B3)', () => {
    const status = new NotificationsManager();
    const host = mount(() => <StatusLine status={status} />);
    const aborted: number[] = [];
    // As the viewer mirrors the plugin's tasks: every task a new object on every progress event.
    const tasks = (i: number) => [
      { id: '1', label: `Loading ${i}`, progress: i / 10, abort: () => aborted.push(i) },
      { id: '2', label: 'Computing' },
    ];
    status.setTasks(tasks(0));
    flush();
    const spinners = () => [...host.querySelectorAll('svg')];
    const cancel = () => host.querySelectorAll<HTMLButtonElement>('button[title="Cancel"]:not([hidden])');
    const first = spinners();
    const button = cancel()[0];
    // Two spinners and the two Cancels' ✕ (the second's hidden).
    expect(first).toHaveLength(4);
    expect(cancel()).toHaveLength(1);
    let idsTold = 0;
    const unlisten = status.lineTaskIds.listen!(() => idsTold++);
    for (let i = 1; i <= 5; i++) {
      status.setTasks(tasks(i));
      flush();
    }
    unlisten();
    // Progress tells the pill, not the list.
    expect(idsTold).toBe(0);
    expect(spinners()).toEqual(first);
    expect(cancel()[0]).toBe(button);
    expect(host.textContent).toContain('Loading 5');
    const bars = host.querySelectorAll<HTMLElement>('[style*="width"]');
    expect(bars[0].style.width).toBe('50%');
    // The second task has no progress and no Cancel (its pill keeps its end padding: no button shown in it).
    // The label's span (the last span reading it: the pill reads it too) and its pill.
    const computing = [...host.querySelectorAll('span')].filter((e) => e.textContent === 'Computing').at(-1);
    expect(computing?.parentElement?.querySelector('button:not([hidden])')).toBeNull();
    // Cancel runs the task's latest abort.
    button.click();
    expect(aborted).toEqual([5]);
    status.setTasks([tasks(6)[1]]);
    flush();
    expect(cancel()).toHaveLength(0);
    expect(host.textContent).not.toContain('Loading');
  });
});
