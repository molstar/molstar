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
import { StatusModel } from '../../../shell/status.js';
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
    const status = new StatusModel();
    const host = mount(() => <StatusLine status={status} />);
    status.post('state-tree could not load', 'error');
    vi.advanceTimersByTime(60_000);
    const line = () => host.querySelector('[data-kind="error"]');
    expect(line()?.textContent).toContain('state-tree could not load');
    (host.querySelector('button[aria-label="Dismiss"]') as HTMLButtonElement).click();
    expect(line()).toBeNull();
  });
});
