/** @vitest-environment jsdom */
/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { type Child, flush, render } from '../../dom/index.js';
import { StatusChips, StatusLine, StatusPill } from '../status-line.js';
import { StatusModel } from '../status-model.js';
import { axeViolations } from '../../testing/axe.js';

// The status column's lines (status §2): a finished one goes on its own, held while hovered (a ring around its ✕
// counts it down); an error stays until its ✕; a task's pill is made once and follows its progress.

beforeEach(() => vi.useFakeTimers());
afterEach(() => {
  vi.useRealTimers();
  document.body.innerHTML = '';
});

function mount(child: () => Child) {
  const host = document.createElement('div');
  document.body.append(host);
  render(child, host);
  return host;
}

describe('StatusPill', () => {
  it('dismisses a line on its own after its time, held while the pointer is over it', () => {
    const dismissed = vi.fn();
    const host = mount(() => <StatusPill icon={null} text="load 1tqn" onDismiss={dismissed} autoDismiss={1000} />);
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

  it('is held while the focus is in it', () => {
    const dismissed = vi.fn();
    const host = mount(() => (
      <StatusPill icon={null} text="done" onClick={() => {}} onDismiss={dismissed} autoDismiss={1000} />
    ));
    const button = host.querySelector('button') as HTMLButtonElement;
    button.dispatchEvent(new FocusEvent('focusin', { bubbles: true }));
    vi.advanceTimersByTime(5000);
    expect(dismissed).not.toHaveBeenCalled();
    button.dispatchEvent(new FocusEvent('focusout', { bubbles: true }));
    vi.advanceTimersByTime(1000);
    expect(dismissed).toHaveBeenCalledTimes(1);
  });
});

describe('StatusLine', () => {
  it('speaks each message in its polite region and shows it until it fades', () => {
    const status = new StatusModel();
    const host = mount(() => <StatusLine status={status} />);
    status.post('Loaded 1tqn', 'success');
    flush();
    const region = host.querySelector('[role=status]')!;
    expect(region.getAttribute('aria-live')).toBe('polite');
    expect(region.textContent).toBe('Loaded 1tqn');
    expect(host.querySelector('[data-kind="success"]')?.textContent).toContain('Loaded 1tqn');
    vi.advanceTimersByTime(4000);
    flush();
    expect(host.querySelector('[data-kind="success"]')).toBeNull();
  });

  it('keeps an error until its ✕', () => {
    const status = new StatusModel();
    const host = mount(() => <StatusLine status={status} />);
    status.post('state could not load', 'error');
    vi.advanceTimersByTime(60_000);
    flush();
    const line = () => host.querySelector('[data-kind="error"]');
    expect(line()?.textContent).toContain('state could not load');
    (host.querySelector('button[aria-label="Dismiss"]') as HTMLButtonElement).click();
    expect(line()).toBeNull();
  });

  it('opens the log at the problems from an error, at everything from another message', () => {
    const status = new StatusModel();
    const opened = vi.fn();
    const host = mount(() => <StatusLine status={status} onOpenLog={opened} />);
    status.post('volume streaming failed', 'error');
    flush();
    (host.querySelector('button[data-tip="Show the problems in the log"]') as HTMLButtonElement).click();
    expect(opened).toHaveBeenCalledWith('error');
    expect(status.log.filter.value).toEqual({ kinds: ['warning', 'error'] });
    status.post('loaded 1iep', 'info');
    flush();
    (host.querySelector('button[data-tip="Show the log"]') as HTMLButtonElement).click();
    expect(opened).toHaveBeenCalledTimes(2);
    expect(status.log.filter.value).toEqual({});
  });

  it('keeps each task pill while its task progresses: the spinner, the bar and Cancel stay', () => {
    const status = new StatusModel();
    const host = mount(() => <StatusLine status={status} />);
    const aborted: number[] = [];
    // As an app mirrors its tasks: every task a new object on every progress event.
    const tasks = (i: number) => [
      { id: '1', label: `Loading ${i}`, progress: i / 10, abort: () => aborted.push(i) },
      { id: '2', label: 'Computing' },
    ];
    status.setTasks(tasks(0));
    flush();
    const pills = () => [...host.querySelectorAll('[data-task]')];
    const cancel = () => host.querySelectorAll<HTMLButtonElement>('button[aria-label^="Cancel"]:not([hidden])');
    const first = pills();
    const button = cancel()[0];
    expect(first).toHaveLength(2);
    expect(cancel()).toHaveLength(1);
    for (let i = 1; i <= 5; i++) {
      status.setTasks(tasks(i));
      flush();
    }
    expect(pills()).toEqual(first);
    expect(cancel()[0]).toBe(button);
    // Named after its task (two running tasks are told apart), following its label.
    expect(button.getAttribute('aria-label')).toBe('Cancel Loading 5');
    expect(host.textContent).toContain('Loading 5');
    const bar = first[0].querySelector('[role=progressbar]') as HTMLElement;
    expect(bar.getAttribute('aria-valuenow')).toBe('50');
    expect((bar.firstElementChild as HTMLElement).style.width).toBe('50%');
    expect((first[1].querySelector('[role=progressbar]') as HTMLElement).style.display).toBe('none');
    button.click();
    expect(aborted).toEqual([5]);
    status.setTasks([tasks(6)[1]]);
    flush();
    expect(cancel()).toHaveLength(0);
    expect(host.textContent).not.toContain('Loading');
  });
});

describe('StatusLine labels', () => {
  it('takes its text from `labels` (the rest English): the log tips, the ✕, each task’s Cancel, a chip’s ✕', () => {
    const status = new StatusModel();
    const host = mount(() => (
      <div>
        <StatusLine
          status={status}
          onOpenLog={() => {}}
          labels={{ showProblems: 'Probleme', dismiss: 'Schließen', cancel: (t) => `${t} abbrechen` }}
        />
        <StatusChips status={status} labels={{ clear: (c) => `${c} leeren`, clearTip: 'Leeren' }} />
      </div>
    ));
    status.post('Kaputt', 'error');
    status.setTasks([
      { id: '1', label: 'Laden', abort: () => {} },
      { id: '2', label: 'Rechnen', abort: () => {} },
    ]);
    status.setChips([{ id: 'sel', label: '3 Reste', dismiss: () => {} }]);
    flush();
    expect(host.querySelector('button[data-tip="Probleme"]')).not.toBeNull();
    expect(host.querySelector('button[aria-label="Schließen"]')).not.toBeNull();
    const cancels = [...host.querySelectorAll('[data-task] button')].map((b) => b.getAttribute('aria-label'));
    expect(cancels).toEqual(['Laden abbrechen', 'Rechnen abbrechen']);
    const clear = host.querySelector('button[aria-label="3 Reste leeren"]')!;
    expect(clear.getAttribute('data-tip')).toBe('Leeren');
    status.retire();
  });

  it('a pill on its own takes its ✕’s name from `labels`', () => {
    const host = mount(() => (
      <StatusPill icon={null} text="Fertig" onDismiss={() => {}} labels={{ dismiss: 'Ausblenden' }} />
    ));
    const x = host.querySelector('button')!;
    expect([x.getAttribute('aria-label'), x.getAttribute('data-tip')]).toEqual(['Ausblenden', 'Ausblenden']);
  });

  it('names each task’s Cancel after its task by default', () => {
    const status = new StatusModel();
    const host = mount(() => <StatusLine status={status} onOpenLog={() => {}} />);
    status.post('Broken', 'error');
    status.setTasks([{ id: '1', label: 'Loading', abort: () => {} }]);
    flush();
    expect(host.querySelector('button[data-tip="Show the problems in the log"]')).not.toBeNull();
    expect(host.querySelector('[data-task] button')?.getAttribute('aria-label')).toBe('Cancel Loading');
    status.retire();
  });
});

describe('StatusLine scopes (a line a pane, say)', () => {
  it('a scope’s line shows its messages; the line taking the rest shows the others, the tasks, and speaks', () => {
    const status = new StatusModel();
    const host = mount(() => (
      <div>
        <div data-line="b">
          <StatusLine status={status} scope={{ id: 'b' }} />
        </div>
        <div data-line="rest">
          <StatusLine status={status} scope={{ rest: true }} />
        </div>
      </div>
    ));
    const line = (id: string) => host.querySelector(`[data-line="${id}"]`)!;
    status.post('Camera reset: Stage B', 'info', { scope: 'b' });
    flush();
    expect(line('b').textContent).toContain('Camera reset: Stage B');
    expect(line('rest').querySelector('[data-kind]')).toBeNull();
    status.post('Saved the session');
    flush();
    expect(line('b').querySelector('[data-kind]')).toBeNull();
    expect(line('rest').querySelector('[data-kind]')?.textContent).toContain('Saved the session');
    // One live region in the app (the rest's); the tasks there too.
    expect(host.querySelectorAll('[role=status]')).toHaveLength(1);
    expect(line('rest').querySelector('[role=status]')).not.toBeNull();
    status.setTasks([{ id: '1', label: 'Loading' }]);
    flush();
    expect(line('rest').querySelector('[data-task]')).not.toBeNull();
    expect(line('b').querySelector('[data-task]')).toBeNull();
  });
});

describe('StatusChips', () => {
  it('activates a chip and clears it by its ✕', () => {
    const status = new StatusModel();
    const activate = vi.fn();
    const host = mount(() => <StatusChips status={status} />);
    status.setChips([{ id: 'sel', label: '3 residues', activate, dismiss: () => status.setChips([]) }]);
    flush();
    (host.querySelector('button') as HTMLButtonElement).click();
    expect(activate).toHaveBeenCalledTimes(1);
    (host.querySelector('button[aria-label="Clear 3 residues"]') as HTMLButtonElement).click();
    flush();
    expect(host.textContent).not.toContain('3 residues');
  });

  it('is accessible (axe): a message, an error, a task and a chip', async () => {
    vi.useRealTimers();
    const status = new StatusModel();
    const host = mount(() => (
      <div>
        <StatusLine status={status} onOpenLog={() => {}} />
        <StatusChips status={status} />
      </div>
    ));
    status.post('Loaded 1tqn', 'success');
    status.post('No such file', 'error');
    status.setTasks([{ id: '1', label: 'Loading', progress: 0.4, abort: () => {} }]);
    status.setChips([{ id: 'sel', label: '3 residues', activate: () => {}, dismiss: () => {} }]);
    flush();
    expect(await axeViolations(host)).toBe('');
  });
});
