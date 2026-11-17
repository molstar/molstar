/** @vitest-environment jsdom */
/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { State } from '@molstar/core/reactive';
import { afterEach, describe, expect, it, vi } from 'vitest';
import { flush, Motion, onCleanup, render } from '../../../dom/index.js';
import { Expandable } from '../expandable.js';

afterEach(() => {
  document.body.innerHTML = '';
});

function mount(node: () => Node) {
  const container = document.createElement('div');
  document.body.append(container);
  const dispose = render(node, container);
  flush();
  const header = () => container.querySelector<HTMLButtonElement>('button')!;
  const body = () => container.querySelector<HTMLElement>('[role="region"]')!;
  return { container, dispose, header, body };
}

describe('Expandable', () => {
  it('is a disclosure: a button in a heading, expanded state, controlling a region named by it', () => {
    const m = mount(() => <Expandable title="Background">{() => <p>body</p>}</Expandable>);
    const heading = m.container.querySelector('[role="heading"]')!;
    expect(heading.getAttribute('aria-level')).toBe('3');
    expect(heading.contains(m.header())).toBe(true);
    expect(m.header().getAttribute('aria-expanded')).toBe('false');
    expect(m.header().getAttribute('aria-controls')).toBe(m.body().id);
    expect(m.body().getAttribute('aria-labelledby')).toBe(m.header().id);
    expect(m.body().hidden).toBe(true);
    m.dispose();
  });

  it('renders the body on first open (lazy), then keeps it, inert and paused while closed', () => {
    let renders = 0;
    const label = new State('one');
    const m = mount(() => (
      <Expandable title="Lighting">
        {() => {
          renders++;
          return <p>{label}</p>;
        }}
      </Expandable>
    ));
    expect(renders).toBe(0);
    m.header().click();
    flush();
    expect(renders).toBe(1);
    expect(m.body().hidden).toBe(false);
    expect(m.body().textContent).toBe('one');
    expect(m.container.firstElementChild!.hasAttribute('data-open')).toBe(true);
    m.header().click();
    flush();
    expect(m.body().hidden).toBe(true);
    expect(m.body().querySelector('p')!.hasAttribute('inert')).toBe(true);
    // Paused: closed, it writes nothing; opened again, the current value, without rendering anew.
    label.set('two');
    flush();
    expect(m.body().textContent).toBe('one');
    m.header().click();
    flush();
    expect(m.body().textContent).toBe('two');
    expect(renders).toBe(1);
    m.dispose();
  });

  it("a State is the host's: pressing sets it, setting it opens and closes; a plain value starts it", () => {
    const open = new State(false);
    const toggles: boolean[] = [];
    const m = mount(() => (
      <Expandable title="Camera" open={open} onToggle={(o) => toggles.push(o)}>
        {() => <p>body</p>}
      </Expandable>
    ));
    m.header().click();
    flush();
    expect(open.value).toBe(true);
    expect(toggles).toEqual([true]);
    expect(m.body().hidden).toBe(false);
    open.set(false);
    flush();
    expect(m.header().getAttribute('aria-expanded')).toBe('false');
    expect(m.body().hidden).toBe(true);
    m.dispose();
    const started = mount(() => (
      <Expandable title="Camera" open>
        {() => <p>body</p>}
      </Expandable>
    ));
    expect(started.body().hidden).toBe(false);
    // The first render stays still: no reveal until the user toggles.
    expect(started.body().hasAttribute('data-reveal')).toBe(false);
    started.dispose();
  });

  it('closing plays back first: the chevron turns at once, the body hides when the closing ends', async () => {
    const saved = Motion.close;
    let finish!: () => void;
    const anim = { finished: new Promise<void>((r) => (finish = r)), cancel: vi.fn() } as unknown as Animation;
    Motion.close = vi.fn(() => anim);
    try {
      const m = mount(() => (
        <Expandable title="Background" open>
          {() => <p>body</p>}
        </Expandable>
      ));
      m.header().click();
      flush();
      expect(Motion.close).toHaveBeenCalledWith(m.body(), 'grow');
      expect(m.header().getAttribute('aria-expanded')).toBe('false');
      expect(m.body().hidden).toBe(false);
      finish();
      await anim.finished;
      await Promise.resolve();
      flush();
      expect(m.body().hidden).toBe(true);
      // Reopened while closing: the closing stops, it stays open.
      m.header().click();
      flush();
      m.header().click();
      flush();
      m.header().click();
      flush();
      expect(anim.cancel).toHaveBeenCalled();
      expect(m.body().hidden).toBe(false);
      m.dispose();
    } finally {
      Motion.close = saved;
    }
  });

  it('eager renders at once; its cleanups run with the panel', () => {
    const cleaned = vi.fn();
    const m = mount(() => (
      <Expandable title="Background" eager>
        {() => {
          onCleanup(cleaned);
          return <p>body</p>;
        }}
      </Expandable>
    ));
    expect(m.body().querySelector('p')).not.toBeNull();
    m.dispose();
    expect(cleaned).toHaveBeenCalledTimes(1);
  });
});
