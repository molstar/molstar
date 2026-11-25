/** @vitest-environment jsdom */
/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { afterEach, describe, expect, it, vi } from 'vitest';
import { flush, render } from '../../../dom/index.js';
import { state } from '../../../state.js';
import { type Choice, Choices } from '../choices.js';

afterEach(() => {
  document.body.innerHTML = '';
});

const styles: Choice[] = [
  { id: 'default', label: 'Default' },
  { id: 'cartoon', label: 'Cartoon' },
  { id: 'spacefill', label: 'Spacefill', disabled: 'Too many atoms' },
  { id: 'surface', label: 'Surface' },
];

function mount(content: () => Node) {
  const container = document.createElement('div');
  document.body.append(container);
  render(content, container);
  return container;
}

const key = (el: Element, k: string) => el.dispatchEvent(new KeyboardEvent('keydown', { key: k, bubbles: true }));

describe('Choices with an active choice: a radio group', () => {
  const setup = (...initial: [string | undefined] | []) => {
    const active = state<string | undefined>(initial.length ? initial[0] : 'cartoon');
    const onPick = vi.fn((id: string) => active.set(id));
    const container = mount(() => <Choices label="Quick styles" choices={styles} active={active} onPick={onPick} />);
    const radios = [...container.querySelectorAll<HTMLButtonElement>('[role=radio]')];
    const checked = () => radios.map((r) => r.getAttribute('aria-checked'));
    const stops = () => radios.map((r) => r.tabIndex);
    return { active, onPick, container, radios, checked, stops };
  };

  it('is a radio group named by its label; the checked choice is the one tab stop', () => {
    const { container, radios, checked, stops } = setup();
    const group = container.querySelector('[role=radiogroup]')!;
    const label = document.getElementById(group.getAttribute('aria-labelledby')!);
    expect(label?.textContent).toBe('Quick styles');
    expect(radios).toHaveLength(4);
    expect(checked()).toEqual(['false', 'true', 'false', 'false']);
    expect(stops()).toEqual([-1, 0, -1, -1]);
  });

  it('a click picks; the marks follow the value, through the two choices it concerns', async () => {
    const { active, radios, checked, stops, onPick } = setup();
    radios[3].click();
    expect(onPick).toHaveBeenCalledWith('surface');
    expect(checked()).toEqual(['false', 'false', 'false', 'true']);
    expect(stops()).toEqual([-1, -1, -1, 0]);
    const touched = new Set<Node>();
    const observer = new MutationObserver((records) => {
      for (const r of records) touched.add(r.target);
    });
    for (const r of radios) observer.observe(r, { attributes: true, attributeFilter: ['aria-checked', 'tabindex'] });
    active.set('default');
    flush();
    await Promise.resolve();
    expect(touched).toEqual(new Set([radios[0], radios[3]]));
    observer.disconnect();
  });

  it('the arrows, Home and End move and pick, wrapping and skipping a disabled choice', () => {
    const { radios, active } = setup();
    key(radios[1], 'ArrowRight');
    expect(active.value).toBe('surface');
    expect(document.activeElement).toBe(radios[3]);
    key(radios[3], 'ArrowDown');
    expect(active.value).toBe('default');
    key(radios[0], 'ArrowLeft');
    expect(active.value).toBe('surface');
    key(radios[3], 'ArrowUp');
    expect(active.value).toBe('cartoon');
    key(radios[1], 'End');
    expect(active.value).toBe('surface');
    key(radios[3], 'Home');
    expect(active.value).toBe('default');
  });

  it('a disabled choice is marked, carries its reason, and a click does not pick it', () => {
    const { radios, onPick } = setup();
    expect(radios[2].getAttribute('aria-disabled')).toBe('true');
    expect(radios[2].dataset.tip).toBe('Too many atoms');
    radios[2].click();
    expect(onPick).not.toHaveBeenCalled();
  });

  it('with nothing checked, the first enabled choice is the tab stop until one is', () => {
    const { active, checked, stops } = setup(undefined);
    expect(checked()).toEqual(['false', 'false', 'false', 'false']);
    expect(stops()).toEqual([0, -1, -1, -1]);
    active.set('surface');
    flush();
    expect(stops()).toEqual([-1, -1, -1, 0]);
    active.set('spacefill');
    flush();
    expect(checked()).toEqual(['false', 'false', 'true', 'false']);
    expect(stops()).toEqual([0, -1, -1, -1]);
  });
});

describe('Choices without an active choice: a row of actions', () => {
  it('is a group of plain buttons; a click runs one', () => {
    const onPick = vi.fn();
    const container = mount(() => (
      <Choices
        aria-label="Wiggle"
        choices={[
          { id: 'dynamics', label: 'Dynamics' },
          { id: 'clear', label: 'Clear', title: 'Stop the wiggle' },
        ]}
        onPick={onPick}
      />
    ));
    const group = container.querySelector('[role=group]')!;
    expect(group.getAttribute('aria-label')).toBe('Wiggle');
    const buttons = [...group.querySelectorAll('button')];
    expect(buttons.map((b) => b.getAttribute('role'))).toEqual([null, null]);
    expect(buttons.map((b) => b.hasAttribute('aria-checked'))).toEqual([false, false]);
    expect(buttons[1].dataset.tip).toBe('Stop the wiggle');
    // A UI label is never cut: no overflow tooltip, no fade.
    expect(buttons[0].querySelector('[data-tip-overflow]')).toBeNull();
    buttons[0].click();
    expect(onPick).toHaveBeenCalledWith('dynamics');
  });

  it('icon-only: each choice named and tipped by its label; columns are set when asked', () => {
    const container = mount(() => (
      <Choices
        iconOnly
        columns={2}
        choices={[
          { id: 'show', label: 'Show', icon: 'eye' },
          { id: 'hide', label: 'Hide', icon: 'eye-off', title: 'Hide all' },
        ]}
        onPick={() => {}}
      />
    ));
    const [show, hide] = container.querySelectorAll('button');
    expect(show.getAttribute('aria-label')).toBe('Show');
    expect(show.dataset.tip).toBe('Show');
    expect(hide.dataset.tip).toBe('Hide all');
    expect(show.textContent).toBe('');
    const row = container.querySelector<HTMLElement>('[role=group]')!;
    expect(row.style.getPropertyValue('--msui-choices-columns')).toBe('2');
  });
});
