/** @vitest-environment jsdom */
/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { BehaviorSubject, State } from '@molstar/core/reactive';
import { afterEach, describe, expect, it, vi } from 'vitest';
import { flush, render } from '../../../dom/index.js';
import { state } from '../../../state.js';
import { Segmented, type SegmentedOption } from '../segmented.js';

afterEach(() => {
  document.body.innerHTML = '';
});

const types: SegmentedOption<string>[] = [
  { value: 'cartoon', label: 'Cartoon', icon: 'cartoon' },
  { value: 'ball-and-stick', label: 'Ball & stick', icon: 'ball-and-stick' },
  { value: 'spacefill', label: 'Spacefill', icon: 'spacefill', disabled: 'Too many atoms' },
  { value: 'surface', label: 'Surface', icon: 'surface' },
];

function mount(content: () => Node) {
  const container = document.createElement('div');
  document.body.append(container);
  render(content, container);
  return container;
}

const key = (el: Element, k: string) => el.dispatchEvent(new KeyboardEvent('keydown', { key: k, bubbles: true }));

function setup(value: State<string> = state('ball-and-stick'), onChange = vi.fn()) {
  const container = mount(() => <Segmented label="Type" options={types} value={value} onChange={onChange} iconOnly />);
  const radios = [...container.querySelectorAll<HTMLButtonElement>('[role=radio]')];
  const checked = () => radios.map((r) => r.getAttribute('aria-checked'));
  const stops = () => radios.map((r) => r.tabIndex);
  const heading = () =>
    container.querySelector('[role=radiogroup]')!.previousElementSibling!.lastElementChild!.textContent;
  return { value, onChange, container, radios, checked, stops, heading };
}

describe('Segmented', () => {
  it('is a radio group named by its heading label; icon-only segments are named by their labels', () => {
    const { container, radios, checked, stops } = setup();
    const group = container.querySelector('[role=radiogroup]')!;
    expect(document.getElementById(group.getAttribute('aria-labelledby')!)?.textContent).toBe('Type');
    expect(radios.map((r) => r.getAttribute('aria-label'))).toEqual([
      'Cartoon',
      'Ball & stick',
      'Spacefill',
      'Surface',
    ]);
    expect(radios[0].dataset.tip).toBe('Cartoon');
    expect(radios[0].querySelector('svg')).not.toBeNull();
    expect(radios[0].textContent).toBe('');
    expect(checked()).toEqual(['false', 'true', 'false', 'false']);
    expect(stops()).toEqual([-1, 0, -1, -1]);
  });

  it('a click sets the state and tells onChange; the marks and the heading follow', () => {
    const { value, onChange, radios, checked, stops, heading } = setup();
    expect(heading()).toBe('Ball & stick');
    radios[3].click();
    expect(value.value).toBe('surface');
    expect(onChange).toHaveBeenCalledWith('surface');
    expect(checked()).toEqual(['false', 'false', 'false', 'true']);
    expect(stops()).toEqual([-1, -1, -1, 0]);
    expect(heading()).toBe('Surface');
  });

  it('the heading and the marks follow a value set from outside', () => {
    const { value, checked, heading } = setup();
    value.set('cartoon');
    flush();
    expect(heading()).toBe('Cartoon');
    expect(checked()).toEqual(['true', 'false', 'false', 'false']);
  });

  it('the arrows, Home and End move and pick, wrapping and skipping a disabled segment', () => {
    const { value, radios } = setup();
    key(radios[1], 'ArrowRight');
    expect(value.value).toBe('surface');
    expect(document.activeElement).toBe(radios[3]);
    key(radios[3], 'ArrowDown');
    expect(value.value).toBe('cartoon');
    key(radios[0], 'ArrowLeft');
    expect(value.value).toBe('surface');
    key(radios[3], 'ArrowUp');
    expect(value.value).toBe('ball-and-stick');
    key(radios[1], 'Home');
    expect(value.value).toBe('cartoon');
    key(radios[0], 'End');
    expect(value.value).toBe('surface');
  });

  it('a disabled segment is greyed, says why, and is never picked', () => {
    const { value, onChange, radios } = setup();
    expect(radios[2].getAttribute('aria-disabled')).toBe('true');
    expect(radios[2].dataset.tip).toBe('Too many atoms');
    radios[2].click();
    expect(value.value).toBe('ball-and-stick');
    expect(onChange).not.toHaveBeenCalled();
  });

  it('with no enabled segment checked, the first enabled one is the tab stop', () => {
    const { value, stops } = setup(state('spacefill'));
    expect(stops()).toEqual([0, -1, -1, -1]);
    value.set('surface');
    flush();
    expect(stops()).toEqual([-1, -1, -1, 0]);
  });

  it('a plain value is where it starts and it keeps the pick; a stream is followed, onChange updates its source', () => {
    const plain = vi.fn();
    const c1 = mount(() => <Segmented label="Type" options={types} value="cartoon" onChange={plain} />);
    const r1 = [...c1.querySelectorAll<HTMLButtonElement>('[role=radio]')];
    expect(r1[0].textContent).toBe('Cartoon');
    r1[1].click();
    expect(plain).toHaveBeenCalledWith('ball-and-stick');
    expect(r1[1].getAttribute('aria-checked')).toBe('true');

    const source = new BehaviorSubject('cartoon');
    const c2 = mount(() => (
      <Segmented aria-label="Type" options={types} value={source} onChange={(v) => source.next(v)} />
    ));
    const r2 = [...c2.querySelectorAll<HTMLButtonElement>('[role=radio]')];
    expect(c2.querySelector('[role=radiogroup]')!.getAttribute('aria-label')).toBe('Type');
    r2[3].click();
    flush();
    expect(source.value).toBe('surface');
    expect(r2[3].getAttribute('aria-checked')).toBe('true');
    expect(r2[0].getAttribute('aria-checked')).toBe('false');
  });

  it('a read-only state is not set by a pick', () => {
    const source = new State('cartoon');
    const onChange = vi.fn();
    const c = mount(() => (
      <Segmented aria-label="Type" options={types} value={source.asReadonly} onChange={onChange} />
    ));
    c.querySelectorAll<HTMLButtonElement>('[role=radio]')[3].click();
    expect(onChange).toHaveBeenCalledWith('surface');
    expect(source.value).toBe('cartoon');
  });
});
