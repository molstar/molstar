/** @jest-environment jsdom */
/** @jsxImportSource @molstar/plugin-ui-next/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { map, Subject } from 'rxjs';
import { observableState } from '../../reactive/observable-state.js';
import type { Child } from '../jsx-types.js';
import { render } from '../render.js';
import { createRoot, getScope, onCleanup } from '../scope.js';

function mount(fn: () => Child) {
  const container = document.createElement('div');
  const dispose = render(fn, container);
  return { container, dispose, html: () => container.innerHTML.replace(/<!---->/g, '') };
}

describe('elements', () => {
  it('creates real DOM with attributes, properties, text and nested children', () => {
    const { container, html } = mount(() => (
      <section id="s" class="a b" data-x="1" aria-label="Section" tabIndex={2}>
        <h1>Title {3}</h1>
        {[<i>x</i>, null, false, 'y']}
      </section>
    ));
    expect(html()).toBe(
      '<section id="s" class="a b" data-x="1" aria-label="Section" tabindex="2"><h1>Title 3</h1><i>x</i>y</section>',
    );
    expect((container.firstElementChild as HTMLElement).tabIndex).toBe(2);
  });

  it('sets input values and booleans as properties', () => {
    const { container } = mount(() => <input value="v" disabled={true} />);
    const input = container.querySelector('input')!;
    expect(input.value).toBe('v');
    expect(input.disabled).toBe(true);
  });

  it('creates SVG in the SVG namespace', () => {
    const { container } = mount(() => (
      <svg viewBox="0 0 24 24">
        <path d="M0 0" stroke-width={2} />
      </svg>
    ));
    const path = container.querySelector('path')!;
    expect(path.namespaceURI).toBe('http://www.w3.org/2000/svg');
    expect(path.getAttribute('stroke-width')).toBe('2');
    expect(container.querySelector('svg')!.getAttribute('viewBox')).toBe('0 0 24 24');
  });

  it('adds listeners for camel-case and exact event names, and calls ref with the element', () => {
    const seen: string[] = [];
    let ref: HTMLButtonElement | undefined;
    const { container } = mount(() => (
      <button
        onClick={() => seen.push('click')}
        onKeyDown={(e) => seen.push(e.key)}
        on:my-event={() => seen.push('custom')}
        ref={(el) => {
          ref = el;
        }}
      />
    ));
    const button = container.querySelector('button')!;
    button.click();
    button.dispatchEvent(new KeyboardEvent('keydown', { key: 'Enter' }));
    button.dispatchEvent(new Event('my-event'));
    expect(seen).toEqual(['click', 'Enter', 'custom']);
    expect(ref).toBe(button);
  });

  it('renders fragments and components', () => {
    const Item = (props: { label: string; children?: Child }) => (
      <>
        <b>{props.label}</b>
        {props.children}
      </>
    );
    const { html } = mount(() => (
      <p>
        <Item label="a">!</Item>
        <Item label="b" />
      </p>
    ));
    expect(html()).toBe('<p><b>a</b>!<b>b</b></p>');
  });
});

describe('bindings', () => {
  it('updates text, attributes, classes and styles in place when states change', () => {
    const label = observableState('one');
    const active = observableState(false);
    const color = observableState<string | undefined>('red');
    const { container } = mount(() => (
      <div
        title={label}
        class="row"
        classList={{ active, fixed: true }}
        style={{ color, '--msp-gap': '4px', marginTop: 2 }}
      >
        {label}
      </div>
    ));
    const div = container.querySelector('div')!;
    const text = div.firstChild?.nextSibling;
    expect(div.title).toBe('one');
    expect(div.className).toBe('row fixed');
    expect(div.style.color).toBe('red');
    expect(div.style.getPropertyValue('--msp-gap')).toBe('4px');
    expect(div.style.marginTop).toBe('2px');
    label.set('two');
    active.set(true);
    color.set(undefined);
    expect(div.title).toBe('two');
    expect(div.textContent).toBe('two');
    // The same text node, updated in place.
    expect(div.firstChild?.nextSibling).toBe(text);
    expect(div.className).toBe('row active fixed');
    expect(div.style.color).toBe('');
  });

  it('replaces non-text content, rendering a function value in its own scope', () => {
    const mode = observableState<'a' | 'b'>('a');
    const cleaned: string[] = [];
    const view = (m: string) => () => {
      onCleanup(() => cleaned.push(m));
      return <em>{m}</em>;
    };
    const { html } = mount(() => <div>{mode.asReadonly.pipe(map(view))}</div>);
    expect(html()).toBe('<div><em>a</em></div>');
    mode.set('b');
    expect(html()).toBe('<div><em>b</em></div>');
    expect(cleaned).toEqual(['a']);
  });

  it('unsubscribes everything when the root is disposed', () => {
    const source = new Subject<string>();
    const label = observableState('x');
    const { html, dispose } = mount(() => (
      <p title={label}>
        {source}
        {label}
      </p>
    ));
    source.next('s');
    expect(html()).toBe('<p title="x">sx</p>');
    expect(source.observed).toBe(true);
    dispose();
    expect(source.observed).toBe(false);
    label.set('y');
    expect(html()).toBe('');
  });

  it('refuses bindings outside a render scope instead of leaking them', () => {
    expect(getScope()).toBeUndefined();
    expect(() => <p>{observableState('x')}</p>).toThrow(/outside a render scope/);
    // A root makes one.
    createRoot((dispose) => {
      expect(() => <p>{observableState('x')}</p>).not.toThrow();
      dispose();
    });
  });
});
