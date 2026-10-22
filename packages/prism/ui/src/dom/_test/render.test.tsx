/** @vitest-environment jsdom */
/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it, vi } from 'vitest';
import { Subject } from 'rxjs';
import { observableState } from '../../reactive/observable-state.js';
import type { Child } from '../jsx-types.js';
import { Show } from '../control-flow.js';
import { render } from '../render.js';
import { createRoot, getScope, onCleanup } from '../scope.js';

import { derive } from '../../reactive/state.js';
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

  it("sets a select's value after its options exist", () => {
    const choice = observableState('b');
    const { container } = mount(() => (
      <select value={choice}>
        <option value="a">A</option>
        <option value="b">B</option>
      </select>
    ));
    const select = container.querySelector('select')!;
    expect(select.value).toBe('b');
    choice.set('a');
    expect(select.value).toBe('a');
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
    const { html } = mount(() => <div>{derive(mode.asReadonly, view)}</div>);
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

describe('delegated events', () => {
  it('give each handler its own element as currentTarget, and stop at stopPropagation', () => {
    const container = document.createElement('div');
    const seen: string[] = [];
    render(
      () => (
        <div class="outer" onClick={(e) => seen.push(`outer:${(e.currentTarget as Element).className}`)}>
          <div class="middle" onClick={(e) => seen.push(`middle:${(e.currentTarget as Element).className}`)}>
            <button
              type="button"
              class="inner"
              onClick={(e) => seen.push(`inner:${(e.currentTarget as Element).className}`)}
            />
            <button
              type="button"
              class="stopper"
              onClick={(e) => {
                seen.push('stopper');
                e.stopPropagation();
              }}
            />
          </div>
        </div>
      ),
      container,
    );
    (container.querySelector('.inner') as HTMLElement).click();
    expect(seen).toEqual(['inner:inner', 'middle:middle', 'outer:outer']);
    seen.length = 0;
    (container.querySelector('.stopper') as HTMLElement).click();
    expect(seen).toEqual(['stopper']);
  });

  it('leaves the pointer moves native: no root walks a move, a handler gets its own (perf audit 11/P6)', () => {
    const container = document.createElement('div');
    document.body.append(container);
    const seen: string[] = [];
    render(
      () => (
        <div class="list" onPointerOver={(e) => seen.push(`over:${(e.currentTarget as Element).className}`)}>
          <button type="button" class="row" onPointerMove={(e) => seen.push(`move:${(e.currentTarget as Element).className}`)} />
          <span class="plain" />
        </div>
      ),
      container,
    );
    const paths = vi.spyOn(Event.prototype, 'composedPath');
    const plain = container.querySelector('.plain') as HTMLElement;
    plain.dispatchEvent(new Event('pointermove', { bubbles: true }));
    plain.dispatchEvent(new Event('mousemove', { bubbles: true }));
    expect(paths).not.toHaveBeenCalled();
    (container.querySelector('.row') as HTMLElement).dispatchEvent(new Event('pointermove', { bubbles: true }));
    plain.dispatchEvent(new Event('pointerover', { bubbles: true }));
    expect(seen).toEqual(['move:row', 'over:list']);
    expect(paths).not.toHaveBeenCalled();
    paths.mockRestore();
    container.remove();
  });

  it('runs each handler once across nested roots, and reaches what is rendered outside a root', () => {
    const outer = document.createElement('div');
    document.body.append(outer);
    const seen: string[] = [];
    render(() => <section onClick={() => seen.push('outer')} />, outer);
    const inner = outer.querySelector('section') as HTMLElement;
    render(() => <button type="button" onClick={() => seen.push('inner')} />, inner);
    (inner.querySelector('button') as HTMLElement).click();
    expect(seen).toEqual(['inner', 'outer']);
    // Made in a root, then moved under the body: the document carries on.
    seen.length = 0;
    const loose = document.createElement('div');
    render(() => <button type="button" class="loose" onClick={() => seen.push('loose')} />, loose);
    const button = loose.querySelector('.loose') as HTMLElement;
    document.body.append(button);
    button.click();
    expect(seen).toEqual(['loose']);
    button.remove();
    outer.remove();
  });

  it('keeps walking the path as dispatched when a handler removes its own subtree or an inner root goes', () => {
    const container = document.createElement('div');
    document.body.append(container);
    const open = observableState(true);
    const seen: string[] = [];
    render(
      () => (
        <div class="outer" onClick={() => seen.push('outer')}>
          <Show when={open}>
            {() => (
              <button
                type="button"
                class="closer"
                onClick={() => {
                  seen.push('closer');
                  // Closing the Show removes this button (and its parent chain) mid-dispatch.
                  open.set(false);
                }}
              />
            )}
          </Show>
        </div>
      ),
      container,
    );
    (container.querySelector('.closer') as HTMLElement).click();
    expect(seen).toEqual(['closer', 'outer']);

    // An inner root removed by its own handler: the outer root still continues above it.
    seen.length = 0;
    const host = container.querySelector('.outer') as HTMLElement;
    const inner = document.createElement('section');
    host.append(inner);
    render(
      () => (
        <button
          type="button"
          class="leaver"
          onClick={() => {
            seen.push('leaver');
            inner.remove();
          }}
        />
      ),
      inner,
    );
    (inner.querySelector('.leaver') as HTMLElement).click();
    expect(seen).toEqual(['leaver', 'outer']);
    container.remove();
  });

  it('give an event that does not bubble to its target only', () => {
    const container = document.createElement('div');
    const seen: string[] = [];
    render(
      () => (
        <div onInput={() => seen.push('parent')}>
          <input onInput={() => seen.push('input')} />
        </div>
      ),
      container,
    );
    container.querySelector('input')!.dispatchEvent(new Event('input'));
    expect(seen).toEqual(['input']);
    seen.length = 0;
    container.querySelector('input')!.dispatchEvent(new Event('input', { bubbles: true }));
    expect(seen).toEqual(['input', 'parent']);
  });
});
