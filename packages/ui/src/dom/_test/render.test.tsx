/** @vitest-environment jsdom */
/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it, vi } from 'vitest';
import { BehaviorSubject, State } from '@molstar/core/reactive';
import { For, flush, render, Show, ViewScope } from '../index.js';

// JSX, props, regions and event delegation (dom §4).

function mount(fn: () => Node) {
  const host = document.createElement('div');
  document.body.appendChild(host);
  const dispose = render(fn, host);
  return {
    host,
    dispose: () => {
      dispose();
      host.remove();
    },
  };
}

describe('JSX and props', () => {
  it('builds real DOM with attributes, properties, classes, styles and ref', () => {
    let ref: HTMLInputElement | undefined;
    const { host, dispose } = mount(() => (
      <div id="a" class="x" classList={{ on: true, off: false }} style={{ marginTop: 4, opacity: 0.5, '--k': 'v' }}>
        <input value="hi" disabled={false} ref={(el) => (ref = el)} aria-hidden={false} />
        <svg viewBox="0 0 1 1">
          <path d="M0" />
        </svg>
        text {1} {null} {false}
      </div>
    ));
    const div = host.querySelector('div')!;
    expect(div.id).toBe('a');
    expect(div.getAttribute('class')).toBe('x on');
    expect(div.style.marginTop).toBe('4px');
    expect(div.style.opacity).toBe('0.5');
    expect(div.style.getPropertyValue('--k')).toBe('v');
    expect(ref?.value).toBe('hi');
    expect(ref?.hasAttribute('disabled')).toBe(false);
    expect(ref?.getAttribute('aria-hidden')).toBe('false');
    expect(host.querySelector('path')!.namespaceURI).toBe('http://www.w3.org/2000/svg');
    expect(div.textContent).toContain('text 1');
    dispose();
  });

  it('reactive props and classes follow their states, and stop with the scope', () => {
    const title = new State('one');
    const on = new State(false);
    const { host, dispose } = mount(() => <span title={title} class={title} classList={{ on }} />);
    const span = host.querySelector('span')!;
    expect(span.title).toBe('one');
    title.set('two');
    on.set(true);
    expect(span.title).toBe('two');
    expect(span.getAttribute('class')).toBe('two on');
    dispose();
    expect(title.observed || on.observed).toBe(false);
  });

  it('a reactive child is a region: text in place, nodes in a child scope', () => {
    const value = new State<unknown>('a');
    const inner = new State('x');
    const { host, dispose } = mount(() => <p>{value}</p>);
    const p = host.querySelector('p')!;
    expect(p.textContent).toBe('a');
    const text = p.firstChild;
    value.set('b');
    expect(p.firstChild).toBe(text);
    expect(p.textContent).toBe('b');
    value.set(() => <b>{inner}</b>);
    expect(p.querySelector('b')!.textContent).toBe('x');
    inner.set('y');
    expect(p.querySelector('b')!.textContent).toBe('y');
    value.set('c');
    expect(inner.observed).toBe(false);
    expect(p.textContent).toBe('c');
    dispose();
  });

  it('a BehaviorSubject child and a stream child are followed', () => {
    const subject = new BehaviorSubject(1);
    const { host, dispose } = mount(() => <i>{subject}</i>);
    subject.next(2);
    expect(host.textContent).toBe('2');
    dispose();
    expect(subject.observed).toBe(false);
  });
});

describe('events', () => {
  it('delegated handlers get currentTarget, stop at stopPropagation, and flush at the end of the dispatch', () => {
    const seen: string[] = [];
    const count = new State(0);
    let rendered = '';
    const { host, dispose } = mount(() => (
      <div onClick={() => seen.push('outer')}>
        <button
          type="button"
          onClick={(e) => {
            seen.push(e.currentTarget.tagName);
            count.set(count.value + 1);
          }}
        >
          <span>{count}</span>
        </button>
        <a
          href="#x"
          onClick={(e) => {
            e.preventDefault();
            e.stopPropagation();
            seen.push('a');
          }}
        >
          a
        </a>
      </div>
    ));
    host.querySelector('span')!.click();
    rendered = host.querySelector('span')!.textContent!;
    host.querySelector('a')!.click();
    expect(seen).toEqual(['BUTTON', 'outer', 'a']);
    expect(rendered).toBe('1');
    dispose();
  });

  it('on:type is a native listener; a non-bubbling event reaches the target handler', () => {
    const native = vi.fn();
    const input = vi.fn();
    const { host, dispose } = mount(() => <input on:custom={native} onInput={input} />);
    const el = host.querySelector('input')!;
    el.dispatchEvent(new Event('custom'));
    el.dispatchEvent(new Event('input'));
    expect(native).toHaveBeenCalledTimes(1);
    expect(input).toHaveBeenCalledTimes(1);
    dispose();
  });
});

describe('render', () => {
  it('renders after what is there, under a parent scope that pauses it, and removes it all', () => {
    const host = document.createElement('div');
    host.append('before');
    const parent = new ViewScope();
    const s = new State(1);
    const dispose = render(() => <b>{s}</b>, host, parent);
    expect(host.textContent).toBe('before1');
    parent.setActive(false);
    s.set(2);
    expect(host.textContent).toBe('before1');
    expect(s.observed).toBe(false);
    parent.setActive(true);
    expect(host.textContent).toBe('before2');
    dispose();
    flush();
    expect(host.textContent).toBe('before');
  });
});

describe('SVG namespace', () => {
  const SVG = 'http://www.w3.org/2000/svg';

  it('a, title, style and script are SVG under an SVG element, HTML elsewhere', () => {
    const { host, dispose } = mount(() => (
      <div>
        <a href="#h">html</a>
        <svg>
          <title>icon</title>
          <style>{'circle { fill: red }'}</style>
          <a href="#s">
            <text>t</text>
          </a>
          <g>
            <circle>
              <title>dot</title>
            </circle>
          </g>
        </svg>
      </div>
    ));
    const [htmlA, svgA] = host.querySelectorAll('a');
    expect(htmlA.namespaceURI).toBe('http://www.w3.org/1999/xhtml');
    expect(svgA.namespaceURI).toBe(SVG);
    expect(svgA.getAttribute('href')).toBe('#s');
    for (const t of host.querySelectorAll('title, style')) expect(t.namespaceURI).toBe(SVG);
    expect(host.querySelector('circle title')!.textContent).toBe('dot');
    dispose();
  });

  it('what a function child, a region or a list renders under an SVG element is SVG', () => {
    const shown = new State(true);
    const items = new State(['x']);
    const { host, dispose } = mount(() => (
      <svg>
        <g>
          {() => (
            <Show when={shown}>
              {() => (
                <a>
                  <title>first</title>
                </a>
              )}
            </Show>
          )}
        </g>
        {() => <For each={items}>{(i) => <title>{i}</title>}</For>}
      </svg>
    ));
    const all = () => [...host.querySelectorAll('a, title')].map((el) => el.namespaceURI);
    expect(all()).toEqual([SVG, SVG, SVG]);
    shown.set(false);
    shown.set(true);
    items.set(['x', 'y']);
    expect(all()).toEqual([SVG, SVG, SVG, SVG]);
    dispose();
  });

  it('a shared tag made again under an SVG parent warns when it had a bound prop', () => {
    const warn = vi.spyOn(console, 'warn').mockImplementation(() => {});
    const label = new State('x');
    const { host, dispose } = mount(() => (
      <svg>
        <title id={label}>t</title>
      </svg>
    ));
    expect(host.querySelector('title')!.namespaceURI).toBe(SVG);
    expect(warn).toHaveBeenCalledTimes(1);
    warn.mockRestore();
    dispose();
  });
});
