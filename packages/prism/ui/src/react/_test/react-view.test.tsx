/** @vitest-environment jsdom */
/** @jsxImportSource react */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { afterEach, describe, expect, it } from 'vitest';
import { act } from 'react';
import { createRoot } from 'react-dom/client';
import { defineElements, MspElement } from '../../dom/element.js';
import { Scope } from '../../dom/scope.js';
import { observableState } from '../../reactive/observable-state.js';
import { MspElementHost } from '../element-host.js';
import { reactView } from '../react-view.js';
import { useObservable } from '../use-observable.js';

(globalThis as typeof globalThis & { IS_REACT_ACT_ENVIRONMENT: boolean }).IS_REACT_ACT_ENVIRONMENT = true;

class Greeting extends MspElement<{ name: string }> {
  static tag = 'msp-test-greeting';
  protected render(model: { name: string }) {
    const span = document.createElement('span');
    span.textContent = `Hello ${model.name}`;
    return span;
  }
}
defineElements([Greeting]);

afterEach(() => {
  document.body.innerHTML = '';
});

describe('React adapter', () => {
  it('mounts a React view over a view model and unmounts it', () => {
    const label = observableState('one');
    const View = ({ ctx }: { ctx: { label: typeof label } }) => <b>{useObservable(ctx.label.asReadonly)}</b>;
    const host = document.createElement('div');
    document.body.append(host);
    let dispose: () => void = () => {};
    act(() => {
      dispose = reactView(View, { createRoot }).mount(host, { label }, new Scope());
    });
    expect(host.textContent).toBe('one');
    act(() => label.set('two'));
    expect(host.textContent).toBe('two');
    act(() => dispose());
    expect(host.innerHTML).toBe('');
  });

  it('tells the view whether it is active, following the scope it is mounted in', () => {
    const View = ({ active }: { ctx: undefined; active?: boolean }) => <b>{String(active)}</b>;
    const host = document.createElement('div');
    document.body.append(host);
    const scope = new Scope();
    let dispose: () => void = () => {};
    act(() => {
      dispose = reactView(View, { createRoot }).mount(host, undefined, scope);
    });
    expect(host.textContent).toBe('true');
    act(() => scope.setActive(false));
    expect(host.textContent).toBe('false');
    act(() => scope.setActive(true));
    expect(host.textContent).toBe('true');
    act(() => dispose());
    scope.setActive(false);
    expect(host.innerHTML).toBe('');
  });

  it('hosts a prism-ui element with its model set as a property', () => {
    const container = document.createElement('div');
    document.body.append(container);
    const root = createRoot(container);
    const model = { name: 'Mol*' };
    act(() => root.render(<MspElementHost tag="msp-test-greeting" properties={{ model }} />));
    const el = container.querySelector('msp-test-greeting') as Greeting;
    expect(el.model).toBe(model);
    expect(el.textContent).toBe('Hello Mol*');
    act(() => root.unmount());
  });
});
