/**
 * @jest-environment jsdom
 * @jest-environment-options {"customExportConditions": ["require", "default"]}
 */
/** @jsxImportSource preact */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { createRoot } from 'preact/compat/client';
import { act } from 'preact/test-utils';
import { defineElements, MspElement } from '../../dom/element.js';
import { observableState } from '../../reactive/observable-state.js';

// What a Preact app's bundler does (`react` aliased to `preact/compat`): the adapter then runs on Preact as written.
// The repo's transform does not hoist `jest.mock`, so the adapter is required after it.
jest.mock('react', () => jest.requireActual('preact/compat'));
const { MspElementHost } = require('../element-host.js') as typeof import('../element-host.js');
const { reactView } = require('../react-view.js') as typeof import('../react-view.js');
const { useObservable } = require('../use-observable.js') as typeof import('../use-observable.js');
const { useObservableProp } = require('../use-observable-prop.js') as typeof import('../use-observable-prop.js');

class Greeting extends MspElement<{ name: string }> {
  static tag = 'msp-test-preact-greeting';
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

describe('React adapter on Preact (react aliased to preact/compat)', () => {
  it('mounts a Preact view over a view model, follows it, and unmounts it', () => {
    const label = observableState('one');
    const View = ({ ctx }: { ctx: { label: typeof label } }) => <b>{useObservable(ctx.label.asReadonly)}</b>;
    const host = document.createElement('div');
    document.body.append(host);
    let dispose: () => void = () => {};
    act(() => {
      dispose = reactView(View as never, { createRoot: createRoot as never }).mount(host, { label });
    });
    expect(host.textContent).toBe('one');
    act(() => label.set('two'));
    expect(host.textContent).toBe('two');
    act(() => dispose());
    expect(host.innerHTML).toBe('');
  });

  it('re-renders for the selected part only', () => {
    const record = observableState({ name: 'a', count: 0 });
    let renders = 0;
    const View = () => {
      renders++;
      return <i>{useObservableProp(record.asReadonly, (r) => r.name)}</i>;
    };
    const host = document.createElement('div');
    document.body.append(host);
    const root = createRoot(host);
    act(() => root.render(<View />));
    const first = renders;
    act(() => record.set({ name: 'a', count: 1 }));
    expect(renders).toBe(first);
    act(() => record.set({ name: 'b', count: 1 }));
    expect(host.textContent).toBe('b');
    act(() => root.unmount());
  });

  it('hosts a prism-ui element with its model set as a property', () => {
    const container = document.createElement('div');
    document.body.append(container);
    const root = createRoot(container);
    const model = { name: 'Mol*' };
    act(() => root.render(<MspElementHost tag="msp-test-preact-greeting" properties={{ model }} />));
    const el = container.querySelector('msp-test-preact-greeting') as Greeting;
    expect(el.model).toBe(model);
    expect(el.textContent).toBe('Hello Mol*');
    act(() => root.unmount());
  });
});
