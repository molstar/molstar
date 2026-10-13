/** @vitest-environment jsdom */
/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { afterEach, describe, expect, it } from 'vitest';
import { type ObservableState, observableState } from '../../reactive/observable-state.js';
import { ReactiveModel } from '../../reactive/reactive-model.js';
import { defineElements, MspElement } from '../element.js';
import { render } from '../render.js';
import { onCleanup } from '../scope.js';
import { domView, elementView } from '../view.js';

interface Counter {
  count: ObservableState<number>;
  log: string[];
}

class TestCounter extends MspElement<Counter> {
  static tag = 'msp-test-counter';
  protected render(model: Counter) {
    model.log.push('render');
    onCleanup(() => model.log.push('dispose'));
    return <span>{model.count.asReadonly}</span>;
  }
}

class OwnedModel extends ReactiveModel {
  readonly events: string[] = [];
  protected onMount() {
    this.events.push('mount');
  }
  protected onUnmount() {
    this.events.push('unmount');
  }
}

class TestOwner extends MspElement<OwnedModel> {
  static tag = 'msp-test-owner';
  protected createModel() {
    return new OwnedModel();
  }
  protected render() {
    return <b>owned</b>;
  }
}

declare module '../jsx-types.js' {
  interface MspElements {
    'msp-test-counter': import('../element.js').MspElementProps<TestCounter>;
  }
}

defineElements([TestCounter, TestOwner]);
const tick = () => new Promise<void>((resolve) => queueMicrotask(resolve));

function counter(): Counter {
  return { count: observableState(1), log: [] };
}

afterEach(() => {
  document.body.innerHTML = '';
});

describe('MspElement', () => {
  it('renders its model once connected, and bindings update it', () => {
    const model = counter();
    const el = document.createElement('msp-test-counter') as TestCounter;
    el.model = model;
    expect(el.innerHTML).toBe('');
    document.body.append(el);
    expect(el.textContent).toBe('1');
    model.count.set(2);
    expect(el.textContent).toBe('2');
    expect(model.log).toEqual(['render']);
  });

  it('keeps its rendering when moved, and disposes it when removed', async () => {
    const model = counter();
    const a = document.createElement('div');
    const b = document.createElement('div');
    document.body.append(a, b);
    const el = document.createElement('msp-test-counter') as TestCounter;
    el.model = model;
    a.append(el);
    const span = el.firstChild;
    b.append(el);
    await tick();
    expect(el.firstChild).toBe(span);
    expect(model.log).toEqual(['render']);
    el.remove();
    await tick();
    expect(model.log).toEqual(['render', 'dispose']);
    expect(el.innerHTML).toBe('');
  });

  it('renders again for another model', () => {
    const first = counter();
    const second = counter();
    second.count.set(5);
    const el = document.createElement('msp-test-counter') as TestCounter;
    el.model = first;
    document.body.append(el);
    el.model = second;
    expect(el.textContent).toBe('5');
    expect(first.log).toEqual(['render', 'dispose']);
  });

  it('takes a model set before the element was defined', () => {
    class Late extends MspElement<Counter> {
      static tag = 'msp-test-late';
      protected render(model: Counter) {
        return <i>{model.count.asReadonly}</i>;
      }
    }
    const el = document.createElement('msp-test-late') as any;
    el.model = counter();
    document.body.append(el);
    defineElements([Late]);
    expect(el.textContent).toBe('1');
  });

  it('mounts a model it created while rendered', async () => {
    const el = document.createElement('msp-test-owner') as TestOwner;
    document.body.append(el);
    const model = el.model!;
    expect(model.events).toEqual(['mount']);
    el.remove();
    await tick();
    expect(model.events).toEqual(['mount', 'unmount']);
  });

  it('takes its model as a property from JSX', () => {
    const model = counter();
    const container = document.createElement('div');
    document.body.append(container);
    render(() => <msp-test-counter model={model} />, container);
    expect(container.querySelector('msp-test-counter')!.textContent).toBe('1');
  });
});

describe('view renderers', () => {
  it('mounts and disposes a JSX view or an element view', async () => {
    const host = document.createElement('div');
    document.body.append(host);
    const cleaned: string[] = [];
    const view = domView<{ name: string }>(({ ctx }) => {
      onCleanup(() => cleaned.push(ctx.name));
      return <p>{ctx.name}</p>;
    });
    const dispose = view.mount(host, { name: 'tab' });
    expect(host.textContent).toBe('tab');
    dispose();
    expect(host.innerHTML).toBe('');
    expect(cleaned).toEqual(['tab']);

    const model = counter();
    const disposeElement = elementView<Counter>('msp-test-counter', 'model').mount(host, model);
    expect(host.textContent).toBe('1');
    disposeElement();
    await tick();
    expect(host.innerHTML).toBe('');
    expect(model.log).toEqual(['render', 'dispose']);
  });
});
