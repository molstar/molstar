/** @jest-environment jsdom */
/** @jsxImportSource react */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { act, StrictMode } from 'react';
import { createRoot, type Root } from 'react-dom/client';
import { observableState } from '../../reactive/observable-state.js';
import { ReactiveModel } from '../../reactive/reactive-model.js';
import type { ReadonlyState } from '../../reactive/readonly-state.js';
import { useObservable } from '../use-observable.js';
import { useObservableProp } from '../use-observable-prop.js';
import { useReactiveModel } from '../use-reactive-model.js';

(globalThis as typeof globalThis & { IS_REACT_ACT_ENVIRONMENT: boolean }).IS_REACT_ACT_ENVIRONMENT = true;

let container: HTMLDivElement;
let root: Root;

beforeEach(() => {
  container = document.createElement('div');
  document.body.append(container);
  root = createRoot(container);
});

afterEach(async () => {
  await act(async () => root.unmount());
  container.remove();
});

function Readout({ source }: { source: ReadonlyState<string> }) {
  const value = useObservable(source);
  return <span>{value}</span>;
}

describe('useObservable', () => {
  it('renders the current value and updates after replacement and remount', async () => {
    const writable = observableState('one');
    await act(async () => {
      root.render(<Readout source={writable.asReadonly} />);
    });
    expect(container.textContent).toBe('one');
    await act(async () => {
      writable.set('two');
    });
    expect(container.textContent).toBe('two');
    await act(async () => {
      root.render(<Readout source={writable.asReadonly} />);
    });
    expect(container.textContent).toBe('two');
  });

  it('does not reconcile again from the initial behavior emission', async () => {
    const writable = observableState('one');
    let renders = 0;
    function Probe() {
      renders += 1;
      const value = useObservable(writable.asReadonly);
      return <span>{value}</span>;
    }
    await act(async () => {
      root.render(<Probe />);
    });
    expect(container.textContent).toBe('one');
    expect(renders).toBe(1);
  });
});

describe('useObservableProp', () => {
  it('re-renders only when the selected slice changes', async () => {
    const source = observableState({ hot: 0, cold: 'a' });
    let renders = 0;
    function Probe() {
      renders += 1;
      const cold = useObservableProp(source.asReadonly, (value) => value.cold);
      return <span>{cold}</span>;
    }
    await act(async () => {
      root.render(<Probe />);
    });
    expect(renders).toBe(1);
    await act(async () => {
      source.set({ hot: 1, cold: 'a' });
    });
    expect(renders).toBe(1);
    await act(async () => {
      source.set({ hot: 2, cold: 'b' });
    });
    expect(renders).toBe(2);
    expect(container.textContent).toBe('b');
  });

  it('uses shallow compare for selected object slices', async () => {
    const list = { available: [] as readonly string[], conversationId: undefined as string | undefined };
    const source = observableState({ list, view: { count: 0 } });
    let renders = 0;
    function Probe() {
      renders += 1;
      const selected = useObservableProp(source.asReadonly, (value) => value.list);
      return <span>{selected.conversationId ?? 'none'}</span>;
    }
    await act(async () => {
      root.render(<Probe />);
    });
    expect(renders).toBe(1);
    await act(async () => {
      source.set({ list, view: { count: 1 } });
    });
    expect(renders).toBe(1);
    await act(async () => {
      source.set({ list: { ...list, conversationId: 'c-1' }, view: { count: 2 } });
    });
    expect(renders).toBe(2);
    expect(container.textContent).toBe('c-1');
  });
});

describe('useReactiveModel', () => {
  it('never retires a model during real Strict Mode effect replay or component cleanup', async () => {
    const order: string[] = [];
    // A stand-in for an expensive long-lived resource, such as a Mol* plugin context.
    const resource = { dispose: () => order.push('resource dispose') };
    class Example extends ReactiveModel {
      constructor() {
        super();
        this.registerPermanentDispose(() => order.push('permanent cleanup'));
      }

      protected override onMount(): void {
        order.push('mount');
        this.registerDispose(() => order.push('mount cleanup'));
      }

      protected override onRetire(): void {
        order.push('retire');
        resource.dispose();
      }
    }
    const model = new Example();
    function Probe() {
      useReactiveModel(model);
      return null;
    }

    await act(async () =>
      root.render(
        <StrictMode>
          <Probe />
        </StrictMode>,
      ),
    );
    expect(order).toEqual(['mount', 'mount cleanup', 'mount']);
    await act(async () => root.render(null));
    expect(order).toEqual(['mount', 'mount cleanup', 'mount', 'mount cleanup']);
    await act(async () =>
      root.render(
        <StrictMode>
          <Probe />
        </StrictMode>,
      ),
    );
    expect(order).toEqual(['mount', 'mount cleanup', 'mount', 'mount cleanup', 'mount', 'mount cleanup', 'mount']);
    await act(async () => root.render(null));
    expect(order).not.toContain('retire');
    expect(order).not.toContain('resource dispose');
    expect(order).not.toContain('permanent cleanup');

    // Only the creator ends permanent ownership, after the final host release.
    model.retire();
    model.retire();
    model.mount();
    expect(order.slice(-4)).toEqual(['mount cleanup', 'retire', 'resource dispose', 'permanent cleanup']);
    expect(order.filter((step) => step === 'retire')).toHaveLength(1);
    expect(order.filter((step) => step === 'resource dispose')).toHaveLength(1);
  });

  it('mounts and unmounts across StrictMode-style remounts', async () => {
    const order: string[] = [];
    class Example extends ReactiveModel {
      protected override onMount(): void {
        order.push('mount');
        this.registerDispose(() => {
          order.push('unmount');
        });
      }
    }
    const model = new Example();
    function Probe() {
      useReactiveModel(model);
      return null;
    }
    await act(async () => {
      root.render(<Probe />);
    });
    await act(async () => {
      root.unmount();
    });
    root = createRoot(container);
    await act(async () => {
      root.render(<Probe />);
    });
    expect(order).toEqual(['mount', 'unmount', 'mount']);
  });

  it('returns undefined for nullish models', async () => {
    let seen: string | undefined = 'pending';
    function Probe({ model }: { model: ReactiveModel | null }) {
      const value = useReactiveModel(model);
      seen = value === undefined ? 'undefined' : 'present';
      return null;
    }
    await act(async () => {
      root.render(<Probe model={null} />);
    });
    expect(seen).toBe('undefined');
  });
});
