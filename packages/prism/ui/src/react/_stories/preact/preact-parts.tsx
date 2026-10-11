/** @jsxImportSource preact */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { createRoot } from 'preact/compat/client';
import type { ViewRenderer } from '../../../dom/view.js';
import type { ObservableState } from '../../../reactive/observable-state.js';
import { reactView } from '../../react-view.js';
import { useObservable } from '../../use-observable.js';

// The Preact side of Internals › React and Preact. Storybook's config resolves `react` to `preact/compat` for the
// modules of this folder and the adapter modules they import (a copy of each, so the React story keeps real React), as
// a Preact app's bundler does for the whole app.

const panel = {
  display: 'grid',
  gap: '8px',
  padding: 'var(--msp-space-3) var(--msp-space-4)',
  fontSize: 'var(--msp-font-size-sm)',
  color: 'var(--msp-color-text)',
};

const button = {
  justifySelf: 'start',
  height: '28px',
  padding: '0 12px',
  border: 0,
  borderRadius: 'var(--msp-radius-md)',
  background: 'var(--msp-color-control)',
  color: 'inherit',
  font: 'inherit',
};

interface PanelProps {
  count: ObservableState<number>;
  post(text: string): void;
}

/** A tab written in Preact: the same adapter (`useObservable`, `reactView`) over the same prism-ui state. */
function PreactPanel(props: PanelProps) {
  const count = useObservable(props.count.asReadonly);
  return (
    <div style={panel}>
      <p style={{ margin: 0, color: 'var(--msp-color-text-muted)' }}>
        This tab is a Preact component, mounted by `reactView` with Preact's `createRoot`. The count is a prism-ui
        state, shared with the React tab.
      </p>
      <strong>Count: {count}</strong>
      <button
        type="button"
        style={button}
        onClick={() => {
          props.count.set(count + 1);
          props.post(`Preact added one (${count + 1})`);
        }}
      >
        Add one
      </button>
    </div>
  );
}

/** The Preact tab's view, for a shell tab's `render`. */
export function preactPanelView(props: PanelProps): ViewRenderer<unknown> {
  // `reactView` and `createRoot` are typed for React; here both are Preact's.
  return reactView((() => <PreactPanel {...props} />) as never, { createRoot: createRoot as never });
}
