/** @jsxImportSource react */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type CSSProperties, version } from 'react';
import type { ObservableState } from '../../reactive/observable-state.js';
import { MspElementHost } from '../element-host.js';
import { useObservable } from '../use-observable.js';

// The React side of Internals › React and Preact: a panel written in React over a prism-ui state, and a React app that
// hosts prism-ui in an element.

export const panelStyle: CSSProperties = {
  display: 'grid',
  gap: 8,
  padding: 'var(--msp-space-3) var(--msp-space-4)',
  fontSize: 'var(--msp-font-size-sm)',
  color: 'var(--msp-color-text)',
};

export const buttonStyle: CSSProperties = {
  justifySelf: 'start',
  height: 28,
  padding: '0 var(--msp-space-4)',
  border: 0,
  borderRadius: 'var(--msp-radius-md)',
  background: 'var(--msp-color-control)',
  color: 'inherit',
  font: 'inherit',
};

/** A tab written in React: it reads a prism-ui state with `useObservable` and writes to it. */
export function ReactPanel(props: { count: ObservableState<number>; post(text: string): void }) {
  const count = useObservable(props.count.asReadonly);
  return (
    <div style={panelStyle}>
      <p style={{ margin: 0, color: 'var(--msp-color-text-muted)' }}>
        This tab is a React component (React {version}), mounted by `reactView`. The count is a prism-ui state, shared
        with the Preact tab.
      </p>
      <strong>Count: {count}</strong>
      <button
        type="button"
        style={buttonStyle}
        onClick={() => {
          props.count.set(count + 1);
          props.post(`React added one (${count + 1})`);
        }}
      >
        Add one
      </button>
    </div>
  );
}

/** A React app with prism-ui inside it: a React header over the same state, the shell in a prism-ui element. */
export function ReactApp(props: { count: ObservableState<number>; shell: unknown }) {
  const count = useObservable(props.count.asReadonly);
  return (
    <div style={{ display: 'grid', gridTemplateRows: 'auto 1fr', height: '100vh', background: 'var(--msp-color-bg)' }}>
      <header
        style={{
          ...panelStyle,
          gridAutoFlow: 'column',
          justifyContent: 'start',
          alignItems: 'center',
          borderBottom: '1px solid var(--msp-color-border)',
        }}
      >
        <strong>A React app</strong>
        <span style={{ color: 'var(--msp-color-text-muted)' }}>its header is React; below is prism-ui</span>
        <span>Count: {count}</span>
        <button type="button" style={buttonStyle} onClick={() => props.count.set(count + 1)}>
          Add one
        </button>
      </header>
      <MspElementHost tag="msp-story-root" properties={{ model: props.shell }} style={{ minHeight: 0 }} />
    </div>
  );
}
