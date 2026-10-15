/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Fixtures for the React UI's component tests (UX3): a real context over a document by value (no network, no `node:fs`
// in the browser), the manual clock, and a canvas stub in place of `mountSection` (a focusable element sized like a
// surface, the section's viewport sized as the engine would).

import { act, type ComponentType, Profiler, type ReactNode, StrictMode } from 'react';
import { createRoot } from 'react-dom/client';
import { afterEach } from 'vitest';
import type { SeqViewerContext } from '../../src/index.js';
import {
  SeqViewerLayout,
  SeqViewerProvider,
  type SeqViewerUISpec,
  type SlotComponents,
  type SlotName,
} from '../../src/react/index.js';
import { openedContext, stubSurface, WIDTH } from '../dom/fixtures.js';

// The UI fixture, the canvas stub and pointer records are the DOM tests' too (_test/dom/fixtures.ts).
export {
  at,
  HEIGHT,
  headlessSpec,
  keys,
  lines,
  openedContext,
  pointer,
  stubSurface,
  subscriptions,
  uiDoc,
  WIDTH,
  xOf,
} from '../dom/fixtures.js';

(globalThis as { IS_REACT_ACT_ENVIRONMENT?: boolean }).IS_REACT_ACT_ENVIRONMENT = true;

/** Render counts per slot: each default (or given) component under a counting profiler. */
export function counted(defaults: Partial<SlotComponents>): {
  components: SeqViewerUISpec['components'];
  counts: Record<string, number>;
  reset(): void;
} {
  const counts: Record<string, number> = {};
  const components: Record<string, unknown> = {};
  for (const [name, Part] of Object.entries(defaults) as [SlotName, ComponentType<object>][]) {
    counts[name] = 0;
    // A profiler counts every commit of the part's subtree, its own state updates included.
    components[name] = (props: object) => (
      <Profiler id={name} onRender={() => counts[name]++}>
        <Part {...props} />
      </Profiler>
    );
  }
  return {
    components: components as SeqViewerUISpec['components'],
    counts,
    reset: () => {
      for (const k in counts) counts[k] = 0;
    },
  };
}

export async function mount(ctx: SeqViewerContext, ui: SeqViewerUISpec, strict = false) {
  const container = document.createElement('div');
  container.style.width = `${WIDTH}px`;
  document.body.appendChild(container);
  const root = createRoot(container);
  const tree = (
    <SeqViewerProvider context={ctx} ui={ui}>
      <SeqViewerLayout />
    </SeqViewerProvider>
  );
  await act(async () => root.render(strict ? <StrictMode>{tree}</StrictMode> : tree));
  await settle();
  return {
    container,
    unmount: async () => {
      await act(async () => root.unmount());
      container.remove();
    },
  };
}

/** Lets pending graph work, state deliveries and React commits finish. */
export async function settle() {
  for (let i = 0; i < 3; i++) await act(async () => new Promise((r) => setTimeout(r, 0)));
}

/**
 * A harness for a test file: `opened(ui)` makes the UI fixture's context and mounts the default arrangement on it in a
 * container, both cleaned up after each test (`cleanups` takes more).
 */
export function harness() {
  const cleanups: (() => Promise<void> | void)[] = [];
  afterEach(async () => {
    for (const c of cleanups.splice(0).reverse()) await c();
  });
  return {
    cleanups,
    /** A root that renders what `view(i)` returns, again on each `render()`: a parent that re-renders (an inline `ui`). */
    parent(view: (i: number) => ReactNode, strict = false) {
      const el = document.createElement('div');
      el.style.width = `${WIDTH}px`;
      document.body.appendChild(el);
      const root = createRoot(el);
      let i = 0;
      const render = async () => {
        const tree = view(i++);
        await act(async () => root.render(strict ? <StrictMode>{tree}</StrictMode> : tree));
        await settle();
      };
      cleanups.push(async () => {
        await act(async () => root.unmount());
        el.remove();
      });
      return { render, root, container: el };
    },
    async opened(ui: SeqViewerUISpec = {}, strict = false) {
      const t = await openedContext();
      const m = await mount(t.ctx, { surface: stubSurface, ...ui }, strict);
      cleanups.push(() => t.ctx.dispose(), m.unmount);
      return { ...t, ...m };
    },
  };
}
