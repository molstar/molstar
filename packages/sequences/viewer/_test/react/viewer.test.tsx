/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The provider and `SeqViewer` (UX3; seqspec:interaction#the-default-ui; UI1's React cases): an inline `ui` keeps every
// mounted surface, its focus and the parts' renders; a changed slot or chrome remakes the scope but not the surfaces, a
// new `surface` does; `{ Menu: undefined }` leaves the default menu; `SeqViewer` makes a context per mount and disposes
// it, applies `document` and `data` by identity, one effect each, and reads `spec` and `ports` at mount.

import { act } from 'react';
import { describe, expect, it, vi } from 'vitest';
import '../../src/skin.css';
import { SeqViewerConfig, type SeqViewerContext } from '../../src/index.js';
import {
  Legend,
  SeqViewer,
  SeqViewerLayout,
  SeqViewerProvider,
  type SeqViewerUISpec,
  useSeqViewer,
} from '../../src/react/index.js';
import { fixtureMeasurer } from '../../src/testing/index.js';
import { harness, headlessSpec, openedContext, settle, stubSurface, uiDoc } from './fixtures.js';

const { cleanups, parent } = harness();

/** The canvas stub, counting what mounts and what is disposed. */
function countingSurface() {
  const n = { mounted: 0, disposed: 0 };
  const surface: SeqViewerUISpec['surface'] = (ctx, section, host) => {
    n.mounted++;
    const s = stubSurface(ctx, section, host);
    return {
      element: s.element,
      dispose: () => {
        n.disposed++;
        s.dispose();
      },
    };
  };
  return { n, surface };
}
async function until(cond: () => unknown, ms = 4000) {
  for (let t = 0; t < ms && !cond(); t += 25) await act(async () => new Promise((r) => setTimeout(r, 25)));
  expect(cond()).toBeTruthy();
}
const surfaceEl = () => document.querySelector('.seqv-surface [tabindex]') as HTMLElement;

describe('the provider keeps its surfaces', () => {
  it("an inline `ui` object, with inline `components`, keeps every surface, its focus and the parts' renders", async () => {
    const { ctx } = await openedContext();
    cleanups.push(() => ctx.dispose());
    const { n, surface } = countingSurface();
    // What a part reads from the provider: one object while the `ui` is equal in two levels.
    const seen = new Set<unknown>();
    const Probe = () => {
      seen.add(useSeqViewer());
      return null;
    };
    const Legend2 = (p: { context: SeqViewerContext }) => <Legend {...p} />;
    const p = parent(() => (
      // New objects each render: the common case of a host component.
      <SeqViewerProvider context={ctx} ui={{ surface, className: 'host', components: { Legend: Legend2 } }}>
        <Probe />
        <SeqViewerLayout />
      </SeqViewerProvider>
    ));
    await p.render();
    expect([n.mounted, n.disposed]).toEqual([1, 0]);
    surfaceEl().focus();
    const was = surfaceEl();
    for (let k = 0; k < 3; k++) await p.render();
    expect([n.mounted, n.disposed, surfaceEl() === was, document.activeElement === was]).toEqual([1, 0, true, true]);
    expect(seen.size).toBe(1);
  });

  it('a changed slot remakes the scope, not the surfaces; a new `surface` function remounts them', async () => {
    const { ctx } = await openedContext();
    cleanups.push(() => ctx.dispose());
    const { n, surface } = countingSurface();
    const A = () => <div data-testid="a">a</div>;
    const B = () => <div data-testid="b">b</div>;
    const state: { Tabs: () => React.JSX.Element; surface: SeqViewerUISpec['surface'] } = { Tabs: A, surface };
    const p = parent(() => (
      <SeqViewerProvider context={ctx} ui={{ surface: state.surface, components: { ViewTabs: state.Tabs } }}>
        <SeqViewerLayout />
      </SeqViewerProvider>
    ));
    await p.render();
    expect(document.querySelector('[data-testid=a]')).not.toBeNull();
    state.Tabs = B;
    await p.render();
    expect([document.querySelector('[data-testid=b]') !== null, n.mounted, n.disposed]).toEqual([true, 1, 0]);
    state.surface = (c, s, h) => surface?.(c, s, h) as ReturnType<NonNullable<typeof surface>>;
    await p.render();
    expect([n.mounted, n.disposed]).toEqual([2, 1]);
  });

  it("`chrome` is a UI option: it changes the parts of this UI and never writes the context's config", async () => {
    const { ctx } = await openedContext();
    cleanups.push(() => ctx.dispose());
    ctx.report('seqv.test', 'a problem');
    const { n, surface } = countingSurface();
    const state = { chrome: 'full' as 'full' | 'clean' };
    const p = parent(() => (
      <SeqViewerProvider context={ctx} ui={{ surface, chrome: state.chrome }}>
        <SeqViewerLayout />
      </SeqViewerProvider>
    ));
    await p.render();
    expect(document.querySelector('.seqv-badge')).not.toBeNull();
    state.chrome = 'clean';
    await p.render();
    expect([document.querySelector('.seqv-badge'), document.querySelector('.seqv-tabs')]).toEqual([null, null]);
    expect([
      ctx.config.get(SeqViewerConfig.UI.Chrome),
      n.mounted,
      document.querySelector('.seqv-player') !== null,
    ]).toEqual(['full', 1, true]);
  });

  it('two providers on one context keep a surfaces map each', async () => {
    const { ctx } = await openedContext();
    cleanups.push(() => ctx.dispose());
    const a = countingSurface();
    const b = countingSurface();
    const p = parent(() => (
      <>
        <SeqViewerProvider context={ctx} ui={{ surface: a.surface }}>
          <SeqViewerLayout />
        </SeqViewerProvider>
        <SeqViewerProvider context={ctx} ui={{ surface: b.surface }}>
          <SeqViewerLayout />
        </SeqViewerProvider>
      </>
    ));
    await p.render();
    await p.render();
    expect([a.n.mounted, b.n.mounted]).toEqual([1, 1]);
  });
});

describe('a slot key set to `undefined`', () => {
  it("`{ Menu: undefined }` leaves the default menu: a switch's menu opens", async () => {
    const { ctx } = await openedContext();
    cleanups.push(() => ctx.dispose());
    const p = parent(() => (
      <SeqViewerProvider context={ctx} ui={{ surface: stubSurface, components: { Menu: undefined, Hint: undefined } }}>
        <SeqViewerLayout />
      </SeqViewerProvider>
    ));
    await p.render();
    const trigger = document.querySelector('[data-seqv-key="switch:section.style"]') as HTMLButtonElement;
    expect(trigger).not.toBeNull();
    await act(async () => trigger.click());
    expect(document.querySelectorAll('[role=menuitemradio]').length).toBeGreaterThan(2);
    expect(ctx.issues.current.value.map((i) => i.code)).not.toContain('seqv.ui.slot');
  });

  it('a value that is not a part falls back to the default and is reported once', async () => {
    const { ctx } = await openedContext();
    cleanups.push(() => ctx.dispose());
    const bad = { Legend: 'nope' } as unknown as SeqViewerUISpec['components'];
    const p = parent(() => (
      <SeqViewerProvider context={ctx} ui={{ surface: stubSurface, components: bad }}>
        <SeqViewerLayout />
      </SeqViewerProvider>
    ));
    await p.render();
    await p.render();
    expect(document.querySelector('.seqv-legend')).not.toBeNull();
    expect(ctx.issues.current.value.filter((i) => i.code === 'seqv.ui.slot').length).toBe(1);
  });
});

describe('SeqViewer', () => {
  const ports = () => ({ text: fixtureMeasurer() });

  it('makes a context per mount from the props of that mount, tells onReady, and disposes it on unmount', async () => {
    const made: SeqViewerContext[] = [];
    const p = parent(() => (
      <SeqViewer
        spec={headlessSpec()}
        ports={ports()}
        document={uiDoc()}
        ui={{ surface: stubSurface }}
        onReady={(c) => made.push(c)}
      />
    ));
    await p.render();
    await until(() => document.querySelector('.seqv-surface'));
    expect(made.length).toBe(1);
    expect(made[0].disposed).toBe(false);
    expect(made[0].registry.get('actions', 'view.exportSvg')).toBeDefined();
    expect(document.querySelector('.seqv-export button')?.textContent).toContain('Export');
    await act(async () => p.root.unmount());
    expect(made[0].disposed).toBe(true);
    cleanups.splice(0, 0);
  });

  it("StrictMode's second mount gets its own context, and the first is disposed", async () => {
    const made: SeqViewerContext[] = [];
    const p = parent(
      () => (
        <SeqViewer
          spec={headlessSpec()}
          ports={ports()}
          document={uiDoc()}
          ui={{ surface: stubSurface }}
          onReady={(c) => made.push(c)}
        />
      ),
      true,
    );
    await p.render();
    await until(() => document.querySelector('.seqv-surface'));
    expect(made.length).toBe(2);
    expect(made.map((c) => c.disposed)).toEqual([true, false]);
    expect(document.querySelectorAll('.seqv-root').length).toBe(1);
  });

  it('a new `document` object goes in alone, a new `data` object alone, and neither applies the other again', async () => {
    let ctx: SeqViewerContext | null = null;
    const spies = { setDocument: vi.fn(), setData: vi.fn() };
    const doc = uiDoc();
    const data = { inputs: { molecules: [] } };
    const props = { document: doc, data };
    const p = parent(() => (
      <SeqViewer
        spec={headlessSpec()}
        ports={ports()}
        document={props.document}
        data={props.data}
        ui={{ surface: stubSurface }}
        onReady={(c) => {
          ctx = c;
          // Counted, and passed on: `setData` over an empty set would replace the document.
          spies.setDocument = vi.spyOn(c, 'setDocument');
          spies.setData = vi.spyOn(c, 'setData').mockResolvedValue({ ok: true, issues: [] } as never);
        }}
      />
    ));
    await p.render();
    await until(() => document.querySelector('.seqv-surface'));
    const c = ctx as unknown as SeqViewerContext;
    const calls = () => [spies.setDocument.mock.calls.length, spies.setData.mock.calls.length];
    expect(calls()).toEqual([1, 1]);
    // An edit the reader made to the document, in the context.
    await act(async () => c.actions.run('track.hide', 't.dom'));
    await settle();
    expect(document.querySelector('.seqv-chips')?.textContent).toBe('Show Domains');
    // The parent renders again with the same objects, with a new inline `ui`: nothing is applied.
    await p.render();
    await p.render();
    expect(calls()).toEqual([1, 1]);
    // A new `data` object: `setData` alone; the document is not applied again, so the edit stays.
    props.data = { inputs: { molecules: [] } };
    await p.render();
    expect(calls()).toEqual([1, 2]);
    expect(document.querySelector('.seqv-chips')?.textContent).toBe('Show Domains');
    // A new `document` object: `setDocument` alone.
    props.document = uiDoc();
    await p.render();
    expect(calls()).toEqual([2, 2]);
    expect(document.querySelector('.seqv-chips')).toBeNull();
  });

  it('`spec` and `ports` are read at mount: later ones change nothing, and the surfaces stay with an inline `ui`', async () => {
    const made: SeqViewerContext[] = [];
    const { n, surface } = countingSurface();
    const p = parent(() => (
      <SeqViewer
        spec={headlessSpec()}
        ports={ports()}
        document={uiDoc()}
        ui={{ surface, components: {} }}
        onReady={(c) => made.push(c)}
      />
    ));
    await p.render();
    await until(() => n.mounted > 0);
    // `document` is a new object each render here, so each one is applied: the context and the surfaces stay.
    await p.render();
    await p.render();
    expect([made.length, made[0].disposed, n.mounted, n.disposed]).toEqual([1, false, 1, 0]);
  });
});
