/** @vitest-environment jsdom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { afterEach, describe, expect, it, vi } from 'vitest';
import {
  createPrismViewerElementClass,
  definePrismViewerElement,
  type PrismViewerElement,
  specWithAttributes,
} from '../element.js';
import type { CreatePrismViewerOptions } from '../page.js';
import type { PrismViewerSpec } from '../spec.js';
import type { PrismViewer } from '../viewer.js';

// `<prism-viewer>` (app.md §3) over a stand-in factory: connecting makes the viewer with the element as its target,
// the attributes become options, `pdb` loads, leaving the page disposes.

const BaseSpec = { plugin: { registry: [], behaviors: [] } } as PrismViewerSpec;

function stubViewer(loads = true) {
  const run = vi.fn(async (_id: string, _values: unknown) => ({ ok: true }) as const);
  const dispose = vi.fn(async () => {});
  const viewer = {
    api: { run, describe: (id: string) => (loads && id === 'load.entry' ? {} : undefined) },
    dispose,
  } as unknown as PrismViewer;
  return { viewer, run, dispose };
}

let n = 0;
function define(create: (o: CreatePrismViewerOptions) => Promise<PrismViewer>) {
  const tag = `test-prism-viewer-${++n}`;
  customElements.define(tag, createPrismViewerElementClass({ create, defaultSpec: () => BaseSpec }));
  return (attributes: Record<string, string> = {}) => {
    const el = document.createElement(tag) as PrismViewerElement;
    for (const [k, v] of Object.entries(attributes)) el.setAttribute(k, v);
    return el;
  };
}

const settle = () => new Promise((r) => setTimeout(r, 0));

afterEach(() => {
  document.body.innerHTML = '';
});

describe('<prism-viewer>', () => {
  it('makes its viewer when connected, the element its target with a shadow root, and loads `pdb`', async () => {
    const stub = stubViewer();
    const create = vi.fn(async (_: CreatePrismViewerOptions) => stub.viewer);
    const make = define(create);
    const el = make({ pdb: '1cbs' });
    expect(create).not.toHaveBeenCalled();
    const ready = vi.fn();
    el.addEventListener('ready', ready);
    document.body.append(el);
    expect(await el.ready).toBe(stub.viewer);
    expect(create).toHaveBeenCalledOnce();
    expect(create.mock.calls[0][0].target).toBe(el);
    expect(create.mock.calls[0][0].spec.plugin).toBe(BaseSpec.plugin);
    expect(el.shadowRoot).not.toBeNull();
    expect(el.viewer).toBe(stub.viewer);
    expect(ready).toHaveBeenCalledOnce();
    await settle();
    expect(stub.run).toHaveBeenCalledWith('load.entry', { source: 'pdb', id: '1cbs' });

    // A new entry loads; the same one does not again.
    el.setAttribute('pdb', '1tqn');
    el.setAttribute('pdb', '1tqn');
    await settle();
    expect(stub.run.mock.calls.map((c) => (c[1] as { id: string }).id)).toEqual(['1cbs', '1tqn']);
  });

  it('turns its attributes into options: theme and hide-controls', async () => {
    const create = vi.fn(async (_: CreatePrismViewerOptions) => stubViewer().viewer);
    const el = define(create)({ theme: 'dark', 'hide-controls': '' });
    document.body.append(el);
    await el.ready;
    const spec = create.mock.calls[0][0].spec;
    expect(spec.ui?.theme).toBe('dark');
    expect(spec.options?.ui?.collapsed).toBe(true);

    // An unknown theme is ignored; hide-controls="false" is off.
    expect(specWithAttributes(BaseSpec, { theme: 'neon', hideControls: 'false' })).toEqual({
      ...BaseSpec,
      options: { ui: { collapsed: false } },
    });
    expect(specWithAttributes(BaseSpec, {})).toEqual(BaseSpec);
  });

  it('takes its spec from the `spec` property, and makes the viewer again when it changes', async () => {
    const first = stubViewer();
    const second = stubViewer();
    const create = vi
      .fn(async (_: CreatePrismViewerOptions) => first.viewer)
      .mockResolvedValueOnce(first.viewer)
      .mockResolvedValueOnce(second.viewer);
    const el = define(create)();
    const own = { ...BaseSpec, overrides: {} } as PrismViewerSpec;
    el.spec = own;
    document.body.append(el);
    await el.ready;
    expect(create.mock.calls[0][0].spec).toEqual(own);
    const next = { ...own, ui: { theme: 'blue' as const } };
    el.spec = next;
    await vi.waitFor(() => expect(el.viewer).toBe(second.viewer));
    expect(first.dispose).toHaveBeenCalledOnce();
    expect(create.mock.calls[1][0].spec.ui?.theme).toBe('blue');
  });

  it('disposes its viewer when it leaves the page, and makes a new one when it comes back', async () => {
    const first = stubViewer();
    const second = stubViewer();
    const create = vi.fn().mockResolvedValueOnce(first.viewer).mockResolvedValueOnce(second.viewer);
    const el = define(create)();
    document.body.append(el);
    await el.ready;
    el.remove();
    await vi.waitFor(() => expect(first.dispose).toHaveBeenCalledOnce());
    expect(el.viewer).toBeUndefined();

    document.body.append(el);
    expect(await el.ready).toBe(second.viewer);
    expect(create).toHaveBeenCalledTimes(2);
  });

  it('keeps its viewer when moved within the page', async () => {
    const stub = stubViewer();
    const create = vi.fn(async () => stub.viewer);
    const el = define(create)();
    const a = document.createElement('div');
    const b = document.createElement('div');
    document.body.append(a, b);
    a.append(el);
    await el.ready;
    b.append(el);
    await settle();
    expect(stub.dispose).not.toHaveBeenCalled();
    expect(create).toHaveBeenCalledOnce();
  });

  it('disposes a viewer that arrives after the element left', async () => {
    const stub = stubViewer();
    let finish!: (v: PrismViewer) => void;
    const create = vi.fn(() => new Promise<PrismViewer>((r) => (finish = r)));
    const el = define(create)();
    document.body.append(el);
    el.remove();
    await settle();
    finish(stub.viewer);
    await vi.waitFor(() => expect(stub.dispose).toHaveBeenCalledOnce());
    expect(el.viewer).toBeUndefined();
  });

  it('several on a page each own their viewer', async () => {
    const viewers = [stubViewer(), stubViewer()];
    let i = 0;
    const create = vi.fn(async () => viewers[i++].viewer);
    const make = define(create);
    const a = make();
    const b = make();
    document.body.append(a, b);
    expect(await a.ready).toBe(viewers[0].viewer);
    expect(await b.ready).toBe(viewers[1].viewer);
  });

  it('is defined once per tag', () => {
    const tag = `test-prism-viewer-defined-${++n}`;
    expect(definePrismViewerElement({ tag, create: async () => stubViewer().viewer })).toBe(true);
    expect(definePrismViewerElement({ tag })).toBe(false);
    expect(customElements.get(tag)).toBeDefined();
  });
});
