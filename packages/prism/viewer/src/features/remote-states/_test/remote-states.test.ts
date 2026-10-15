/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { afterEach, describe, expect, it, vi } from 'vitest';
import { exampleFile, expectOk, headlessViewer, until } from '../../../_test/harness.js';
import type { PrismViewer } from '../../../viewer.js';
import { LazyRemoteStates } from '../../lazy/index.js';
import { LoadersFeature } from '../../loaders/index.js';
import { RemoteStatesFeature } from '../index.js';
import { createRemoteStatesModel } from '../model.js';

// The `remote-states` feature: `session.remote.*` over a stubbed `fetch` that plays the state server (nothing reaches
// the network; the plugin downloads through `fetch` in Node).

const SERVER = 'https://state.test';
const LIST = [
  { id: 'b', name: 'Second', description: 'later', timestamp: 2000, isSticky: false },
  { id: 'a', name: 'First', description: '', timestamp: 1000, isSticky: false },
  { id: 'p', name: 'Pinned', description: 'always', timestamp: 3000, isSticky: true },
  { id: 'c', name: 'First', description: 'a twin', timestamp: 4000, isSticky: false },
];

interface Served {
  list?: unknown;
  /** The snapshot of `get/a`. */
  snapshot?: unknown;
}

/** A state server: answers by path, else 404. Returns the stub, to read the requests. */
function serve(served: Served, status: Record<string, number> = {}) {
  const stub = vi.fn(async (input: string | URL | Request, init?: RequestInit) => {
    const url = String(typeof input === 'object' && 'url' in input ? input.url : input);
    const path = url.slice(SERVER.length + 1).split('?')[0];
    if (status[path]) return new Response('failed', { status: status[path] });
    if (path === 'list' && served.list) return new Response(JSON.stringify(served.list));
    if (path === 'get/a' && served.snapshot) return new Response(JSON.stringify({ data: served.snapshot }));
    if (path === 'set' && init?.method === 'POST') return new Response('ok');
    if (path.startsWith('remove/')) return new Response('ok');
    if (url.includes('1cbs.bcif'))
      return new Response(new Uint8Array(await exampleFile('1cbs_full.bcif').arrayBuffer()));
    return new Response('not found', { status: 404 });
  });
  vi.stubGlobal('fetch', stub);
  return stub;
}
const requests = (stub: ReturnType<typeof serve>) => stub.mock.calls.map((c) => String(c[0]));

afterEach(() => {
  vi.unstubAllGlobals();
});

async function viewer() {
  const v = await headlessViewer({
    extensions: [LoadersFeature, RemoteStatesFeature],
    // Behind their option (viewer-spec.md §3.1).
    spec: { options: { data: { stateServer: SERVER }, ui: { showRemoteStates: true } } },
  });
  // No canvas here: states leave the camera out.
  v.plugin.state.setSnapshotParams({ camera: false });
  return v;
}
const steps = (v: PrismViewer) => v.managers.history.steps.value.length;
const structures = (v: PrismViewer) => v.plugin.managers.structure.hierarchy.current.structures;

/** A saved state of one loaded structure, as the server holds it. */
async function savedState() {
  serve({});
  const source = await viewer();
  expectOk(await source.api.load.entry!({ source: 'pdb', id: '1cbs' }));
  return JSON.parse(await (await source.plugin.managers.snapshot.serialize({ type: 'molj' })).text());
}

describe('remote states', () => {
  it('adds the session.remote operations and the commands', async () => {
    const v = await viewer();
    for (const [member, id] of [
      ['remoteList', 'session.remote.list'],
      ['remoteOpen', 'session.remote.open'],
      ['remoteUpload', 'session.remote.upload'],
      ['remoteRemove', 'session.remote.remove'],
    ] as const) {
      expect(v.api.session[member], member).toBeDefined();
      expect(v.managers.commands.get(id), id).toBeDefined();
    }
    expect(v.managers.commands.get('remote-states.open')).toBeDefined();
    expect(v.api.describe('session.remote.list')?.definition).toMatchObject({
      read: true,
      modules: ['network'],
      approval: 'network',
    });
    expect(v.api.describe('session.remote.open')?.definition).toMatchObject({
      approval: 'destructive',
      modules: ['network'],
    });
    expect(v.api.describe('session.remote.upload')?.definition).toMatchObject({
      approval: 'upload',
      history: false,
      undoable: false,
    });
    expect(v.api.describe('session.remote.remove')?.definition).toMatchObject({
      approval: 'destructive',
      undoable: false,
    });
  });

  it('lists the entries: the pinned first, then by time, with their addresses', async () => {
    const stub = serve({ list: LIST });
    const v = await viewer();
    const r = await v.api.session.remoteList!();
    expectOk(r);
    expect(requests(stub)).toEqual([`${SERVER}/list`]);
    expect(r.value!.server).toBe(SERVER);
    expect(r.value!.entries.map((e) => e.id)).toEqual(['p', 'a', 'b', 'c']);
    expect(r.value!.entries[1]).toMatchObject({
      name: 'First',
      timestamp: 1000,
      isSticky: false,
      url: `${SERVER}/get/a`,
      removeUrl: `${SERVER}/remove/a`,
    });
    expect(steps(v)).toBe(0);
  });

  it('lists another server when given one', async () => {
    const stub = vi.fn(async (_url: string) => new Response('[]'));
    vi.stubGlobal('fetch', stub);
    const v = await viewer();
    const r = await v.api.session.remoteList!({ server: 'https://other.test/states/' });
    expectOk(r);
    expect(r.value).toEqual({ server: 'https://other.test/states', entries: [] });
    expect(String(stub.mock.calls[0][0])).toBe('https://other.test/states/list');
  });

  it('fails with the server in the hint when it is down or not a state server', async () => {
    serve({}, {});
    const v = await viewer();
    const r = await v.api.session.remoteList!();
    expect(r).toMatchObject({ ok: false, code: 'failed' });
    expect(!r.ok && r.hint).toContain(SERVER);

    vi.stubGlobal(
      'fetch',
      vi.fn(async () => Promise.reject(new TypeError('fetch failed'))),
    );
    const down = await v.api.session.remoteList!();
    expect(down).toMatchObject({ ok: false, code: 'failed' });
    expect(!down.ok && down.hint).toContain(SERVER);

    serve({ list: { not: 'a list' } });
    expect(await v.api.session.remoteList!()).toMatchObject({ ok: false, code: 'failed' });
  });

  it('opens an entry by id, replacing the scene in one step that Undo restores', async () => {
    const snapshot = await savedState();
    const stub = serve({ list: LIST, snapshot });
    const v = await viewer();
    expect(structures(v)).toHaveLength(0);
    const r = await v.api.session.remoteOpen!({ entry: 'a' });
    expectOk(r);
    expect(requests(stub)).toEqual([`${SERVER}/list`, `${SERVER}/get/a`, expect.stringContaining('1cbs.bcif')]);
    await until(() => structures(v).length === 1);
    expect(steps(v)).toBe(1);
    expectOk(await v.api.history.undo());
    expect(structures(v)).toHaveLength(0);
  });

  it('opens an entry by its name, and asks which when several share it', async () => {
    const snapshot = await savedState();
    serve({ list: LIST, snapshot });
    const v = await viewer();
    expect(await v.api.session.remoteOpen!({ entry: 'first' })).toMatchObject({ ok: false, code: 'ambiguous' });
    expect(await v.api.session.remoteOpen!({ entry: 'nope' })).toMatchObject({ ok: false, code: 'not-found' });
    expect(await v.api.session.remoteOpen!({})).toMatchObject({ ok: false, code: 'invalid' });
    expect(steps(v)).toBe(0);
    // `Second` is unique but the stub serves only `get/a`: its name resolves, the fetch then fails.
    const r = await v.api.session.remoteOpen!({ entry: 'second' });
    expect(r).toMatchObject({ ok: false, code: 'failed' });
    expect(!r.ok && r.hint).toContain(SERVER);
    expect(structures(v)).toHaveLength(0);
  });

  it('fails on an entry that is not a state, leaving the scene', async () => {
    serve({ list: LIST, snapshot: { not: 'a state' } });
    const v = await viewer();
    expect(await v.api.session.remoteOpen!({ entry: 'a' })).toMatchObject({ ok: false, code: 'failed' });
    expect(structures(v)).toHaveLength(0);
  });

  it('uploads the current state with its name and description in the query', async () => {
    const stub = serve({ list: LIST });
    const v = await viewer();
    // Nothing loaded: nothing to upload.
    expect(await v.api.session.remoteUpload!({ name: 'Mine' })).toMatchObject({ ok: false, code: 'blocked' });
    expectOk(await v.api.load.entry!({ source: 'pdb', id: '1cbs' }));
    const before = steps(v);
    const r = await v.api.session.remoteUpload!({ name: 'My scene', description: 'a & b', playOnLoad: true });
    expectOk(r);
    const post = stub.mock.calls.find((c) => c[1]?.method === 'POST')!;
    expect(String(post[0])).toBe(`${SERVER}/set?name=My%20scene&description=a%20%26%20b`);
    const body = JSON.parse(String(post[1]!.body));
    expect(body).toMatchObject({ name: 'My scene', description: 'a & b' });
    expect(body.entries.length).toBeGreaterThan(0);
    expect(steps(v)).toBe(before);
  });

  it('fails an upload the server refuses, with the server in the hint', async () => {
    serve({ list: LIST }, { set: 500 });
    const v = await viewer();
    expectOk(await v.api.load.entry!({ source: 'pdb', id: '1cbs' }));
    const r = await v.api.session.remoteUpload!({ name: 'x' });
    expect(r).toMatchObject({ ok: false, code: 'failed' });
    expect(!r.ok && r.hint).toContain(SERVER);
  });

  it('removes an entry by calling its removeUrl, except a pinned one', async () => {
    const stub = serve({ list: LIST });
    const v = await viewer();
    expectOk(await v.api.session.remoteRemove!({ entry: 'b' }));
    expect(requests(stub)).toContain(`${SERVER}/remove/b`);
    expect(await v.api.session.remoteRemove!({ entry: 'p' })).toMatchObject({ ok: false, code: 'blocked' });
    expect(requests(stub)).not.toContain(`${SERVER}/remove/p`);
    expect(await v.api.session.remoteRemove!({ entry: 'zzz' })).toMatchObject({ ok: false, code: 'not-found' });
    expect(steps(v)).toBe(0);
  });

  it('fails a removal the server refuses', async () => {
    serve({ list: LIST }, { 'remove/b': 404 });
    const v = await viewer();
    const r = await v.api.session.remoteRemove!({ entry: 'b' });
    expect(r).toMatchObject({ ok: false, code: 'failed' });
    expect(!r.ok && r.hint).toContain(SERVER);
  });
});

describe('the Remote states model', () => {
  it('lists, adopts the answering server, uploads and removes, refreshing after', async () => {
    const stub = serve({ list: LIST });
    const v = await viewer();
    expectOk(await v.api.load.entry!({ source: 'pdb', id: '1cbs' }));
    const posted: string[] = [];
    const model = createRemoteStatesModel({ api: v.api, notify: (t, k) => posted.push(`${k}:${t}`) });
    model.refresh();
    await until(() => model.states.value.length === 4);
    expect(model.serverUrl.value).toBe(SERVER);
    expect(model.states.value.map((s) => s.id)).toEqual(['p', 'a', 'b', 'c']);
    expect(model.states.value[0]).toMatchObject({ label: 'Pinned', date: expect.stringMatching(/1970/) });

    model.upload!('Mine');
    await until(() => stub.mock.calls.some((c) => c[1]?.method === 'POST'));
    await until(() => requests(stub).filter((u) => u.endsWith('/list')).length === 2);
    model.remove!('b');
    await until(() => requests(stub).includes(`${SERVER}/remove/b`));
    await until(() => !model.loading.value);
    expect(posted.filter((p) => p.startsWith('error'))).toEqual([]);

    // A server that fails posts the error and empties the list.
    vi.stubGlobal(
      'fetch',
      vi.fn(async () => new Response('x', { status: 500 })),
    );
    model.setServerUrl('https://down.test');
    await until(() => posted.some((p) => p.startsWith('error') && p.includes('https://down.test')));
    expect(model.states.value).toEqual([]);
  });
});

describe('the lazy stub', () => {
  it('declares what the feature brings, and the feature loads', async () => {
    const feature = await LazyRemoteStates.load();
    expect(feature.builtIn).toBe(LazyRemoteStates.builtIn);
    const v = await headlessViewer({
      extensions: [LazyRemoteStates],
      spec: { options: { ui: { showRemoteStates: true } } },
    });
    expect(v.extensions.status(LazyRemoteStates.id)?.value).toBe('stub');
    expect(await v.extensions.load(LazyRemoteStates.id)).toBe(true);
    expect(v.extensions.status(LazyRemoteStates.id)?.value).toBe('loaded');
    for (const op of LazyRemoteStates.stub.operations ?? []) {
      expect(v.api.describe(op.id)?.definition.title).toBe(op.title);
    }
    for (const c of LazyRemoteStates.stub.commands ?? []) expect(v.managers.commands.get(c.id)?.title).toBe(c.title);
  });
});
