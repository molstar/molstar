/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { readFileSync } from 'node:fs';
import { observableState } from '@molstar/prism-ui/reactive/observable-state';
import { createAgentPolicy } from '@molstar/prism-viewer/agent/policies';
import { DefaultPrismViewerSpec } from '@molstar/prism-viewer/default-spec';
import { createPrismViewerHeadless } from '@molstar/prism-viewer/headless';
import type { LazyViewerExtension, ViewerExtension } from '@molstar/prism-viewer/spec';
import type { PrismViewer } from '@molstar/prism-viewer/viewer';
import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { type ZenodoRecordInfo, ZenodoExtension } from '../prism-viewer.js';
import { LazyZenodo, ZenodoCommands, ZenodoOperations, zenodoRecordId } from '../prism-viewer-stub.js';

// Zenodo in the Prism viewer on a headless viewer (parity.md §6). No network: `fetch` (what the plugin's download uses
// in Node) answers from a record written here in Zenodo's API shape, and the file is 1CBS from data/examples.

const Id = '1234567';
const Api = `https://zenodo.org/api/records/${Id}`;
const cif = readFileSync(new URL('../../../../../data/examples/1cbs_updated.cif', import.meta.url));

const file = (key: string, size: number) => ({
  id: `0f1e2d3c-${key}`,
  key,
  size,
  checksum: 'md5:5d41402abc4b2a76b9719d911017c592',
  links: { self: `${Api}/files/${key}/content` },
});

const Record = {
  created: '2024-03-01T10:00:00.000000+00:00',
  id: Number(Id),
  doi: `10.5281/zenodo.${Id}`,
  conceptrecid: '1234566',
  metadata: { title: 'Cellular retinoic acid binding protein: structure and trajectory', doi: `10.5281/zenodo.${Id}` },
  files: [file('1cbs_updated.cif', cif.length), file('md.xtc', 4_200_000), file('README.txt', 812)],
  links: { self: Api },
};

const fetchMock = vi.fn(async (input: string | URL | Request) => {
  const url = String(input instanceof Request ? input.url : input);
  if (url === Api) return new Response(JSON.stringify(Record), { headers: { 'content-type': 'application/json' } });
  if (url === `${Api}/files/1cbs_updated.cif/content`) return new Response(cif);
  return new Response('{"status": 404, "message": "PID does not exist."}', { status: 404 });
});

beforeEach(() => {
  fetchMock.mockClear();
  vi.stubGlobal('fetch', fetchMock);
});

const viewers: PrismViewer[] = [];
afterEach(async () => {
  for (const v of viewers.splice(0)) await v.dispose();
  vi.unstubAllGlobals();
});

async function viewer(ext: ViewerExtension | LazyViewerExtension) {
  const base = DefaultPrismViewerSpec();
  const v = await createPrismViewerHeadless({ spec: { ...base, extensions: [...(base.extensions ?? []), ext] } });
  viewers.push(v);
  return v;
}

const structures = (v: PrismViewer) => v.plugin.managers.structure.hierarchy.current.structures.length;

describe('the zenodo extension in the Prism viewer', () => {
  it('reads a record id from an id, a DOI or a link', () => {
    expect(zenodoRecordId(' 1234567 ')).toBe(Id);
    expect(zenodoRecordId('10.5281/zenodo.1234567')).toBe(Id);
    expect(zenodoRecordId('https://doi.org/10.5281/zenodo.1234567')).toBe(Id);
    expect(zenodoRecordId('https://zenodo.org/records/1234567')).toBe(Id);
    expect(zenodoRecordId('https://zenodo.org/record/1234567#files')).toBe(Id);
    expect(zenodoRecordId('1cbs')).toBeUndefined();
  });

  it('is a stub until first used, and a call loads it and runs', async () => {
    const v = await viewer(LazyZenodo);
    expect(v.extensions.status('zenodo')?.value).toBe('stub');
    for (const op of Object.values(ZenodoOperations)) expect(v.api.describe(op.id)?.definition.title).toBe(op.title);
    for (const id of Object.values(ZenodoCommands)) expect(v.managers.commands.get(id)).toBeDefined();
    const r = await v.api.run('zenodo.list', { record: `10.5281/zenodo.${Id}` });
    expect(r.ok, r.ok ? '' : r.message).toBe(true);
    expect(v.extensions.status('zenodo')?.value).toBe('loaded');
    // The loaded operations are the stub's, flags and all.
    for (const op of LazyZenodo.stub.operations ?? []) {
      const def = v.api.describe(op.id)!.definition;
      expect([def.title, def.read, def.approval, def.modules]).toEqual([op.title, op.read, op.approval, op.modules]);
    }
    expect((await LazyZenodo.load()).skills?.map((s) => s.id)).toEqual(LazyZenodo.stub.skills?.map((s) => s.id));
  });

  it("lists a record's files, saying which the viewer reads", async () => {
    const v = await viewer(ZenodoExtension);
    const r = await v.api.zenodo!.list({ record: `https://zenodo.org/records/${Id}` });
    expect(r.ok, r.ok ? '' : r.message).toBe(true);
    const value = (r.ok ? r.value : undefined) as ZenodoRecordInfo;
    expect(value.title).toMatch(/retinoic acid/);
    expect(value.files.map((f) => [f.name, f.kind, f.format])).toEqual([
      ['1cbs_updated.cif', 'structure', 'mmcif'],
      ['md.xtc', 'coordinates', 'xtc'],
      ['README.txt', undefined, undefined],
    ]);
    expect(r.ok && r.did).toMatch(/3 files, 2 the viewer reads/);
    expect(fetchMock).toHaveBeenCalledTimes(1);
    // A read: no history step.
    expect(v.managers.history.steps.value.length).toBe(0);
  });

  it('refuses what is not a record, and says when Zenodo has none', async () => {
    const v = await viewer(ZenodoExtension);
    const bad = await v.api.zenodo!.list({ record: 'not a record' });
    expect(bad.ok ? 'ok' : bad.code).toBe('invalid');
    expect(fetchMock).not.toHaveBeenCalled();
    const missing = await v.api.zenodo!.list({ record: '999' });
    expect(missing.ok ? 'ok' : missing.code).toBe('not-found');
  });

  it('loads a file as one step, which Undo takes back; again is already so', async () => {
    const v = await viewer(ZenodoExtension);
    const steps = v.managers.history.steps.value.length;
    const r = await v.api.zenodo!.load({ record: Id, file: '1cbs_updated.cif' });
    expect(r.ok, r.ok ? '' : r.message).toBe(true);
    expect(r.ok && r.did).toBe(`Loaded 1cbs_updated.cif from Zenodo record ${Id}.`);
    expect(structures(v)).toBe(1);
    expect(v.managers.history.steps.value.length).toBe(steps + 1);
    // The record was read once, then the file.
    expect(fetchMock.mock.calls.map(([u]) => String(u))).toEqual([Api, `${Api}/files/1cbs_updated.cif/content`]);
    const again = await v.api.zenodo!.load({ record: Id, file: '1cbs_updated.cif' });
    expect(again.ok && again.unchanged).toBe(true);
    expect(structures(v)).toBe(1);
    await v.api.history.undo();
    expect(structures(v)).toBe(0);
  });

  it('refuses a file the record does not have, one the viewer does not read, and coordinates alone', async () => {
    const v = await viewer(ZenodoExtension);
    const code = async (values: Record<string, string>) => {
      const r = await v.api.zenodo!.load({ record: Id, ...values });
      return r.ok ? 'ok' : r.code;
    };
    expect(await code({ file: 'nope.pdb' })).toBe('not-found');
    expect(await code({ file: 'README.txt' })).toBe('invalid');
    expect(await code({ file: 'md.xtc' })).toBe('invalid');
    expect(await code({ file: '' })).toBe('invalid');
    expect(await code({ file: 'md.xtc', coordinates: '1cbs_updated.cif' })).toBe('invalid');
    expect(structures(v)).toBe(0);
    expect(v.managers.history.steps.value.length).toBe(0);
  });

  it('reaches only zenodo.org, a known service: the agent goes without asking, and asks for any other address', async () => {
    const v = await viewer(ZenodoExtension);
    const list = v.api.describe('zenodo.list')!.definition;
    const load = v.api.describe('zenodo.load')!.definition;
    expect(list.addresses!({ record: `10.5281/zenodo.${Id}` })).toEqual([Api]);
    expect(load.addresses!({ record: Id, file: 'a b.cif', coordinates: '' })).toEqual([
      Api,
      `${Api}/files/a%20b.cif/content`,
    ]);
    const policy = createAgentPolicy({ checkpointing: observableState(true).asReadonly });
    expect(policy.approves(list, { record: Id })).toBe(true);
    expect(policy.approves(load, { record: Id, file: '1cbs_updated.cif', coordinates: '' })).toBe(true);
    // The stub says the same before the extension loads (viewer-spec.md §3).
    expect(LazyZenodo.stub.operations?.[0].addresses?.({ record: Id })).toEqual([Api]);
  });

  it('brings its commands and its skill', async () => {
    const v = await viewer(ZenodoExtension);
    for (const id of Object.values(ZenodoCommands)) expect(v.managers.commands.get(id)).toBeDefined();
    expect(ZenodoExtension.skills?.[0].tools.map((t) => t.operation)).toEqual(['zenodo.list', 'zenodo.load']);
  });
});
