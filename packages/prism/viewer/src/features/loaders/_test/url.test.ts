/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { afterEach, describe, expect, it, vi } from 'vitest';
import { exampleFile, headlessViewer, recordRuns, until } from '../../../_test/harness.js';
import { LoadersFeature, LoadUrlParameters, loadsFromUrl, parseSearch, runUrlLoads } from '../index.js';

// The URL's load parameters (app.md §4): each row of the table becomes an operation call, DOM-free.

afterEach(() => {
  vi.unstubAllGlobals();
});

describe('loadsFromUrl', () => {
  it('has nothing for a URL without load parameters', () => {
    expect(loadsFromUrl('')).toEqual([]);
    expect(loadsFromUrl('?hide-controls=1&pixel-scale=2')).toEqual([]);
  });

  it.each([
    ['pdb', 'pdb', '1tqn'],
    ['pdb-ihm', 'pdb-ihm', '8zzc'],
    ['pdb-dev', 'pdb-ihm', '8zzc'],
    ['afdb', 'afdb', 'Q8W3K0'],
    ['emdb', 'emdb', 'EMD-1832'],
    ['model-archive', 'model-archive', 'ma-bak-cepc-0003'],
  ])('reads %s as load.entry of %s', (name, source, id) => {
    expect(loadsFromUrl(`?${name}=${id}`)).toEqual([{ operation: 'load.entry', values: { source, id } }]);
  });

  it('reads url with its format and binary flag; a missing format is read from the address', () => {
    expect(loadsFromUrl('?url=https://x.org/a.cif&url-format=mmcif&url-is-binary=0')).toEqual([
      { operation: 'load.url', values: { url: 'https://x.org/a.cif', format: 'mmcif', binary: false } },
    ]);
    expect(loadsFromUrl('?url=https://x.org/a.bcif&url-format=mmcif&url-is-binary=1')[0].values).toMatchObject({
      binary: true,
    });
    expect(loadsFromUrl('?url=https://x.org/a.bcif')[0].values).toMatchObject({ format: 'auto', binary: false });
  });

  it('reads structure-url with its format (default mmcif) and binary flag', () => {
    expect(loadsFromUrl('?structure-url=https://x.org/a.pdb&structure-url-format=PDB')).toEqual([
      { operation: 'load.url', values: { url: 'https://x.org/a.pdb', format: 'pdb', binary: false } },
    ]);
    expect(loadsFromUrl('?structure-url=https://x.org/a.bcif&structure-url-is-binary=1')[0].values).toEqual({
      url: 'https://x.org/a.bcif',
      format: 'mmcif',
      binary: true,
    });
  });

  it('reads mvs-url, mvs-data and mvs-format', () => {
    expect(loadsFromUrl('?mvs-url=https://x.org/a.mvsj')).toEqual([
      { operation: 'load.mvs', values: { url: 'https://x.org/a.mvsj', format: 'mvsj' } },
    ]);
    expect(loadsFromUrl('?mvs-url=https://x.org/a.mvsx&mvs-format=mvsx')[0].values).toMatchObject({ format: 'mvsx' });
    expect(loadsFromUrl('?mvs-data=%7B%22a%22%3A1%7D')[0].values).toEqual({ data: '{"a":1}', format: 'mvsj' });
  });

  it('keeps a + in mvs-data (base64), as v5 does', () => {
    expect(loadsFromUrl('?mvs-data=base64,a+b/c%3D&mvs-format=mvsx')[0].values).toMatchObject({
      data: 'base64,a+b/c=',
    });
  });

  it('passes mvs-url and mvs-data together on, for the operation to refuse', () => {
    expect(loadsFromUrl('?mvs-url=https://x.org/a.mvsj&mvs-data=%7B%7D')).toEqual([
      { operation: 'load.mvs', values: { url: 'https://x.org/a.mvsj', data: '{}', format: 'mvsj' } },
    ]);
  });

  it('reads snapshot-url (its type only when given: the operation follows the extension) and snapshot-id', () => {
    expect(loadsFromUrl('?snapshot-url=https://x.org/s.molj')).toEqual([
      { operation: 'load.snapshot', values: { url: 'https://x.org/s.molj' } },
    ]);
    expect(loadsFromUrl('?snapshot-url=https://x.org/s.molx&snapshot-url-type=MOLX')[0].values).toMatchObject({
      type: 'molx',
    });
    expect(loadsFromUrl('?snapshot-id=abc123')).toEqual([{ operation: 'load.snapshot', values: { id: 'abc123' } }]);
  });

  it('runs them in v5 order, takes the first of a repeated name, ignores empty values and the rest', () => {
    const ops = loadsFromUrl(
      '?pdb=1cbs&snapshot-id=s1&pdb=2abc&afdb=&url=https://x.org/a.cif&hide-controls=1&emdb=EMD-1',
    );
    expect(ops.map((o) => `${o.operation} ${o.values.source ?? o.values.id ?? o.values.url}`)).toEqual([
      'load.snapshot s1',
      'load.entry pdb',
      'load.entry emdb',
      'load.url https://x.org/a.cif',
    ]);
    expect(ops[1].values.id).toBe('1cbs');
  });

  it('knows its parameter names, case-insensitively', () => {
    expect(LoadUrlParameters).toContain('model-archive');
    expect([...parseSearch('?PDB=1cbs&x').keys()]).toEqual(['pdb', 'x']);
    expect(loadsFromUrl('?PDB=1cbs')).toHaveLength(1);
  });
});

describe('runUrlLoads', () => {
  it('runs the loads in order with the url origin', async () => {
    const bcif = new Uint8Array(await exampleFile('1cbs_full.bcif').arrayBuffer());
    const seen: string[] = [];
    vi.stubGlobal(
      'fetch',
      vi.fn(async (input: string | URL | Request) => {
        const url = String(input);
        seen.push(url);
        return url.includes('1cbs.bcif') ? new Response(bcif, { status: 200 }) : new Response('no', { status: 404 });
      }),
    );
    const v = await headlessViewer({ extensions: [LoadersFeature] });
    const runs = recordRuns(v);
    const results = await runUrlLoads(v.api, '?pdb=1cbs&url=https://x.test/none.cif&url-format=mmcif');
    await until(() => runs.length === 2);
    expect(results.map((r) => [r.operation, r.result.ok])).toEqual([
      ['load.entry', true],
      ['load.url', false],
    ]);
    expect(runs.map((r) => [r.operation, r.origin.kind])).toEqual([
      ['load.entry', 'url'],
      ['load.url', 'url'],
    ]);
    expect(seen[0]).toContain('1cbs.bcif');
    expect(v.plugin.managers.structure.hierarchy.current.structures).toHaveLength(1);
    // A failed load does not stop the rest, and is told to the person (the status line).
    expect(v.managers.notifications.current.value?.kind).toBe('warning');
  });

  it('does nothing without load parameters', async () => {
    const v = await headlessViewer({ extensions: [LoadersFeature] });
    expect(await runUrlLoads(v.api, '?hide-controls=1')).toEqual([]);
  });

  it('refuses mvs-url with mvs-data together, loading neither', async () => {
    const v = await headlessViewer({ extensions: [LoadersFeature] });
    const [r] = await runUrlLoads(v.api, '?mvs-url=https://x.test/a.mvsj&mvs-data=%7B%7D');
    expect(r.result).toMatchObject({ ok: false, code: 'invalid' });
  });
});
