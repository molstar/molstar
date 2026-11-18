/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The headless host's asset policy and the context's bounds (_spec/README.md, Context, The asset policy; the
// security bounds): files under a root and the hosts the caller lists, no network by default; the spec's bounds a
// composition sets.

import { Task } from '@molstar/core/task';
import type { AssetManager } from '@molstar/core/util/assets';
import type { Node } from '@molstar/core/repr/tree';
import { SeqLanguageVersion } from '@molstar/seqx-core/lang';
import { describe, expect, it, vi } from 'vitest';
import {
  AssetRefused,
  checkLocation,
  createSeqContext,
  headlessAssetPolicy,
  HeadlessAssetPolicy,
  policyAssets,
  rootUrl,
  SeqCoreConfig,
} from '../index.js';

const signal = new AbortController().signal;
const refused = (r: ReturnType<typeof checkLocation>) => 'refused' in r;

describe('the headless policy', () => {
  it('reads nothing without a root and hosts', () => {
    const p = headlessAssetPolicy({});
    expect(p.schemes).toEqual([]);
    for (const uri of ['https://x.org/a', 'file:///etc/passwd', 'a.seq', 'http://localhost/a', '//x.org/a'])
      expect(refused(checkLocation(uri, p))).toBe(true);
    expect(refused(checkLocation('file:///etc/passwd', HeadlessAssetPolicy))).toBe(true);
  });

  it('reads files under its root, relative ones too, and nothing above it', () => {
    const p = headlessAssetPolicy({ root: '/srv/specs' });
    expect(p.root).toBe('file:///srv/specs/');
    expect(checkLocation('file:///srv/specs/a.seq', p)).toEqual({ url: 'file:///srv/specs/a.seq' });
    expect(checkLocation('data/a.seq', p)).toEqual({ url: 'file:///srv/specs/data/a.seq' });
    expect(checkLocation('file:///srv/specs/sub/../a.seq', p)).toEqual({ url: 'file:///srv/specs/a.seq' });
    for (const uri of [
      'file:///etc/passwd',
      'file:///srv/specs/../secret',
      'file:///srv/specs/%2e%2e/secret',
      'file:///srv/specs/%2E%2E/secret',
      'file:///srv/specs/a%2f..%2fb',
      'file:///srv/specs/a%5c..%5cb',
      'file:///srv/specs',
      'file:///srv/specs-other/a',
      'file://remote/srv/specs/a',
      '../secret',
      '/etc/passwd',
      '//remote/srv/specs/a',
      '\\\\remote\\share',
      'https://x.org/a',
      'http://x.org/a',
    ])
      expect(refused(checkLocation(uri, p)), uri).toBe(true);
  });

  it('reads the listed hosts over https only, and no network by default', () => {
    const p = headlessAssetPolicy({ hosts: ['data.example.org'] });
    expect(checkLocation('https://data.example.org/a.seq', p)).toEqual({ url: 'https://data.example.org/a.seq' });
    for (const uri of [
      'https://other.example.org/a',
      'https://data.example.org.evil.org/a',
      'https://data.example.org@evil.org/a',
      'http://data.example.org/a',
      'file:///srv/specs/a',
      'a.seq',
    ])
      expect(refused(checkLocation(uri, p)), uri).toBe(true);
    const both = headlessAssetPolicy({ root: '/srv/specs', hosts: ['data.example.org'] });
    expect(refused(checkLocation('https://data.example.org/a', both))).toBe(false);
    expect(refused(checkLocation('file:///srv/specs/a', both))).toBe(false);
    expect(refused(checkLocation('https://internal.local/a', both))).toBe(true);
  });

  it('never sends credentials or follows redirects, and bounds a read as the browser’s does', () => {
    const p = headlessAssetPolicy({ root: '/srv/specs' });
    expect(p.credentials).toBeUndefined();
    expect(p.redirect).toBeUndefined();
    expect(p.maxBytes).toBe(256 * 1024 * 1024);
    expect(p.timeoutMs).toBe(30_000);
    expect(headlessAssetPolicy({ root: '/x' }, { maxBytes: 10 }).maxBytes).toBe(10);
  });

  it('takes a root as a path, a Windows path or a file URL, and nothing else', () => {
    expect(rootUrl('/srv/a b/')).toBe('file:///srv/a%20b/');
    expect(rootUrl('/srv/#1/')).toBe('file:///srv/%231/');
    expect(rootUrl('C:\\data\\specs')).toBe('file:///C:/data/specs/');
    expect(rootUrl('file:///srv/x')).toBe('file:///srv/x/');
    expect(() => rootUrl('specs')).toThrow();
    expect(() => rootUrl('file://remote/x')).toThrow();
    expect(() => rootUrl('https://x.org/')).toThrow();
  });

  it('is refused before the manager is asked: a read outside it is an AssetRefused', async () => {
    const urls: string[] = [];
    const manager = {
      resolve: (a: { url: string }) => {
        urls.push(a.url);
        return Task.create('read', async () => ({ data: 'x', dispose: () => {} }));
      },
      dispose: vi.fn(),
    } as unknown as AssetManager;
    const reader = policyAssets(manager, headlessAssetPolicy({ root: '/srv/specs', hosts: ['data.example.org'] }));
    await expect(reader.readText('https://169.254.169.254/latest/meta-data', signal)).rejects.toBeInstanceOf(
      AssetRefused,
    );
    await expect(reader.readText('file:///etc/passwd', signal)).rejects.toBeInstanceOf(AssetRefused);
    await expect(reader.readText('../../etc/passwd', signal)).rejects.toBeInstanceOf(AssetRefused);
    expect(urls).toEqual([]);
    await expect(reader.readText('a.seq', signal)).resolves.toBe('x');
    expect(urls).toEqual(['file:///srv/specs/a.seq']);
  });
});

describe('the context’s bounds', () => {
  const nest = (n: number): Node => {
    let node: Node = { kind: 'ext:note', id: 'c' };
    for (let i = 1; i < n; i++) node = { kind: 'ext:note', id: 'c', children: [node] };
    return { kind: 'spec', params: { version: SeqLanguageVersion }, children: [node] };
  };

  it('are the defaults, and the composition’s item over them', () => {
    const plain = createSeqContext();
    expect(plain.bounds.maxDepth).toBe(256);
    // the chain is no valid spec: the description's issues are there, the bound's only past it
    const bound = (c: ReturnType<typeof createSeqContext>, n: number) =>
      c.check(nest(n)).map((i) => i.code).filter((code) => code.startsWith('seq.lang.bound'));
    expect(bound(plain, 256)).toEqual([]);
    expect(bound(plain, 257)).toEqual(['seq.lang.bound.depth']);
    const tight = createSeqContext({ config: [[SeqCoreConfig.bounds, { maxDepth: 3 }]] });
    expect(tight.bounds).toMatchObject({ maxDepth: 3, maxBytes: 64 * 1024 * 1024 });
    expect(bound(tight, 2)).toEqual([]);
    expect(bound(tight, 3)).toEqual([]);
    expect(tight.check(nest(4)).map((i) => i.code)).toEqual(['seq.lang.bound.depth']);
  });
});
