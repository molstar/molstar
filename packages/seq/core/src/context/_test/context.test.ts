/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Task } from '@molstar/core/task';
import type { AssetManager } from '@molstar/core/util/assets';
import { SizeLimitExceeded } from '@molstar/core/util/assets';
import { type FormatRuntime, SeqRuntime } from '@molstar/seqx-core/data';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { DefaultRenderConfig } from '@molstar/seqx-core/render';
import { describe, expect, it, vi } from 'vitest';
import {
  AssetRefused,
  BrowserAssetPolicy,
  checkLocation,
  configItem,
  createSeqContext,
  HeadlessAssetPolicy,
  policyAssets,
  utf8Length,
  SeqCoreConfig,
} from '../index.js';

const noSignal = new AbortController().signal;

function fakeManager(text: string, urls: string[] = []) {
  const dispose = vi.fn();
  const manager = {
    resolve: (asset: { url: string }) => {
      urls.push(asset.url);
      return Task.create('read', async () => ({ data: text, dispose: () => {} }));
    },
    dispose,
  } as unknown as AssetManager;
  return { manager, dispose };
}

function tinySpec() {
  const b = createSeqSpec({ title: 'tiny' });
  const p = b.data.sequence('p', { alphabet: 'protein', letters: 'MEEPQSDP' });
  b.views.view('main', (v) => {
    v.section('protein', { of: p }).track('letters').layer('letters', { type: 'letters' });
  });
  return b.build();
}

describe('asset policy', () => {
  const base = 'https://example.org/app/page.html';

  it('reads nothing by default, relative locations under the base and https hosts by the browser policy', () => {
    expect(checkLocation('a.seq', HeadlessAssetPolicy)).toHaveProperty('refused');
    expect(checkLocation('https://x.org/a', HeadlessAssetPolicy)).toHaveProperty('refused');
    const p = { ...BrowserAssetPolicy, base };
    expect(checkLocation('data/a.seq', p)).toEqual({ url: 'https://example.org/app/data/a.seq' });
    expect(checkLocation('https://x.org/a', p)).toEqual({ url: 'https://x.org/a' });
    expect(checkLocation('http://x.org/a', p)).toHaveProperty('refused');
    expect(checkLocation('file:///etc/passwd', p)).toHaveProperty('refused');
    expect(checkLocation('javascript:alert(1)', p)).toHaveProperty('refused');
  });

  it('confines relative locations to the base folder and absolute ones to the listed hosts', () => {
    const p = { ...BrowserAssetPolicy, base, hosts: ['x.org'] };
    expect(checkLocation('../secret', p)).toEqual({ refused: 'outside the base folder' });
    expect(checkLocation('//evil.org/a', p)).toHaveProperty('refused');
    expect(checkLocation('https://evil.org/a', p)).toHaveProperty('refused');
    expect(checkLocation('../ok', { ...p, confine: false })).toEqual({ url: 'https://example.org/ok' });
  });

  it('classifies a location by the URL it resolves to: no spelling of another host or folder gets through', () => {
    const p = { ...HeadlessAssetPolicy, relative: true, base: 'https://example.org/app/page.html', confine: true };
    const spellings = [
      ' //evil.com/a',
      '\\\\evil.com/a',
      '/\\evil.com/a',
      '\t//evil.com/a',
      '\n/\t/evil.com/a',
      '/etc/passwd',
      '/app/../x',
      '%2e%2e/secret',
      '.%2e/secret',
      '..%2fsecret',
      'a/../../secret',
      '\thttp://evil.com/',
      ' https://evil.com/',
      'file:///etc/passwd',
      'data:text/plain,hi',
      'blob:https://example.org/abc',
      'javascript:alert(1)',
    ];
    for (const u of spellings) expect(checkLocation(u, p), JSON.stringify(u)).toHaveProperty('refused');
    // with the browser's schemes and a host list the same spellings reach no other host either
    const listed = { ...BrowserAssetPolicy, base: p.base, hosts: ['x.org'] };
    for (const u of [' //evil.com/a', '\\\\evil.com/a', '/\\evil.com/a', '\thttps://evil.com/', 'data:text/plain,hi'])
      expect(checkLocation(u, listed), JSON.stringify(u)).toHaveProperty('refused');
    // what is allowed stays allowed, spelled with the parser's own forgiveness
    expect(checkLocation(' d/a.seq\n', p)).toEqual({ url: 'https://example.org/app/d/a.seq' });
    expect(checkLocation('d/./b/../a.seq', p)).toEqual({ url: 'https://example.org/app/d/a.seq' });
  });

  it('a file base: relative locations stay in its folder, another host is refused', () => {
    const p = { ...HeadlessAssetPolicy, relative: true, base: 'file:///pkg/fixtures/', confine: true };
    expect(checkLocation('p53/a.seq', p)).toEqual({ url: 'file:///pkg/fixtures/p53/a.seq' });
    for (const u of ['../package.json', '//host/share/a', '/etc/passwd', 'file:///etc/passwd', '%2e%2e/x'])
      expect(checkLocation(u, p), u).toHaveProperty('refused');
  });

  it('a refused read rejects before the manager is asked; a read under the policy goes through', async () => {
    const urls: string[] = [];
    const { manager } = fakeManager('hello', urls);
    const reader = policyAssets(manager, { ...BrowserAssetPolicy, base });
    await expect(reader.readText('http://x.org/a', noSignal)).rejects.toBeInstanceOf(AssetRefused);
    expect(urls).toEqual([]);
    expect(await reader.readText('d.seq', noSignal)).toBe('hello');
    expect(urls).toEqual(['https://example.org/app/d.seq']);
  });

  it('refuses an asset over the size limit and a read that takes longer than its limit', async () => {
    const { manager } = fakeManager('x'.repeat(100));
    const small = policyAssets(manager, { ...BrowserAssetPolicy, base, maxBytes: 10 });
    await expect(small.readText('d.seq', noSignal)).rejects.toThrow('larger than the limit');
    // bytes, not characters: 4 letters of 3 bytes are over a limit of 10
    const wide = policyAssets(fakeManager('\u20ac\u20ac\u20ac\u20ac').manager, {
      ...BrowserAssetPolicy,
      base,
      maxBytes: 10,
    });
    await expect(wide.readText('d.seq', noSignal)).rejects.toThrow('larger than the limit');
    expect(utf8Length('a\u00e9\u20ac\ud83d\ude00')).toBe(1 + 2 + 3 + 4);
    vi.useFakeTimers();
    try {
      const slow = {
        resolve: () => Task.create('slow', () => new Promise<never>(() => {})),
      } as unknown as AssetManager;
      const reader = policyAssets(slow, { ...BrowserAssetPolicy, base, timeoutMs: 50 });
      const read = expect(reader.readText('d.seq', noSignal)).rejects.toThrow('took too long');
      await vi.advanceTimersByTimeAsync(60);
      await read;
    } finally {
      vi.useRealTimers();
    }
  });

  it('passes its choices to the manager on every read (the four options)', async () => {
    const seen: any[] = [];
    const manager = {
      resolve: (_: unknown, __: unknown, options: unknown) => {
        seen.push(options);
        return Task.create('read', async () => ({ data: 'hello', dispose: () => {} }));
      },
    } as unknown as AssetManager;
    const c = new AbortController();
    await policyAssets(manager, { ...BrowserAssetPolicy, base }).readText('d.seq', c.signal);
    expect(seen[0]).toMatchObject({ credentials: 'omit', redirect: 'error', maxBytes: BrowserAssetPolicy.maxBytes });
    expect(seen[0].signal).toBeInstanceOf(AbortSignal);
    // the composition's credentials, and no bound when the policy has none
    await policyAssets(manager, { ...HeadlessAssetPolicy, schemes: ['https:'], credentials: 'include' }).readText(
      'https://x.org/a',
      c.signal,
    );
    expect(seen[1].credentials).toBe('include');
    expect(seen[1].redirect).toBe('error');
    expect('maxBytes' in seen[1]).toBe(false);
  });

  it('the caller’s abort and the timeout reach the manager as its signal', async () => {
    let managerSignal: AbortSignal | undefined;
    const hang = {
      resolve: (_: unknown, __: unknown, options: { signal: AbortSignal }) => {
        managerSignal = options.signal;
        return Task.create('slow', () => new Promise<never>(() => {}));
      },
    } as unknown as AssetManager;
    const c = new AbortController();
    const read = expect(
      policyAssets(hang, { ...BrowserAssetPolicy, base }).readText('d.seq', c.signal),
    ).rejects.toThrow('stop');
    await Promise.resolve();
    expect(managerSignal!.aborted).toBe(false);
    c.abort(new Error('stop'));
    await read;
    expect(managerSignal!.aborted).toBe(true);

    vi.useFakeTimers();
    try {
      const timed = expect(
        policyAssets(hang, { ...BrowserAssetPolicy, base, timeoutMs: 50 }).readText('d.seq', noSignal),
      ).rejects.toThrow('took too long');
      await vi.advanceTimersByTimeAsync(60);
      await timed;
      expect(managerSignal!.aborted).toBe(true);
    } finally {
      vi.useRealTimers();
    }
  });

  it('a body over the limit that the manager stopped is a refused read', async () => {
    const full = {
      resolve: () =>
        Task.create('read', async () => {
          throw new SizeLimitExceeded(10);
        }),
    } as unknown as AssetManager;
    await expect(
      policyAssets(full, { ...BrowserAssetPolicy, base, maxBytes: 10 }).readText('d.seq', noSignal),
    ).rejects.toBeInstanceOf(AssetRefused);
  });

  it('a refused redirect is told from a failed download where the platform says, and follow is a policy option', async () => {
    const failing = (message: string) =>
      ({
        resolve: () =>
          Task.create('read', async () => {
            throw new TypeError('fetch failed', { cause: new Error(message) });
          }),
      }) as unknown as AssetManager;
    const refused = policyAssets(failing('unexpected redirect'), { ...BrowserAssetPolicy, base });
    await expect(refused.readText('d.seq', noSignal)).rejects.toBeInstanceOf(AssetRefused);
    await expect(refused.readText('d.seq', noSignal)).rejects.toThrow('a redirect is not followed');
    // a network failure is not a redirect
    const down = policyAssets(failing('ECONNREFUSED'), { ...BrowserAssetPolicy, base });
    await expect(down.readText('d.seq', noSignal)).rejects.toBeInstanceOf(TypeError);
    // a size error of another copy of the class is the same refusal
    const copy = {
      resolve: () =>
        Task.create('read', async () => {
          throw Object.assign(new Error('big'), { name: 'SizeLimitExceeded' });
        }),
    } as unknown as AssetManager;
    await expect(
      policyAssets(copy, { ...BrowserAssetPolicy, base }).readText('d.seq', noSignal),
    ).rejects.toBeInstanceOf(AssetRefused);
    const seen: any[] = [];
    const recording = {
      resolve: (_: unknown, __: unknown, options: unknown) => {
        seen.push(options);
        return Task.create('read', async () => ({ data: 'x', dispose: () => {} }));
      },
    } as unknown as AssetManager;
    await policyAssets(recording, { ...BrowserAssetPolicy, base, redirect: 'follow' }).readText('d.seq', noSignal);
    expect(seen[0].redirect).toBe('follow');
  });

  it('on the default manager: no redirect is followed, no credentials are sent, a long body is cut off', async () => {
    const { DefaultAssetManager } = await import('@molstar/core/util/assets');
    const requests: RequestInit[] = [];
    const original = globalThis.fetch;
    globalThis.fetch = (async (input: RequestInfo | URL, init: RequestInit = {}) => {
      requests.push(init);
      if (String(input).endsWith('/moved') && init.redirect === 'error') throw new TypeError('redirect');
      return new Response('x'.repeat(100), { headers: { 'Content-Length': '100' } });
    }) as typeof fetch;
    try {
      const reader = policyAssets(new DefaultAssetManager(), { ...BrowserAssetPolicy, base, maxBytes: 50 });
      await expect(reader.readText('moved', noSignal)).rejects.toThrow('redirect');
      await expect(reader.readText('d.seq', noSignal)).rejects.toThrow('larger than the limit');
      expect(requests.every((r) => r.credentials === 'omit' && r.redirect === 'error')).toBe(true);
    } finally {
      globalThis.fetch = original;
    }
  });

  it('an aborted read rejects with the signal’s reason', async () => {
    const c = new AbortController();
    const slow = { resolve: () => Task.create('slow', () => new Promise<never>(() => {})) } as unknown as AssetManager;
    const reader = policyAssets(slow, { ...BrowserAssetPolicy, base });
    const read = expect(reader.readText('d.seq', c.signal)).rejects.toThrow('stop');
    c.abort(new Error('stop'));
    await read;
  });
});

describe('createSeqContext', () => {
  it('renders a spec headless, with the bundled text metrics and no reads', async () => {
    const ctx = createSeqContext();
    const r = await ctx.renderSvg(tinySpec(), { size: { width: 400, height: 60 } });
    expect(r.settled).toBe(true);
    expect(r.issues).toEqual([]);
    expect(r.svg.startsWith('<svg')).toBe(true);
    await expect(ctx.assets.readText('a.seq', noSignal)).rejects.toThrow();
    ctx.retire();
  });

  it('resolves with the registered runtimes and checks a spec', () => {
    const ctx = createSeqContext({ registry: [{ kind: 'format-runtime', runtime: SeqRuntime }] });
    expect(ctx.runtimes).toEqual([SeqRuntime]);
    expect(ctx.check(tinySpec())).toEqual([]);
    const model = ctx.resolve(tinySpec());
    expect(ctx.view(model)?.sections.length).toBe(1);
    expect(ctx.view(model, '/views/none')).toBeUndefined();
  });

  it('the registry takes the same runtime again, and refuses another under one name, listing all', () => {
    const other: FormatRuntime = { ...SeqRuntime };
    const entry = { kind: 'format-runtime', runtime: SeqRuntime } as const;
    expect(createSeqContext({ registry: [entry, entry] }).runtimes.length).toBe(1);
    expect(() => createSeqContext({ registry: [entry, { kind: 'format-runtime', runtime: other }] })).toThrow(/seq/);
  });

  it('config items set the theme, the layout and the resolve bounds', async () => {
    const ctx = createSeqContext({
      config: [
        [SeqCoreConfig.theme, 'dark'],
        [SeqCoreConfig.render, { inset: 20 }],
      ],
    });
    expect(ctx.renderConfig.inset).toBe(20);
    expect(ctx.tokens).not.toBe(createSeqContext().tokens);
    expect(configItem<number>('x').key).toBe('x');
  });

  it('the render item takes nested settings field by field (`RenderConfigOverrides`)', () => {
    // a type test too (`typecheck:test`): a nested partial (`values` without `aggregate`) is the item's value type
    const over: NonNullable<(typeof SeqCoreConfig.render)['type']> = { values: { coarser: 1 }, outline: { hover: 3 } };
    const ctx = createSeqContext({ config: [[SeqCoreConfig.render, over]] });
    expect(ctx.renderConfig.values).toEqual({ ...DefaultRenderConfig.values, coarser: 1 });
    expect(ctx.renderConfig.outline.hover).toBe(3);
    expect(ctx.renderConfig.outline.cursor).toBe(DefaultRenderConfig.outline.cursor);
  });

  it('retire disposes the ports it was told to, and leaves a host’s', () => {
    const owned = fakeManager('');
    const kept = fakeManager('');
    const measure = vi.fn(() => 5);
    const dispose = vi.fn();
    const a = createSeqContext({
      ports: { assets: () => ({ manager: owned.manager }), textMetrics: () => ({ measure, dispose }) },
    });
    const b = createSeqContext({
      ports: {
        assets: () => ({ manager: kept.manager, disposeWithViewer: false }),
        textMetrics: () => ({ measure, dispose, disposeWithViewer: false }),
      },
    });
    a.retire();
    a.retire();
    b.retire();
    expect(owned.dispose).toHaveBeenCalledTimes(1);
    expect(kept.dispose).not.toHaveBeenCalled();
    expect(dispose).toHaveBeenCalledTimes(1);
  });

  it('takes a reader as the assets port (a second adapter), disposed with the context unless the host keeps it', async () => {
    const dispose = vi.fn();
    const reader = { readText: vi.fn(async () => 'text') };
    const a = createSeqContext({ ports: { assets: () => ({ reader, dispose }) } });
    expect(await a.assets.readText('x', noSignal)).toBe('text');
    const b = createSeqContext({ ports: { assets: () => ({ reader, dispose, disposeWithViewer: false }) } });
    a.retire();
    b.retire();
    expect(dispose).toHaveBeenCalledTimes(1);
  });

  it('a factory that throws disposes the ports already made, last first, and a structure port with them', () => {
    const order: string[] = [];
    const assets = () => ({ reader: { readText: async () => '' }, dispose: () => order.push('assets') });
    const text = () => ({ measure: () => 1, dispose: () => order.push('text') });
    const structure = () => {
      throw new Error('no structure');
    };
    expect(() => createSeqContext({ ports: { assets, textMetrics: text, structure } })).toThrow('no structure');
    expect(order).toEqual(['text', 'assets']);
    // a structure port made is disposed with the context
    const dispose = vi.fn();
    createSeqContext({ ports: { structure: () => ({ runtime: {}, dispose }) } }).retire();
    expect(dispose).toHaveBeenCalledTimes(1);
  });

  it('measures text through the port', () => {
    const measure = vi.fn((font: string, text: string) => text.length * 7);
    const ctx = createSeqContext({ ports: { textMetrics: () => ({ measure }) } });
    expect(ctx.text.measure('12px sans', 'abc')).toBe(21);
    expect(measure).toHaveBeenCalledWith('12px sans', 'abc');
  });
});
