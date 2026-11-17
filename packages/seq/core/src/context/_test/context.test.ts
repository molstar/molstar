/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Task } from '@molstar/core/task';
import type { AssetManager } from '@molstar/core/util/assets';
import { type FormatRuntime, SeqRuntime } from '@molstar/seqx-core/data';
import { createSeqSpec } from '@molstar/seqx-core/lang';
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
    v.section('protein', { of: p }).track('letters').layer({ type: 'letters' });
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
