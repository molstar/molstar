// @vitest-environment jsdom
/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The stories' search box (support/search-box.ts): built on the public surface alone (its imports are packages'
// exports, never a path into a package's sources: an application can build the same), and its DOM: the text searched
// once it rests, the count, Enter and Shift+Enter stepping, Esc clearing.

import { readFileSync } from 'node:fs';
import { dirname, join } from 'node:path';
import { fileURLToPath } from 'node:url';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { DefaultSeqViewerSpec } from '@molstar/seqx-viewer/default-spec';
import { createSeqViewer, type SeqViewer } from '@molstar/seqx-viewer/viewer';
import { afterEach, describe, expect, it, vi } from 'vitest';
import { mountSearchBox, SearchRestMs } from '../search-box.js';

const here = dirname(fileURLToPath(import.meta.url));
const repo = join(here, '../../../../../../..');

/** Whether `@molstar/<pkg>/<sub>` is an export of its package (an entry, or its `./*` wildcard; never one mapped to null). */
function exported(specifier: string): boolean {
  const m = /^(@molstar\/[^/]+)(\/.*)?$/.exec(specifier);
  if (!m) return false;
  const pkgs: Record<string, string> = {
    '@molstar/core': 'packages/core',
    '@molstar/ui': 'packages/ui',
    '@molstar/seqx-core': 'packages/seq/core',
    '@molstar/seqx-viewer': 'packages/seq/viewer',
  };
  const folder = pkgs[m[1]];
  if (!folder) return false;
  const exports = JSON.parse(readFileSync(join(repo, folder, 'package.json'), 'utf8')).exports as Record<
    string,
    unknown
  >;
  const sub = `.${m[2] ?? ''}`;
  if (Object.hasOwn(exports, sub)) return exports[sub] !== null;
  return sub.split('/').some((_, i, parts) => {
    const prefix = `${parts.slice(0, i).join('/')}/*`;
    return i > 0 && Object.hasOwn(exports, prefix) && exports[prefix] !== null;
  });
}

describe('the search box’s imports', () => {
  it('are packages’ public exports only (and the stories’ own toolbar)', () => {
    const text = readFileSync(join(here, '../search-box.ts'), 'utf8');
    const specifiers = [...text.matchAll(/^import .* from '([^']+)';$/gm)].map((m) => m[1]);
    expect(specifiers.length).toBeGreaterThan(5);
    const wrong = specifiers.filter((s) => (s.startsWith('.') ? s !== './toolbar.js' : !exported(s)));
    expect(wrong).toEqual([]);
    // never the viewer's internals by path
    expect(specifiers.some((s) => s.includes('/src/'))).toBe(false);
  });
});

const viewers: SeqViewer[] = [];
afterEach(() => {
  for (const v of viewers.splice(0)) v.retire();
  document.body.innerHTML = '';
  vi.useRealTimers();
});

function spec() {
  const b = createSeqSpec({ title: 'box' });
  const p = b.data.sequence('p', { alphabet: 'protein', letters: 'MEEPMEEPMEEP' });
  b.views.view('main', (v) => v.section('p', { of: p }).track('l').layer('l', { type: 'letters' }));
  return b.build();
}

async function make(initial?: string) {
  vi.useFakeTimers({ toFake: ['setTimeout', 'clearTimeout'] });
  const v = createSeqViewer({ ...DefaultSeqViewerSpec(), ui: undefined });
  viewers.push(v);
  v.resize(600);
  await v.load(spec());
  await vi.advanceTimersByTimeAsync(64);
  const el = document.createElement('div');
  document.body.append(el);
  const off = mountSearchBox(v, el, initial ? { initial } : {});
  await vi.advanceTimersByTimeAsync(64);
  const input = el.querySelector<HTMLInputElement>('input')!;
  const count = () => el.querySelector('.seqx-search-count')!.textContent;
  const type = async (text: string) => {
    input.value = text;
    input.dispatchEvent(new Event('input', { bubbles: true }));
  };
  /** The clock run while a run's chunks land (their imports are real I/O, never the fake clock's) until it ends. */
  const ran = async () => {
    for (let i = 0; i < 2000; i++) {
      await vi.advanceTimersByTimeAsync(16);
      await new Promise((r) => setImmediate(r));
      const s = v.search!;
      if (s.status.value === 'done' || s.status.value === 'failed') return;
    }
  };
  const key = (k: string, shiftKey = false) =>
    input.dispatchEvent(new KeyboardEvent('keydown', { key: k, shiftKey, bubbles: true }));
  return { v, el, off, input, count, type, key, ran };
}

describe('the search box', () => {
  it('searches the text once it rests 200 ms, says the count, steps with Enter and Shift+Enter', async () => {
    const b = await make();
    await b.type('EP');
    await vi.advanceTimersByTimeAsync(SearchRestMs - 1);
    expect(b.v.search!.query.value).toBeUndefined();
    await b.ran();
    expect(b.count()).toBe('3 hits');
    b.key('Enter');
    await vi.advanceTimersByTimeAsync(16);
    expect(b.count()).toBe('1 of 3');
    b.key('Enter', true);
    await vi.advanceTimersByTimeAsync(16);
    expect(b.count()).toBe('3 of 3');
  });

  it('typing again before the rest cancels the earlier text: one run, of the last', async () => {
    const b = await make();
    await b.type('E');
    await vi.advanceTimersByTimeAsync(100);
    await b.type('EE');
    await vi.advanceTimersByTimeAsync(SearchRestMs);
    await b.ran();
    expect(b.v.search!.query.value?.pattern).toBe('EE');
    expect(b.count()).toBe('3 hits');
  });

  it('an initial pattern is searched once the viewer has its spec; Esc empties the box and clears the search', async () => {
    const b = await make('MEEP');
    await b.ran();
    expect(b.count()).toBe('3 hits');
    // the text input's own Escape empties it (a search field), and the box clears the search
    b.input.dispatchEvent(new KeyboardEvent('keydown', { key: 'Escape', bubbles: true }));
    await vi.advanceTimersByTimeAsync(16);
    expect(b.input.value).toBe('');
    expect(b.v.search!.query.value).toBeUndefined();
    expect(b.count()).toBe('');
  });

  it('a protein-only view shows no strand switch; no hits are said as such', async () => {
    const b = await make();
    const labels = [...b.el.querySelectorAll<HTMLElement>('button')];
    const shown = (text: string) =>
      labels.some((x) => x.textContent === text && !x.closest<HTMLElement>('span[style*="display: none"]'));
    expect(shown('IUPAC')).toBe(true);
    expect(shown('Forward')).toBe(false);
    await b.type('WWW');
    await vi.advanceTimersByTimeAsync(SearchRestMs);
    await b.ran();
    expect(b.count()).toBe('No hits');
  });
});
