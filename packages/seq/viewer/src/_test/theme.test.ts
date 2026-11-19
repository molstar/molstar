// @vitest-environment jsdom
/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { uiHistory } from '../ui-history.js';
import { DarkTokens, LightTokens } from '@molstar/seqx-core/render';
import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { RecordingContext } from '../canvas/_test/recording.js';
import { cssColorToHex, readCssTokens, themeOf } from '../canvas/css-tokens.js';
import { viewerItem } from '../config.js';
import { DefaultSeqViewerSpec } from '../default-spec.js';
import { createSeqViewer, type SeqViewer } from '../viewer.js';
import { proteinOnGene } from './fixture.js';
import { type Frames, installFrames } from './frames.js';

// The token table from CSS (_spec/README.md "Canvas2D backend", Tokens) in a DOM (jsdom, with a stubbed computed
// style), and a theme switch through the attribute.

// the observers counted here are the canvas mount's: the toolbar's are `toolbar.test.ts`'s
const canvasOnly = () => ({ ...DefaultSeqViewerSpec(), ui: { canvas: true as const, history: uiHistory } });

let frames: Frames;
let viewer: SeqViewer | undefined;
/** The custom properties the stubbed computed style holds, by the theme on the document element. */
let css: Record<string, Record<string, string>> = {};
const observers = { observed: 0, disconnected: 0 };
const RealMutationObserver = globalThis.MutationObserver;

class CountingObserver extends RealMutationObserver {
  override observe(target: Node, options?: MutationObserverInit) {
    observers.observed++;
    super.observe(target, options);
  }
  override disconnect() {
    observers.disconnected++;
    super.disconnect();
  }
}

beforeEach(() => {
  frames = installFrames();
  css = {};
  observers.observed = observers.disconnected = 0;
  vi.stubGlobal('MutationObserver', CountingObserver);
  vi.stubGlobal(
    'ResizeObserver',
    class {
      observe() {}
      unobserve() {}
      disconnect() {}
    },
  );
  vi.spyOn(HTMLElement.prototype, 'clientWidth', 'get').mockReturnValue(600);
  vi.spyOn(HTMLCanvasElement.prototype, 'getContext').mockImplementation(function (this: HTMLCanvasElement) {
    return ((this as any).__ctx ??= new RecordingContext());
  });
  vi.spyOn(document.defaultView!, 'getComputedStyle').mockImplementation((() => {
    const theme = document.documentElement.getAttribute('data-msui-theme') ?? 'light';
    const props = css[theme] ?? {};
    return { getPropertyValue: (n: string) => props[n] ?? '', colorScheme: '' } as unknown as CSSStyleDeclaration;
  }) as typeof getComputedStyle);
});

afterEach(() => {
  viewer?.retire();
  viewer = undefined;
  frames.restore();
  vi.restoreAllMocks();
  vi.unstubAllGlobals();
  document.documentElement.removeAttribute('data-msui-theme');
  document.body.innerHTML = '';
});

async function mount(items: [any, unknown][] = []) {
  const target = document.createElement('div');
  document.body.appendChild(target);
  const v = createSeqViewer({ ...canvasOnly(), config: items }, target);
  viewer = v;
  await v.load(proteinOnGene());
  await v.settle();
  frames.flush();
  return v;
}

const tick = () => new Promise<void>((r) => setTimeout(r, 0));

describe('the token reader', () => {
  it('reads the computed style, and falls back per token to the theme the root shows', () => {
    css.light = { '--msui-bg': '#102030', '--seq-track': 'rgb(1, 2, 3)', '--msui-font-mono': ' monospace ' };
    const t = readCssTokens(document.body);
    expect(t.colors['msui-bg']).toBe(0x102030ff);
    expect(t.colors['seq-track']).toBe(0x010203ff);
    expect(t.colors['msui-fg']).toBe(LightTokens.colors['msui-fg']);
    expect(t.fonts.mono).toBe('monospace');
    expect(t.fonts.sans).toBe(LightTokens.fonts.sans);

    document.documentElement.dataset.msuiTheme = 'dark';
    css.dark = {};
    const d = readCssTokens(document.body);
    expect(d.theme).toBe('dark');
    expect(d.version).toBe(DarkTokens.version);
  });

  it('takes the theme of the nearest ancestor over the document element', () => {
    document.documentElement.dataset.msuiTheme = 'light';
    const holder = document.createElement('section');
    holder.dataset.msuiTheme = 'dark';
    const child = document.createElement('div');
    holder.appendChild(child);
    document.body.appendChild(holder);
    expect(themeOf(child)).toBe('dark');
    expect(themeOf(document.body)).toBe('light');
  });

  it('reads the colors a browser resolves to', () => {
    expect(cssColorToHex('rgb(255 0 0 / 50%)')).toBe('#ff000080');
    expect(cssColorToHex('color(srgb 0 1 0)')).toBe('#00ff00ff');
    expect(cssColorToHex('color-mix(in oklab, red, blue)')).toBeUndefined();
  });
});

describe('the mounted viewer', () => {
  it('draws with the root’s table', async () => {
    css.light = { '--msui-bg': '#102030' };
    const v = await mount();
    expect(v.tokens.colors['msui-bg']).toBe(0x102030ff);
  });

  it('a switch of data-msui-theme gives a new table, a frame and a new version; the SVG follows', async () => {
    document.documentElement.dataset.msuiTheme = 'light';
    css.light = {};
    css.dark = { '--msui-bg': '#010203' };
    const v = await mount();
    const before = v.tokens.version;
    const moved: number[] = [];
    v.tokensVersion.listen((x) => moved.push(x));
    const drawn = v.counts.frames;
    document.documentElement.dataset.msuiTheme = 'dark';
    await tick();
    expect(v.tokens.version).not.toBe(before);
    expect(v.tokens.theme).toBe('dark');
    expect(moved).toEqual([v.tokens.version]);
    frames.flush();
    expect(v.counts.frames).toBeGreaterThan(drawn);
    const svg = await v.renderSvg();
    expect(svg.svg.toLowerCase()).toContain('#010203');
  });

  it('an unchanged theme keeps the version and asks for no frame', async () => {
    document.documentElement.dataset.msuiTheme = 'dark';
    css.dark = {};
    const v = await mount();
    const tokens = v.tokens;
    const moved: number[] = [];
    v.tokensVersion.listen((x) => moved.push(x));
    const drawn = v.counts.frames;
    // the attribute is written again with the same theme
    document.documentElement.setAttribute('data-msui-theme', 'dark');
    await tick();
    frames.flush();
    expect(v.tokens).toBe(tokens);
    expect(moved).toEqual([]);
    expect(v.counts.frames).toBe(drawn);
  });

  it('retire disconnects the observer', () => {
    const holder = document.createElement('div');
    holder.dataset.msuiTheme = 'light';
    document.body.appendChild(holder);
    const target = document.createElement('div');
    holder.appendChild(target);
    const v = createSeqViewer(canvasOnly(), target);
    // one observer on the document element and the nearest ancestor that has the attribute
    expect(observers.observed).toBe(2);
    v.retire();
    expect(observers.disconnected).toBe(1);
  });

  it('theme.source: config keeps the config’s table and observes nothing', async () => {
    css.light = { '--msui-bg': '#102030' };
    css.dark = { '--msui-bg': '#102030' };
    const v = await mount([[viewerItem('theme.source'), 'config']]);
    expect(v.tokens).toBe(v.context.tokens);
    expect(observers.observed).toBe(0);
    document.documentElement.dataset.msuiTheme = 'dark';
    await tick();
    expect(v.tokens).toBe(v.context.tokens);
  });
});
