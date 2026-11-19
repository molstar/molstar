/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { printSexpr } from '@molstar/core/repr/codecs/sexpr';
import { type Node, resolveAddress } from '@molstar/core/repr/tree';
import { FeatureViews, P53Views } from '@molstar/seqx-examples';
import { LetterSchemes } from '@molstar/seqx-examples/p53/themes';
import { DomHost } from '@molstar/ui/dom';
import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { composeEditing, type EditingSession } from '../../features/editing/compose.js';
import {
  cancelScheme,
  Letters,
  pickerSchemes,
  pickScheme,
  previewScheme,
  setSignalDomain,
  signalDomains,
} from '../p53/themes-compose.js';

// The p53 themes stories, headless on a fake clock, composed as the editing story composes them: a scheme previewed
// while hovered (no step) and picked in one patch; the shared scale's domain set for both signals in one step and undone
// in one; the export's legends below the data.

const view = FeatureViews.find((v) => v.id === 'themes-editing')!;
let session: EditingSession | undefined;

const raf = DomHost.requestAnimationFrame;
beforeEach(() => {
  DomHost.requestAnimationFrame = (fn) => {
    setTimeout(() => fn(performance.now()), 16);
  };
  vi.useFakeTimers({ toFake: ['setTimeout', 'clearTimeout', 'Date', 'performance'] });
});

afterEach(() => {
  DomHost.requestAnimationFrame = raf;
  session?.dispose();
  session = undefined;
  vi.useRealTimers();
});

async function make() {
  const s = composeEditing(view);
  session = s;
  s.viewer.resize(900);
  await s.viewer.load(view.build());
  await vi.advanceTimersByTimeAsync(400);
  return s;
}

/** The classes track's letters' color as shown (the preview included), printed. */
const lettersColor = (s: EditingSession) => {
  const layer = resolveAddress(s.viewer.spec.shown.value, Letters.split('/').filter(Boolean));
  const c = layer?.params?.color;
  return c === undefined ? '' : printSexpr(c as Node);
};

describe('p53 Themes, headless', () => {
  it('both views are listed, the same spec', () => {
    expect(P53Views.find((v) => v.id === 'themes')?.build).toBe(view.build);
  });

  it('a scheme previewed is no step, picked is one; undo restores the default', async () => {
    const s = await make();
    const offered = pickerSchemes(s);
    for (const name of LetterSchemes) expect(offered).toContain(name);
    expect(offered).not.toContain('nucleotide');

    const before = lettersColor(s);
    expect(previewScheme(s, 'zappo')).toBe(true);
    await vi.advanceTimersByTimeAsync(100);
    expect(lettersColor(s)).toContain('zappo');
    expect(s.stack.steps.value).toHaveLength(0);
    cancelScheme(s);
    await vi.advanceTimersByTimeAsync(100);
    expect(lettersColor(s)).toBe(before);

    previewScheme(s, 'taylor');
    await pickScheme(s, 'taylor');
    await vi.advanceTimersByTimeAsync(100);
    expect(lettersColor(s)).toContain('taylor');
    expect(s.stack.steps.value).toHaveLength(1);
    await s.stack.undo();
    await vi.advanceTimersByTimeAsync(100);
    expect(lettersColor(s)).toBe(before);
  });

  it('the shared scale’s domain: both signals follow in one step, undone in one', async () => {
    const s = await make();
    const [a, b] = signalDomains(s);
    // open at the top: the union of both signals (the accessibility's range is the larger)
    expect(a).toBeDefined();
    expect(a).toEqual(b);
    expect(a![0]).toBe(0);
    expect(a![1]).toBeGreaterThan(1);

    await setSignalDomain(s, [0, 0.5]);
    await vi.advanceTimersByTimeAsync(400);
    expect(signalDomains(s)).toEqual([
      [0, 0.5],
      [0, 0.5],
    ]);
    expect(s.stack.steps.value).toHaveLength(1);

    await s.stack.undo();
    await vi.advanceTimersByTimeAsync(400);
    expect(signalDomains(s)).toEqual([a, b]);
  });

  it('the legends: one per vocabulary, the shared scale’s merged, and the export draws them below', async () => {
    const s = await make();
    const legends = s.viewer.legends.legends.value;
    const signal = legends.filter((l) => l.scale === '/views/signal');
    expect(signal).toHaveLength(1);
    expect(signal[0].layers).toHaveLength(2);
    expect(legends.some((l) => l.field === 'consequence')).toBe(true);
    const p = s.viewer.renderSvg({ idPrefix: 'th' });
    await vi.advanceTimersByTimeAsync(300);
    const svg = (await p).svg;
    expect(svg).toContain('data-seq-scale-');
    expect(svg).toMatch(/missense/);
  });
});
