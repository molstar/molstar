/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { rangesLocation } from '@molstar/seqx-core/query';
import { P53Views } from '@molstar/seqx-examples';
import type { SeqViewer } from '@molstar/seqx-viewer/viewer';
import { DomHost } from '@molstar/ui/dom';
import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { composeViewer } from '../../../features/selection/compose.js';
import { highlightOf, nightingaleData } from '../compose.js';

// The Nightingale story, headless on a fake clock, composed as the story composes it: the viewer's resolved view and
// store feed the elements' shapes; the viewer's hover and selection are the elements' `highlight`.

let viewer: SeqViewer | undefined;
const raf = DomHost.requestAnimationFrame;
beforeEach(() => {
  DomHost.requestAnimationFrame = (fn) => {
    setTimeout(() => fn(performance.now()), 16);
  };
});
afterEach(() => {
  DomHost.requestAnimationFrame = raf;
  viewer?.retire();
  viewer = undefined;
  vi.useRealTimers();
});

async function make() {
  vi.useFakeTimers({ toFake: ['setTimeout', 'clearTimeout', 'Date', 'performance'] });
  const v = composeViewer(P53Views.find((x) => x.id === 'tracks')!, { noRails: true });
  viewer = v;
  v.resize(900);
  await v.load(P53Views.find((x) => x.id === 'tracks')!.build());
  await vi.advanceTimersByTimeAsync(400);
  return v;
}

describe('Nightingale story, headless', () => {
  it('the elements are fed from the viewer, in each theme’s colors', async () => {
    const v = await make();
    const light = await nightingaleData(v, 'light');
    const dark = await nightingaleData(v, 'dark');
    expect(light.read.letters).toHaveLength(393);
    expect(light.features).toHaveLength(6);
    expect(light.lines[0].values).toHaveLength(391);
    expect(light.colors.text).not.toBe(dark.colors.text);
  });

  it('hover and selection are the highlight, on the section’s system only', async () => {
    const v = await make();
    const { read } = await nightingaleData(v, 'light');
    expect(highlightOf(v, read.system)).toBe('');
    v.hover.setLocation(rangesLocation(read.system, [94, 95]), {});
    expect(highlightOf(v, read.system)).toBe('95:95');
    await v.execute('address.select', { address: `${read.system}/(range :from 94 :to 100)` });
    expect(highlightOf(v, read.system).split(',').length).toBe(2);
    v.hover.clear(v.hover.owner);
    expect(highlightOf(v, 'elsewhere')).toBe('');
  });
});
