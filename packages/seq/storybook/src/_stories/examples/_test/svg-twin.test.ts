/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ComplexViews } from '@molstar/seqx-examples';
import { nodeFixtureAssets } from '@molstar/seqx-examples/node-assets';
import { DomHost } from '@molstar/ui/dom';
import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { composeEditing, type EditingSession } from '../../features/editing/compose.js';
import { followSvg } from '../../follow-svg.js';
import { hideLigandsPast } from '../complex/compose.js';

const raf = DomHost.requestAnimationFrame;
let session: EditingSession | undefined;

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

const strip = (x: string) => x.replace(/t\d+-/g, '');

describe('the SVG twin follows the canvas', () => {
  it('Ligands: after hide and after undo the twin is the export of what the canvas shows', async () => {
    const view = ComplexViews.find((v) => v.id === 'ligands')!;
    const s = composeEditing(view, { ports: { assets: () => nodeFixtureAssets() }, delayMs: 0 });
    session = s;
    s.viewer.resize(720, 900);
    const host = { innerHTML: '' } as unknown as HTMLElement;
    const svg = followSvg(s.viewer, host, { prefix: 't' });
    await s.viewer.load(view.build());
    await s.viewer.settle();
    await vi.advanceTimersByTimeAsync(600);
    await svg.draw();
    const fresh = async () =>
      strip((await Promise.all([s.viewer.renderSvg({ idPrefix: 't1-' }), vi.advanceTimersByTimeAsync(400)]))[0].svg);
    const before = strip(host.innerHTML);
    await hideLigandsPast(s, 2);
    // before the view rests (`canvas.settleMs`): the twin is not a rest behind the canvas
    await vi.advanceTimersByTimeAsync(s.viewer.config.canvas.settleMs / 2);
    const hidden = strip(host.innerHTML);
    expect(hidden).not.toBe(before);
    expect(hidden).toBe(await fresh());
    await s.stack.undo();
    await vi.advanceTimersByTimeAsync(s.viewer.config.canvas.settleMs / 2);
    expect(strip(host.innerHTML)).toBe(before);
    expect(strip(host.innerHTML)).toBe(await fresh());
  });
});

describe('the example root’s twin', () => {
  it('Family: after a fold the twin is the export of what the canvas shows', async () => {
    const { openExample } = await import('./draw-checks.js');
    const { GlobinsViews } = await import('@molstar/seqx-examples');
    const v = await openExample(GlobinsViews.find((x) => x.id === 'family')!, { fixtures: true });
    const host = { innerHTML: '' } as unknown as HTMLElement;
    const svg = followSvg(v, host, { prefix: 't' });
    await svg.draw();
    const before = strip(host.innerHTML);
    await v.execute('clades.collapseAll', { track: '/views/main/family/rows' });
    // before the view rests
    await vi.advanceTimersByTimeAsync(v.config.canvas.settleMs / 2);
    const after = strip(host.innerHTML);
    const fresh = strip(
      (await Promise.all([v.renderSvg({ idPrefix: 't1-' }), vi.advanceTimersByTimeAsync(400)]))[0].svg,
    );
    expect(after).not.toBe(before);
    expect(after).toBe(fresh);
    svg.dispose();
    v.retire();
  });
});
