/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { RingPlacement } from '@molstar/seqx-core/render';
import { loadGenome } from '@molstar/seqx-core/resolve';
import { rondoViewer } from '@molstar/seqx-ext-rondo/viewer';
import type { SeqViewer } from '@molstar/seqx-viewer/viewer';
import { afterEach, beforeAll, describe, expect, it, vi } from 'vitest';
import { openExample } from '../../../examples/_test/draw-checks.js';
import { RingSectionId, ringStoryView } from '../views.js';

// Features/Genome Ring › The Genome At The Fit, headless: the ring nobody follows takes Cmd/Ctrl+wheel as the cut
// (_spec/README.md "Input", The ring's zoom is a cut): in over a chromosome opens it, out shows the genome again, over
// the hole it says what the ring shows and changes nothing.

let viewer: SeqViewer | undefined;

beforeAll(async () => {
  await loadGenome();
});
afterEach(() => {
  viewer?.retire();
  viewer = undefined;
  vi.useRealTimers();
});

async function open() {
  viewer = await openExample(ringStoryView('ring'), {
    delayMs: 0,
    wait: 3000,
    fixtures: true,
    extensions: [rondoViewer],
  });
  viewer.draw();
  viewer.input.mac = true;
  return viewer;
}

const address = (v: SeqViewer) =>
  v.view.resolved.value!.sections.find((s) => s.address.endsWith(`/${RingSectionId}`))!.address;
const placement = (v: SeqViewer) => v.placementOf(address(v)) as RingPlacement;

/** The page point on the ring's first band at the middle of block `k`. */
function over(v: SeqViewer, k: number) {
  const p = placement(v);
  const b = v.view.resolved.value!.sections.find((s) => s.address === address(v))!.blocks![k];
  const s = p.strips;
  p.warp.map(p.frameX((b.offset ?? 0) + 1000), s.bandY[0] + s.bandH[0] / 2);
  return { x: p.warp.px, y: p.warp.py };
}

const included = (v: SeqViewer) => {
  const find = (n: { id?: string; children?: unknown[]; params?: Record<string, unknown> }): unknown => {
    if (n.id === RingSectionId) return n.params?.include;
    for (const c of (n.children ?? []) as (typeof n)[]) {
      const r = find(c);
      if (r !== undefined) return r;
    }
    return undefined;
  };
  return find(v.spec.root.value as never) as string[] | undefined;
};

const zoom = (v: SeqViewer, at: { x: number; y: number }, deltaY: number) =>
  v.input.wheel(at.x, at.y, { deltaX: 0, deltaY, deltaMode: 0, metaKey: true, ctrlKey: false });

const rest = async (v: SeqViewer) => {
  await vi.advanceTimersByTimeAsync(1000);
  v.draw();
};

describe('The Genome At The Fit: the ring’s Cmd+wheel', () => {
  it('in over a chromosome opens it, out shows the genome again, over the hole says so once', async () => {
    const v = await open();
    const said: string[] = [];
    v.notice.listen(({ text }) => said.push(text));
    const n = v.view.resolved.value!.sections.find((s) => s.address === address(v))!.blocks!.length;
    expect(n).toBeGreaterThan(20);
    expect(included(v)).toBeUndefined();
    // the hole, the ring's middle: nothing changes
    const { box } = placement(v);
    expect(zoom(v, { x: box.x + box.width / 2 + 6, y: box.y + box.height / 2 + 4 }, -100)).toBe(true);
    expect(included(v)).toBeUndefined();
    expect(said).toHaveLength(1);
    await rest(v);
    // in over chr4 (the fourth block)
    expect(zoom(v, over(v, 3), -100)).toBe(true);
    await rest(v);
    expect(included(v)).toEqual(['chr4']);
    // out over it: every block
    expect(zoom(v, over(v, 0), 100)).toBe(true);
    await rest(v);
    expect(included(v)).toBeUndefined();
    expect(v.view.resolved.value!.sections.find((s) => s.address === address(v))!.blocks!.length).toBe(n);
  });
});
