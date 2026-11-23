/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The ribbons' time bounds (a ring of 23 blocks): 2,000 links laid out (stems and ribbons) under 4 ms and drawn under
// 8 ms; 12,000 links at the fit laid out under 16 ms once, a frame's draw under 8 ms, a hover (a hundred hits and the
// hover pass) under 2 ms. The layout timed is the type's own, with the context the renderer gave it (the ring placed,
// the store as it is); the bounds are the heavy lane's (`itWithin`, best of three), the quick lane runs them unbounded.

import { defineExtension, type LayerType, type LayoutContext, type PassWriter } from '@molstar/seqx-core/extension';
import { DrawLists } from '@molstar/seqx-core/render';
import { itWithin } from '@molstar/test-harness/timing';
import { describe, expect } from 'vitest';
import { RondoHues } from '../hues.js';
import type { RondoLay, RondoPrep } from '../lay.js';
import { RibbonsType } from '../ribbons.js';
import { genomeChromosomes, lay, ringSpec, Size, settle, toyLinks } from './helpers.js';

/** A ring of 23 blocks and `n` links, settled; the ribbons' last layout's arguments kept. */
async function ring(n: number) {
  const last: { prep?: RondoPrep; ctx?: LayoutContext; data?: unknown } = {};
  const Captured: LayerType<RondoPrep, RondoLay> = {
    ...RibbonsType,
    layout(prep, ctx, data, overlay, prev) {
      last.prep = prep;
      last.ctx = ctx;
      last.data = data;
      return RibbonsType.layout(prep, ctx, data, overlay, prev);
    },
  };
  const ext = defineExtension({
    id: 'rondo',
    version: '0.1.0',
    requires: '0.10',
    registry: [
      { kind: 'layer-type', type: Captured },
      { kind: 'scheme', scheme: RondoHues },
    ],
  });
  const chroms = genomeChromosomes(23);
  const s = await settle(ringSpec({ chroms, links: toyLinks(chroms, n, 11) }), [ext]);
  const f = await s.frame({ size: Size, dpr: 2 });
  const l = lay(f)!;
  const layout = () => RibbonsType.layout(last.prep!, last.ctx!, last.data as never, undefined as never, l);
  return { l, layout };
}

function drawAll(l: RondoLay, lists: DrawLists) {
  lists.reset();
  RibbonsType.draw(l, 'data', lists.writer('data') as PassWriter);
  RibbonsType.draw(l, 'labels', lists.writer('labels') as PassWriter);
}

describe('the ribbons in time', () => {
  itWithin('23 blocks and 2,000 links: laid out under 4 ms, drawn under 8 ms', async (within) => {
    const { l, layout } = await ring(2000);
    expect(l.nr).toBeGreaterThan(0);
    const lists = new DrawLists();
    layout();
    drawAll(l, lists);
    within(4, layout);
    within(8, () => drawAll(l, lists));
  });

  itWithin('12,000 links at the fit: laid out under 16 ms, drawn under 8 ms, a hover under 2 ms', async (within) => {
    const { l, layout } = await ring(12_000);
    expect(l.items).toBe(false);
    const lists = new DrawLists();
    within(16, layout);
    drawAll(l, lists);
    within(8, () => drawAll(l, lists));
    const at = new Float64Array(2);
    const hover = l.hits.find((x) => x?.element?.key.toString().startsWith('stem:'));
    within(2, () => {
      for (let i = 0; i < 100; i++) {
        at[0] = l.frame.cx + ((i % 10) / 10 - 0.5) * l.frame.r;
        at[1] = l.frame.cy + (Math.floor(i / 10) / 10 - 0.5) * l.frame.r;
        RibbonsType.hit(l, at);
      }
      RibbonsType.draw(l, 'hover', lists.writer('hover') as PassWriter, { hover });
    });
  });
});
