/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The ribbons' allocation targets ("a ribbon frame"): a warm draw of the data pass, a hit and a hover pass (a stem's
// ribbons lit over the wash, a ribbon outlined) of a ring of 23 blocks and 12,000 links allocate nothing. By hand,
// never in CI: `npx vitest run -c packages/seq/ext-rondo/vitest.perf.config.ts`.

import type { Hit, PassWriter } from '@molstar/seqx-core/extension';
import { DrawLists } from '@molstar/seqx-core/render';
import { quickCase } from '@molstar/test-harness/quick';
import { describe, expect, inject, it } from 'vitest';
import { RibbonsType } from '../ribbons.js';
import { genomeChromosomes, lay, ringSpec, Size, settle, toyLinks } from './helpers.js';

async function ring() {
  const chroms = genomeChromosomes(23);
  const s = await settle(ringSpec({ chroms, links: toyLinks(chroms, 12_000, 11) }));
  return lay(await s.frame({ size: Size, dpr: 2 }))!;
}

describe('the ribbons draw and hover with nothing allocated once warm', () => {
  it('a warm draw of the data pass', async ({ task }) => {
    const l = await ring();
    const lists = new DrawLists();
    const { r } = quickCase(task, 'ribbons draw', () => {
      lists.reset();
      RibbonsType.draw(l, 'data', lists.writer('data') as PassWriter);
    });
    console.info(`ribbons draw: ${r.bytesPerStep.toFixed(0)} bytes a frame (${r.note})`);
    expect(r.bytesPerStep, r.note).toBeLessThanOrEqual(inject('perfThreshold'));
  });

  it('a hit', async ({ task }) => {
    const l = await ring();
    const at = new Float64Array(2);
    const { r } = quickCase(task, 'ribbons hit', (i) => {
      at[0] = l.frame.cx + (((i * 37) % 100) / 100 - 0.5) * l.frame.r;
      at[1] = l.frame.cy + (((i * 53) % 100) / 100 - 0.5) * l.frame.r;
      RibbonsType.hit(l, at);
    });
    console.info(`ribbons hit: ${r.bytesPerStep.toFixed(0)} bytes a hit (${r.note})`);
    expect(r.bytesPerStep, r.note).toBeLessThanOrEqual(inject('perfThreshold'));
  });

  it('the hover pass (stems and ribbons in turn)', async ({ task }) => {
    const l = await ring();
    const lists = new DrawLists();
    const hovers = l.hits.filter((h): h is Hit => !!h);
    const emphasis: { hover: Hit | undefined } = { hover: undefined };
    const { r } = quickCase(task, 'ribbons hover', (i) => {
      emphasis.hover = hovers[i % hovers.length];
      lists.reset();
      RibbonsType.draw(l, 'hover', lists.writer('hover') as PassWriter, emphasis);
    });
    console.info(`ribbons hover: ${r.bytesPerStep.toFixed(0)} bytes a step (${r.note})`);
    expect(r.bytesPerStep, r.note).toBeLessThanOrEqual(inject('perfThreshold'));
  });
});
