/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A liftover on the pointer's path (Features/Across sections/Liftover): a hover held on the duplicated gene's image,
// carried across the chain onto the old gene and the other image, its frame drawn again and its card's Places read,
// allocates nothing once warm (the projector's answers kept for the store's version; the hover's projection for its
// version). By hand, never in CI: `npx vitest run -c packages/seq/storybook/vitest.perf.config.ts`.

import type { Hit } from '@molstar/seqx-core/render';
import { loadBridges } from '@molstar/seqx-core/render/bridge-load';
import { loadGenome } from '@molstar/seqx-core/resolve';
import { NewSection } from '@molstar/seqx-examples/across/liftover';
import { AcrossViews } from '@molstar/seqx-examples/across/views';
import { DefaultSeqViewerSpec } from '@molstar/seqx-viewer/default-spec';
import { createSeqViewer } from '@molstar/seqx-viewer/viewer';
import { quickCase } from '@molstar/test-harness/quick';
import { beforeAll, describe, expect, inject, it } from 'vitest';

beforeAll(async () => {
  await loadGenome();
  await loadBridges();
  await import('@molstar/seqx-viewer/bridge-view');
});

const view = AcrossViews.find((v) => v.id === 'across-liftover')!;

describe('a liftover on the pointer’s path', () => {
  it('a hover held on a lifted image: carried across the map, drawn, its places read', async ({ task }) => {
    const base = DefaultSeqViewerSpec();
    const v = createSeqViewer({
      ...base,
      ui: undefined,
      registry: [...(base.registry ?? []), ...view.registry!({ delayMs: 0 })],
    });
    v.resize(1440);
    await v.load(view.build());
    await v.settle();
    v.view.setRange(NewSection, 1_480_000, 1_560_000);
    await v.settle();
    v.draw();
    const sec = v.frame!.layout.sections.find((s) => s.section.address === NewSection)!;
    let hit: Hit | undefined;
    for (let y = sec.y; y < sec.y + sec.height && !hit; y += 2)
      for (let x = 0; x < 1440 && !hit; x += 4) {
        const h = v.hit(x, y);
        if (h?.element?.key === 'dup') hit = { ...h, section: undefined } as unknown as Hit;
      }
    expect(hit).toBeDefined();
    v.hover.set(hit, 'test');
    const regions = v.store.regions;
    const model = v.view.model.value;
    const resolved = v.view.resolved.value!;
    const own = v.hover.location()!.parts[0];
    const at = [own.ranges[0], own.ranges[1]] as const;
    const system = own.system;
    expect(regions.places?.(model, resolved, system, at[0], at[1])).toMatch(/^2 places/);
    const lit = (v.hover.projectedAll()?.parts ?? []).filter((p) => p.system.startsWith('/data/v2/'));
    expect(lit.length).toBe(2);
    const step = () => {
      v.hover.projectedAll();
      regions.places?.(model, resolved, system, at[0], at[1]);
      v.draw();
    };
    const { r } = quickCase(task, 'liftover held hover', step);
    v.retire();
    console.info(`liftover held hover: ${r.bytesPerStep.toFixed(0)} bytes a step (${r.note})`);
    expect(r.bytesPerStep, r.note).toBeLessThanOrEqual(inject('perfThreshold'));
  });
});
