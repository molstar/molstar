/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The hover of an extension's layer type costs nothing once warm (performance.md): with the `demo` extension's
// lollipop hovered, a frame drawn again (the view outlines the mark the type reports, `LayerType.mark`, into one reused
// box) makes no object. By hand, never in CI: `npx vitest run -c packages/seq/storybook/vitest.perf.config.ts`.

import { createSeqViewer } from '@molstar/seqx-viewer/viewer';
import { measureBytesPerStep, quickCase } from '@molstar/test-harness/quick';
import { describe, expect, inject, it } from 'vitest';
import { DemoExtension, extensibilityComposition, extensibilitySpec } from '../compose.js';

describe('the hover of an extension’s layer type', () => {
  it('a lollipop hovered, the frame drawn again', async ({ task }) => {
    const v = createSeqViewer(extensibilityComposition(false));
    v.resize(900);
    await v.load(extensibilitySpec());
    await v.settle();
    await v.extend(DemoExtension);
    await v.settle();
    v.draw();
    /** A point on the first item of the layers of a type, found by hit-testing down the band. */
    const point = (type: string, key: string) => {
      const layer = v
        .frame!.layout.sections.flatMap((t) => t.tracks.flatMap((k) => k.layers))
        .find((l) => l.type.name === type)!;
      const it = layer.type.items!(layer.lay).find((i) => i.key === key)!;
      const x = (it.x0 + it.x1) / 2;
      const y = [...Array(120).keys()].map((k) => k * 1).find((k) => v.hit(x, k)?.element?.key === key);
      expect(y, `${type} ${key}`).toBeDefined();
      return [x, y!] as const;
    };
    const frame = () => void v.draw();
    // a built-in item hovered (what a frame of this viewer costs with a hover: its own card, announcements), then the
    // lollipop: the outline the view draws for it adds nothing
    v.input.move(...point('ranges', 'TAD1'));
    v.draw();
    const builtIn = Math.min(measureBytesPerStep(frame).bytesPerStep, measureBytesPerStep(frame).bytesPerStep);
    v.input.move(...point('demo:lollipops', 'S6'));
    v.draw();
    const { r } = quickCase(task, 'extension hover frame', frame, {});
    const hovered = Math.min(r.bytesPerStep, measureBytesPerStep(frame).bytesPerStep);
    console.info(
      `extension hover frame: ${hovered.toFixed(0)} bytes hovered, ${builtIn.toFixed(0)} a built-in's (${r.note})`,
    );
    expect(hovered - builtIn, r.note).toBeLessThanOrEqual(inject('perfThreshold'));
    v.retire();
  });
});
