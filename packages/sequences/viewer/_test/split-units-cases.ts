/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The split rows' hover case, for its heap guard and its quick twin (split-units.heap.ts, split-units.quick.ts): the
// pointer moving over the chains' rows of the split-units document (a highlight on one unit, an accent on the shared
// sequence), each move an overlay frame on the recording backend, and how many times the hover was rebuilt. A fresh
// viewer (`vi.resetModules`).

import { vi } from 'vitest';

export async function unitsHover() {
  vi.resetModules();
  const { lit } = await import('./overlay-fixtures.js');
  const { recording, COPY_BACK } = await import('./composer-fixtures.js');
  const { splitUnitsDoc } = await import('./split-units-fixtures.js');
  const { Dirty } = await import('../src/draw/index.js');
  const { Inputs, PointerType } = await import('../src/index.js');
  const t = await lit(splitUnitsDoc({ highlights: true }), { width: 1000, height: 400, dpr: 2 });
  t.ctx.behaviors.get('hoverLabel')?.remove();
  const be = recording(COPY_BACK);
  t.frame(be, Dirty.all);
  const p = Inputs.pointer();
  p.section = t.ctx.layout.sectionIndex(t.section);
  p.type = PointerType.move;
  const ys = [0, 1, 2].map((r) => t.yOf(t.section, 't-chains', r));
  const frame = (i: number) => {
    p.x = t.xOf(t.section, 0, 8 + (i % 44));
    p.y = ys[i % 3];
    t.ctx.input.pointer(p);
    t.frame(be, Dirty.overlay);
  };
  const channels = (t.composer as unknown as { channels: { rebuilds: { hover: number } } }).channels;
  return { frame, hovers: () => channels.rebuilds.hover };
}
