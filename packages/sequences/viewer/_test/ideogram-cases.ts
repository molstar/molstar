/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The ideogram navigator's case (ST5), for its heap guard and its quick twin (ideogram.heap.ts, ideogram.quick.ts):
// hg38's committed cytobands on the streamed genome's first three chromosomes (about 200 bands); a frame redraws the
// strip by hand (the nav pass, as a layout or skin change would) and the window over it, moved between two places.
// Each call loads the viewer afresh (`vi.resetModules`).

import { vi } from 'vitest';

export async function ideogramFrame() {
  vi.resetModules();
  const { composed, recording } = await import('./composer-fixtures.js');
  const { GRCH38, hg38Streams, streamedDoc } = await import('../../core/_test/format/stream-provider.js');
  const { ideogramNavigator } = await import('../src/index.js');
  const { readFileSync } = await import('node:fs');
  const { GENOME } = await import('./stream-fixtures.js');
  const bands = JSON.parse(readFileSync(new URL('./fixtures/hg38-cytobands.json', import.meta.url), 'utf8'));
  const s = hg38Streams({ chromosomes: GRCH38.slice(0, 3), letters: true });
  const t = await composed(streamedDoc({ values: false }), {
    width: 1000,
    height: 300,
    dpr: 2,
    loaders: { load: () => Promise.reject(new Error('none')), stream: s.stream },
    registry: [
      {
        navigators: [ideogramNavigator(bands.chromosomes)],
        presets: [{ name: 'karyo', geometry: 'linear', render: {}, navigator: 'ideogram' }],
      },
    ],
  });
  t.ctx.ui.setStyle(GENOME, 'karyo');
  await t.ctx.graph.update();
  t.vp?.zoomTo(0, 100_000_000, 100_050_000);
  // The nav pass drawn by hand, as the composer does on a layout or skin change (a pan never redraws the strip).
  const { NavView, navDraw, navDrawWindow } = await import('../src/draw/nav.js');
  const av = t.ctx.layout.axis(GENOME);
  const entry = t.ctx.layout.navigatorOf(0);
  if (!av || !entry) throw new Error('no navigator');
  const n = new NavView().from(av);
  const be = recording();
  n.backend = be;
  n.text = t.ctx.text;
  n.parts = entry.parts;
  n.dpr = 2;
  n.top = 0;
  n.height = 24;
  n.gripTop = 0;
  n.gripBottom = 24;
  n.names = av.axis?.names ?? [];
  n.from(av);
  let i = 0;
  const at = Float64Array.of(av.navX0, av.navX0 + 9, av.navX0 + 40, av.navX0 + 49);
  // The window moves between two places over the strip; the strip is redrawn each time.
  const frame = () => {
    n.x0 = at[i & 1];
    n.x1 = at[2 + (i++ & 1)];
    be.reset();
    navDraw(entry, n);
    navDrawWindow(entry, n);
  };
  for (let k = 0; k < 8; k++) frame();
  return { frame, backend: be, dispose: () => t.ctx.dispose() };
}
