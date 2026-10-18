/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The case of the header hover's guard (seqspec:interaction#hover; header-hover.heap.ts) and its quick twin
// (header-hover.quick.ts): the hover alternating between two track headers and a group's header (forty, ten and five
// items), each landing described by the `headerHover` behavior and projected by the channels. Each call loads the viewer
// afresh (`vi.resetModules`), as the draw guards do. The heap guard and the quick twin share it.

import { createSeqViewSpec } from '@molstar/seq-core/builder';
import { vi } from 'vitest';

/** The heap guard's steps in the shorter run; the longer runs four times as many. */
export const N = 50_000;
const S = 'v.main';

function doc() {
  const b = createSeqViewSpec({ title: 'Headers' });
  const a = b.data.sequence({ ref: 'a', value: 'MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ'.repeat(4) });
  const ranges = (ref: string, n: number, step: number) => {
    const start = Array.from({ length: n }, (_, i) => step * i + 1);
    return b.data.annotation({
      ref,
      geometry: 'ranges',
      target: a,
      value: { start, stop: start.map((s) => s + 1), label: start.map((s) => `${ref}${s}`) },
    });
  };
  const v = b.view({ ref: 'v', label: 'V' }).block({ sequence: a });
  v.track({ ref: 'v.forty', label: 'Forty' }).ranges({ data: ranges('d40', 40, 3) });
  v.track({ ref: 'v.ten', label: 'Ten' }).ranges({ data: ranges('d10', 10, 7) });
  v.group({ ref: 'v.g', label: 'Group', collapsed: false } as never)
    .track({ ref: 'v.g.five', label: 'Five' })
    .ranges({ data: ranges('d5', 5, 11) });
  return b.getDocument();
}

/** The headers' document opened afresh, and a step: the hover on a header by turns, described and projected. */
export async function headers() {
  vi.resetModules();
  const { opened } = await import('./layout-fixtures.js');
  const { ChannelMarks, EntryKind, Hit } = await import('../src/index.js');
  const t = await opened(doc(), { size: { [S]: [1000, 500] } });
  t.ctx.behaviors.get('hoverLabel')?.remove();
  const l = t.layout(S);
  const at = (y: number) => {
    const h = Hit.create();
    t.ctx.layout.hitTest(S, 10.5, y, h);
    return h;
  };
  let groupY = -1;
  for (let e = 0; e < l.entryCount; e++)
    if (l.entryKind[e] === EntryKind.group) groupY = l.scrolledTop + l.entryTop[e] + l.entryHeight[e] / 2;
  const hits = [at(t.yOf(S, 'v.forty')), at(t.yOf(S, 'v.ten')), at(groupY)];
  const ch = new ChannelMarks(t.ctx);
  const header = t.ctx.hover.live.current.header;
  const counts = [0, 0, 0];
  const out = { t, ch, counts, sink: 0, step: (_: number) => {} };
  out.step = (i: number) => {
    const k = i % 3;
    t.ctx.hover.setFromHit(hits[k], 'pointer');
    ch.sync();
    counts[k] = header.items;
    out.sink += ch.marks(0)?.hover.count ?? 0;
  };
  return out;
}
