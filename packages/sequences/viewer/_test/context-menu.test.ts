/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The `contextMenu` behavior (seqspec:viewer#behaviors): a secondary-button press on a track asks for its menu
// (`events.menu`, origin `pointer`, with the pointer's place); the menu is the host's. Input delivers, the behavior reacts.

import { describe, expect, it } from 'vitest';
import { case3 } from '../../core/_test/format/cases.js';
import type { ActivateEvent } from '../src/index.js';
import { driven } from './input-fixtures.js';

const S = 'family.main';

describe('the track menu asked for by a secondary press (events.menu)', () => {
  async function opened() {
    const t = await driven(case3(), { size: { [S]: [1000, 300] } });
    const asked: ActivateEvent[] = [];
    t.ctx.events.menu.subscribe((e) => asked.push({ ...e }));
    return { t, asked };
  }
  /** The middle of the first track that is not the ruler, in surface px. */
  const onTrack = (t: Awaited<ReturnType<typeof driven>>) => {
    const l = t.layout(S);
    const e = l.pinnedCount;
    const y = (e < l.pinnedCount ? l.entryTop[e] : l.scrolledTop + l.entryTop[e]) + l.entryHeight[e] / 2;
    return { x: l.dataLeft + 40, y };
  };

  it("asks for the track under the pointer on a secondary press, with the pointer's place; a primary press does not", async () => {
    const { t, asked } = await opened();
    const p = onTrack(t);
    t.pointer('down', S, p.x, p.y, { button: 2 });
    t.pointer('up', S, p.x, p.y, { button: 2 });
    expect(asked).toHaveLength(1);
    expect(asked[0]).toMatchObject({ section: S, origin: 'pointer', x: p.x, y: p.y });
    expect(asked[0].track).toBeTruthy();
    asked.length = 0;
    t.click(S, p.x, p.y);
    expect(asked).toEqual([]);
    t.ctx.dispose();
  });

  it('asks nothing where no track is (the navigator, past the rows), and is a behavior with a region param', async () => {
    const { t, asked } = await opened();
    expect(t.ctx.behaviors.get('contextMenu')).not.toBeNull();
    t.pointer('down', S, 5, t.layout(S).height + 40, { button: 2 });
    t.pointer('up', S, 5, t.layout(S).height + 40, { button: 2 });
    expect(asked).toEqual([]);
    // Off for the data: the press there asks nothing.
    t.ctx.behaviors.get('contextMenu')?.setParams({ regions: ['header'] });
    const p = onTrack(t);
    t.pointer('down', S, p.x, p.y, { button: 2 });
    t.pointer('up', S, p.x, p.y, { button: 2 });
    expect(asked).toEqual([]);
    t.ctx.dispose();
  });

  it('a spec without the behavior asks nothing', async () => {
    const { t, asked } = await opened();
    t.ctx.behaviors.get('contextMenu')?.remove();
    const p = onTrack(t);
    t.pointer('down', S, p.x, p.y, { button: 2 });
    expect(asked).toEqual([]);
    t.ctx.dispose();
  });
});
