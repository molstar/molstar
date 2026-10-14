/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The `contextMenu` behavior (seqspec:interaction#the-track-menu): a press of its `on` binding (the secondary button)
// on a track asks for its menu (`events.menu`, origin `pointer`, with the pointer's place); the menu is the host's.
// Not in `DefaultSeqViewerSpec`: a host adds it (`SeqViewerOptionalBehaviors`). Input delivers, the behavior reacts.

import { describe, expect, it } from 'vitest';
import { case3 } from '../../core/_test/format/cases.js';
import { Hit, HitRegion, type MenuEvent, SeqViewerOptionalBehaviors, SeqViewerSpec } from '../src/index.js';
import { driven } from './input-fixtures.js';

const S = 'family.main';
type Driven = Awaited<ReturnType<typeof driven>>;

describe('the track menu asked for by a secondary press (events.menu)', () => {
  /** case 3 in 1000 × 300: the navigator on top, the pinned alignment rows, one scrolled track, then empty rows. */
  async function opened(withMenu = true) {
    const t = await driven(case3(), {
      size: { [S]: [1000, 300] },
      behaviors: withMenu
        ? (list) => [...list, SeqViewerSpec.Behavior(SeqViewerOptionalBehaviors.ContextMenu)]
        : undefined,
    });
    const asked: MenuEvent[] = [];
    t.ctx.events.menu.subscribe((e) => asked.push({ ...e }));
    return { t, asked };
  }
  /** The middle of the first scrolled track's data, in surface px. */
  const onTrack = (t: Driven) => {
    const l = t.layout(S);
    const e = l.pinnedCount;
    return { x: l.dataLeft + 40, y: l.scrolledTop + l.entryTop[e] + l.entryHeight[e] / 2 };
  };
  /** The first point down the column at `x` whose hit is in `region` with (or without) a track. */
  const find = (t: Driven, x: number, region: number, track: boolean) => {
    const h = Hit.create();
    for (let y = 1; y < t.layout(S).height; y++) {
      t.ctx.layout.hitTest(S, x, y, h);
      if (h.region === region && h.track >= 0 === track) return { x, y };
    }
    throw new Error(`no ${track ? 'track' : 'empty'} point in region ${region}`);
  };
  const press = (t: Driven, p: { x: number; y: number }, opts: { button?: number; ctrl?: boolean } = {}) => {
    t.pointer('down', S, p.x, p.y, { button: 2, ...opts });
    t.pointer('up', S, p.x, p.y, { button: 2, ...opts });
  };

  it('is not in DefaultSeqViewerSpec: a secondary press asks nothing until a host adds the behavior', async () => {
    const { t, asked } = await opened(false);
    expect(t.ctx.behaviors.get('contextMenu')).toBeNull();
    press(t, onTrack(t));
    expect(asked).toEqual([]);
    t.ctx.dispose();
  });

  it("asks for the track under the pointer on a secondary press, with the pointer's place; a primary press does not", async () => {
    const { t, asked } = await opened();
    const p = onTrack(t);
    t.pointer('down', S, p.x, p.y, { button: 2 });
    // On the press, before the release.
    expect(asked).toHaveLength(1);
    t.pointer('up', S, p.x, p.y, { button: 2 });
    expect(asked).toHaveLength(1);
    expect(asked[0]).toMatchObject({ section: S, origin: 'pointer', x: p.x, y: p.y });
    expect(asked[0].track).toBeTruthy();
    asked.length = 0;
    t.click(S, p.x, p.y);
    expect(asked).toEqual([]);
    t.ctx.dispose();
  });

  it('asks nothing where no track is: on the navigator, or on the data area past the last row', async () => {
    const { t, asked } = await opened();
    const x = t.layout(S).dataLeft + 40;
    press(t, find(t, x, HitRegion.navigator, false));
    press(t, find(t, x, HitRegion.data, false));
    expect(asked).toEqual([]);
    t.ctx.dispose();
  });

  it('asks only in its `regions`: on the header with `["header"]`, not on the data', async () => {
    const { t, asked } = await opened();
    t.ctx.behaviors.get('contextMenu')?.setParams({ regions: ['header'] });
    press(t, onTrack(t));
    expect(asked).toEqual([]);
    const header = find(t, 20, HitRegion.header, true);
    press(t, header);
    expect(asked).toHaveLength(1);
    expect(asked[0]).toMatchObject({ section: S, origin: 'pointer', x: header.x, y: header.y });
    t.ctx.dispose();
  });

  it('asks on the presses of its `on` binding: a host adds Ctrl+primary (macOS), a plain primary press still does not', async () => {
    const { t, asked } = await opened();
    const p = onTrack(t);
    press(t, p, { button: 0, ctrl: true });
    expect(asked).toEqual([]);
    t.ctx.behaviors.get('contextMenu')?.setParams({ on: [{ button: 'secondary' }, { button: 'primary', ctrl: true }] });
    press(t, p, { button: 0, ctrl: true });
    expect(asked).toHaveLength(1);
    press(t, p, { button: 0 });
    expect(asked).toHaveLength(1);
    press(t, p);
    expect(asked).toHaveLength(2);
    t.ctx.dispose();
  });

  it('a spec without the behavior asks nothing (removed at runtime)', async () => {
    const { t, asked } = await opened();
    t.ctx.behaviors.get('contextMenu')?.remove();
    press(t, onTrack(t));
    expect(asked).toEqual([]);
    t.ctx.dispose();
  });
});
