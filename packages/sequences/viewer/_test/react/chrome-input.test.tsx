/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The React UI's chrome through the public streams (seqspec:interaction#streams), as the DOM UI's: its root's listener
// reports a click on a part with its kind and key, and a hover on a track's box with its section and track.

import { act } from 'react';
import { describe, expect, it } from 'vitest';
import '../../src/skin.css';
import type { SeqInputEvent } from '../../src/index.js';
import { at, harness, pointer } from './fixtures.js';

const { opened } = harness();
const S = 'chains.main';

describe("the React UI's chrome through ctx.interaction", () => {
  it("a track box's hover and a view tab's click, with their parts", async () => {
    const t = await opened();
    const seen: SeqInputEvent[] = [];
    t.ctx.interaction.click.subscribe((e) => seen.push(e));
    t.ctx.interaction.hover.subscribe((e) => e.origin === 'chrome' && seen.push(e));
    const h = at(t.ctx, S, 't.dom');
    await act(async () => pointer(t.ctx)('move', S, h.header, h.y));
    const box = t.container.querySelector<HTMLElement>('.seqv-track-actions button');
    expect(box).not.toBeNull();
    box?.dispatchEvent(new PointerEvent('pointerover', { bubbles: true }));
    const tab = t.container.querySelector<HTMLElement>('[data-seqv-key="tab:msa"]');
    tab?.dispatchEvent(new MouseEvent('click', { bubbles: true, detail: 1 }));
    expect(tab?.closest('[data-seqv-part]')?.getAttribute('data-seqv-part')).toBe('viewTab');
    expect(box?.closest('[data-seqv-part]')?.getAttribute('data-seqv-part')).toBe('trackAction');
    expect(seen.map((e) => [e.stream, e.part.kind, e.part.key, e.part.section, e.part.track])).toEqual([
      ['hover', 'trackAction', 'track-actions', S, 't.dom'],
      ['click', 'viewTab', 'tab:msa', null, null],
    ]);
  });
});
