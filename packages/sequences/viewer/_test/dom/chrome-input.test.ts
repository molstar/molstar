/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The chrome's input through the public streams (seqspec:interaction#streams), in Chromium over the DOM UI: a click, a
// context-menu request, a hover and a leave on a part (its root's `data-seqv-part`, its class the fallback) report its
// kind, section, track and key, origin `chrome`; the
// surfaces are the engine's; nothing is reported to a stream with no subscriber; a host's node tagged `data-seqv-part`
// reports too.

import { afterEach, describe, expect, it } from 'vitest';
import * as D from '../../src/dom/index.js';
import '../../src/skin.css';
import type { SeqClickEvent, SeqInputEvent } from '../../src/index.js';
import { flush, openedContext, stubSurface } from './fixtures.js';

const cleanups: (() => void)[] = [];
afterEach(() => {
  for (const c of cleanups.splice(0)) c();
});
async function mounted() {
  const t = await openedContext();
  const host = document.createElement('div');
  host.style.cssText = 'position:absolute;left:0;top:0;width:900px';
  document.body.appendChild(host);
  const m = D.mountSeqViewer(host, t.ctx, { surface: stubSurface });
  cleanups.push(() => {
    m.dispose();
    host.remove();
    t.ctx.dispose();
  });
  await flush();
  return { ...t, host };
}
const mouse = (type: string, o: MouseEventInit = {}) =>
  new MouseEvent(type, { bubbles: true, cancelable: true, composed: true, ...o });

describe('the chrome through ctx.interaction', () => {
  it('every part root carries its kind in `data-seqv-part`, which is what is reported', async () => {
    const { ctx, host } = await mounted();
    const kinds = new Set([...host.querySelectorAll('[data-seqv-part]')].map((e) => e.getAttribute('data-seqv-part')));
    for (const k of ['viewTab', 'sectionHeader', 'legend', 'snapshotPlayer', 'description'])
      expect(kinds.has(k), k).toBe(true);
    const seen: string[] = [];
    ctx.interaction.click.subscribe((e) => seen.push(e.part.kind));
    // The tag decides, whatever the class.
    const legend = host.querySelector<HTMLElement>('[data-seqv-part="legend"]');
    legend?.setAttribute('data-seqv-part', 'description');
    legend?.querySelector('summary')?.dispatchEvent(mouse('click', { detail: 1 }));
    // Without a tag, the class is the fallback.
    legend?.removeAttribute('data-seqv-part');
    legend?.querySelector('summary')?.dispatchEvent(mouse('click', { detail: 1 }));
    expect(seen).toEqual(['description', 'legend']);
  });

  it("a click on a view tab: its part and key, origin chrome; a click on the surface is not the chrome's", async () => {
    const { ctx, host } = await mounted();
    const seen: SeqClickEvent[] = [];
    ctx.interaction.click.subscribe((e) => seen.push(e));
    const tab = host.querySelector<HTMLElement>('[data-seqv-key="tab:msa"]');
    expect(tab).not.toBeNull();
    tab?.dispatchEvent(mouse('click', { detail: 1, shiftKey: true, clientX: 5, clientY: 6 }));
    host.querySelector('.seqv-surface')?.dispatchEvent(mouse('click', { detail: 1 }));
    expect(seen).toHaveLength(1);
    expect(seen[0]).toMatchObject({
      origin: 'chrome',
      at: null,
      count: 1,
      part: { kind: 'viewTab', key: 'tab:msa', section: null, track: null },
    });
    expect(seen[0].modifiers.shift).toBe(true);
    expect(seen[0].place.clientX).toBe(5);
  });

  it('a section header knows its section; a context-menu request on it is a menu a subscriber may prevent', async () => {
    const { ctx, host } = await mounted();
    const seen: SeqInputEvent[] = [];
    ctx.interaction.menu.subscribe((e) => {
      seen.push(e);
      e.native.preventDefault();
    });
    const header = host.querySelector<HTMLElement>('.seqv-section-header');
    const section = header?.closest('[data-seqv-section]')?.getAttribute('data-seqv-section');
    expect(section).toBeTruthy();
    const e = mouse('contextmenu', { button: 2 });
    header?.dispatchEvent(e);
    expect(seen.map((x) => [x.part.kind, x.part.section])).toEqual([['sectionHeader', section]]);
    expect(e.defaultPrevented).toBe(true);
  });

  it('a hover when the part or its key changes, a leave when the pointer leaves the chrome', async () => {
    const { ctx, host } = await mounted();
    const seen: string[] = [];
    ctx.interaction.hover.subscribe((e) => seen.push(`hover ${e.part.kind} ${e.part.key}`));
    ctx.interaction.leave.subscribe((e) => seen.push(`leave ${e.part.kind}`));
    const tabs = [...host.querySelectorAll<HTMLElement>('[data-seqv-key^="tab:"]')];
    expect(tabs.length).toBe(2);
    const over = (el: Element) => el.dispatchEvent(new PointerEvent('pointerover', { bubbles: true }));
    over(tabs[0]);
    over(tabs[0]);
    over(tabs[1]);
    over(host.querySelector('.seqv-surface') as Element);
    expect(seen).toEqual(['hover viewTab tab:chains', 'hover viewTab tab:msa', 'leave viewTab']);
  });

  it("nothing is reported to a stream with no subscriber; a host's tagged node reports", async () => {
    const { ctx, host } = await mounted();
    const tab = host.querySelector<HTMLElement>('[data-seqv-key="tab:msa"]');
    const before = ctx.views.current.value;
    tab?.dispatchEvent(mouse('click', { detail: 1 }));
    // The tab did its own work (the view switched) and nothing else ran.
    expect(ctx.views.current.value).not.toBe(before);
    const box = document.createElement('div');
    box.setAttribute('data-seqv-part', 'legend');
    const inner = document.createElement('button');
    inner.setAttribute('data-seqv-key', 'mine');
    box.append(inner);
    document.body.append(box);
    const stop = D.listenChrome(box, ctx);
    const seen: SeqClickEvent[] = [];
    ctx.interaction.click.subscribe((e) => seen.push(e));
    inner.dispatchEvent(mouse('click', { detail: 2 }));
    stop();
    inner.dispatchEvent(mouse('click', { detail: 1 }));
    box.remove();
    expect(seen.map((e) => [e.part.kind, e.part.key, e.count])).toEqual([['legend', 'mine', 2]]);
  });
});
