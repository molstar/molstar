// @vitest-environment jsdom
/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The headers' DOM (_spec/README.md "Ports", Header menu and Header parts) in jsdom: the menu's overlay opened by a
// secondary click (the browser's menu prevented, Shift's left alone), its focus in, Esc, the focus back to the canvas;
// `headerPartsDom()` placing a part at each box of a header whose actions show, a press running its entry (one patch),
// roving focus, Esc back to the canvas and Tab never captured (EC-155).

import { createSeqSpec } from '@molstar/seqx-core/lang';
import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { type HeadersBehavior, HeadersProvider } from '../behaviors/headers.js';
import { RecordingContext } from '../canvas/_test/recording.js';
import { DefaultSeqViewerSpec } from '../default-spec.js';
import { headerMenuOverlay } from '../mount-header-menu.js';
import { headerPartsDom } from '../mount-header-parts.js';
import { type UiHistory, uiHistory } from '../ui-history.js';
import { createSeqViewer, type SeqViewer } from '../viewer.js';
import { type Frames, installFrames } from './frames.js';

const Section = '/views/main/protein';
const Domains = `${Section}/domains`;

function spec() {
  const b = createSeqSpec({ title: 'header overlay' });
  const p = b.data.sequence('p', { alphabet: 'protein', letters: 'MEEPQSDPSVEPPLSQETFSDLWKLLPENN' });
  b.views.view('main', (v) => {
    const s = v.section('protein', { of: p });
    s.track('letters', { label: 'Letters' }).layer('l', { type: 'letters' });
    s.track('domains', { label: 'Domains' }).layer('l', { type: 'letters' });
  });
  return b.build();
}

let frames: Frames;
const viewers: SeqViewer[] = [];
beforeEach(() => {
  frames = installFrames();
  vi.stubGlobal(
    'ResizeObserver',
    class {
      observe() {}
      unobserve() {}
      disconnect() {}
    },
  );
  vi.spyOn(HTMLElement.prototype, 'clientWidth', 'get').mockReturnValue(600);
  vi.spyOn(HTMLCanvasElement.prototype, 'getContext').mockImplementation(function (this: HTMLCanvasElement) {
    return ((this as any).__ctx ??= new RecordingContext());
  });
});
afterEach(() => {
  for (const v of viewers.splice(0)) v.retire();
  frames.restore();
  document.body.textContent = '';
});

async function mount(dom: boolean) {
  const target = document.createElement('div');
  document.body.appendChild(target);
  const base = DefaultSeqViewerSpec();
  const v = createSeqViewer(
    {
      ...base,
      ui: { canvas: true, history: uiHistory },
      behaviors: base.behaviors.map((b) =>
        b.provider === HeadersProvider
          ? {
              provider: HeadersProvider,
              params: { menu: headerMenuOverlay(), ...(dom ? { parts: headerPartsDom() } : {}) },
            }
          : b,
      ),
    },
    target,
  );
  viewers.push(v);
  await v.load(spec());
  await v.settle();
  frames.flush();
  const container = target.querySelector<HTMLElement>('[role="application"]')!;
  return { v, container, headers: v.behaviors.get('headers') as HeadersBehavior };
}

const geo = (v: SeqViewer, track: string) => v.headers!.value.find((g) => g.target.track === track)!;
const mouse = (type: string, x: number, y: number, init: MouseEventInit = {}) =>
  new MouseEvent(type, { clientX: x, clientY: y, bubbles: true, cancelable: true, ...init });

describe('the header menu’s overlay', () => {
  it('a secondary click opens it (the browser’s prevented), the focus in; Esc closes it, the focus back', async () => {
    const { v, container, headers } = await mount(false);
    container.focus();
    const g = geo(v, Domains);
    const at = [g.label.x + 2, g.line.y + g.line.height / 2] as const;
    // Shift: the browser's, whatever the param
    const shifted = mouse('contextmenu', ...at, { shiftKey: true });
    container.dispatchEvent(shifted);
    expect(shifted.defaultPrevented).toBe(false);
    expect(container.querySelector('[role="menu"]')).toBeNull();
    const e = mouse('contextmenu', ...at);
    container.dispatchEvent(e);
    expect(e.defaultPrevented).toBe(true);
    const menu = container.querySelector('[role="menu"]')!;
    const labels = [...menu.querySelectorAll('[role^="menuitem"]')].map((b) => b.textContent);
    // the track's items, then the section's (a key quiet at the row's end)
    expect(labels.map((l) => l!.replace(/Shift\+\w$/, ''))).toEqual([
      'Hide track',
      'Color ▸',
      'Move up',
      'Reverse',
      'Collapse or expand hidden stretches',
    ]);
    expect(menu.contains(document.activeElement)).toBe(true);
    document.activeElement!.dispatchEvent(new KeyboardEvent('keydown', { key: 'Escape', bubbles: true }));
    expect(container.querySelector('[role="menu"]')).toBeNull();
    expect(headers.menu.value).toBeUndefined();
    expect(document.activeElement).toBe(container);
  });
});

describe('headerPartsDom()', () => {
  it('a part at each box where the actions show; a press runs its entry (one patch); roving focus, Esc back, Tab never captured (EC-155)', async () => {
    const { v, container } = await mount(true);
    const layer = container.querySelector<HTMLElement>('.seqx-header-parts')!;
    // at rest: no action parts (a 200-track view holds a few buttons, not 600)
    expect(layer.querySelectorAll('[data-part="seq:hide"]')).toHaveLength(0);
    const g = geo(v, Domains);
    container.dispatchEvent(
      new PointerEvent('pointermove', { clientX: g.label.x + 2, clientY: g.line.y + 4, bubbles: true }),
    );
    frames.flush();
    const hide = layer.querySelector<HTMLElement>('[data-part="seq:hide"]')!;
    const box = geo(v, Domains).slots.find((s) => s.name === 'seq:hide')!.box;
    expect([hide.style.left, hide.style.top]).toEqual([`${box.x}px`, `${box.y}px`]);
    expect(layer.querySelectorAll('[data-part="seq:hide"]')).toHaveLength(1);
    // roving: one tab stop, → to the next part, Esc back to the canvas; Tab is the page's
    const parts = [...hide.parentElement!.querySelectorAll<HTMLElement>('button')];
    expect(parts.map((b) => b.tabIndex)).toEqual([0, ...parts.slice(1).map(() => -1)]);
    hide.focus();
    hide.dispatchEvent(new KeyboardEvent('keydown', { key: 'ArrowRight', bubbles: true }));
    expect(document.activeElement).toBe(parts[1]);
    const tab = new KeyboardEvent('keydown', { key: 'Tab', bubbles: true, cancelable: true });
    parts[1].dispatchEvent(tab);
    expect(tab.defaultPrevented).toBe(false);
    parts[1].dispatchEvent(new KeyboardEvent('keydown', { key: 'Escape', bubbles: true }));
    expect(document.activeElement).toBe(container);
    // a press: the slot's entry, one patch
    const steps = (v.spec.history as UiHistory).stack.steps.value.length;
    layer.querySelector<HTMLElement>('[data-part="seq:hide"]')!.click();
    await new Promise((done) => setTimeout(done, 0));
    frames.flush();
    expect(v.view.resolved.value!.sections[0].tracks.map((t) => t.address)).toEqual([`${Section}/letters`]);
    expect((v.spec.history as UiHistory).stack.steps.value.length).toBe(steps + 1);
  });
});
