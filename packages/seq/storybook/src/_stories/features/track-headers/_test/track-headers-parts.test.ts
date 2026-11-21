// @vitest-environment jsdom
/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ExampleViews } from '@molstar/seqx-examples';
import { partsDomPort } from '@molstar/seqx-viewer/header-parts-dom';
import { uiHistory } from '@molstar/seqx-viewer/ui-history';
import { createSeqViewer, type SeqViewer } from '@molstar/seqx-viewer/viewer';
import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { headersOf, headersSpec } from '../compose.js';

// The DOM parts of the Track headers story, in jsdom: the reference port (`partsDomPort`, what `headerPartsDom()`
// loads) over the headless viewer of the story's view, given the behavior's own parts: a button for each chip, the
// words of the geometry's, a press running the entry. The cursor's focus is checked with the canvas mounted
// (viewer's header-overlay.test.ts).

let viewer: SeqViewer | undefined;
// jsdom has no canvas context: its `getContext` logs "Not implemented" on every call, so a mounted viewer gets none
beforeEach(() => {
  vi.spyOn(HTMLCanvasElement.prototype, 'getContext').mockReturnValue(null);
});

afterEach(() => {
  viewer?.retire();
  viewer = undefined;
  document.body.textContent = '';
});

describe('the story’s DOM parts', () => {
  it('a button a chip with its words (a mark none); the section has no part: its hidden tracks are the application’s', async () => {
    const view = ExampleViews.find((v) => v.id === 'headers')!;
    const v = (viewer = createSeqViewer({
      ...headersSpec({ parts: 'dom' }),
      ui: undefined,
      history: uiHistory(),
    }));
    v.resize(1440);
    await v.load(view.build());
    await v.settle();
    const container = document.createElement('div');
    document.body.append(container);
    const port = partsDomPort()(container, headersOf(v));
    port.update(headersOf(v).partsState());
    const group = (track: string) =>
      container.querySelector<HTMLElement>(`[role="group"][data-track="/views/main/protein/${track}"]`)!;
    // no group for the section (nothing in the view for a hidden track); Domains' two states are marks (the column
    // holds one slot)
    expect(container.querySelector('[role="group"]:not([data-track])')).toBeNull();
    expect(container.querySelector('[data-part="seq:hidden-tracks"]')).toBeNull();
    const marks = [...group('domains').querySelectorAll('button')];
    expect(marks.map((b) => b.dataset.part)).toEqual(['seq:filtered', 'seq:hidden-layers']);
    expect(marks.map((b) => b.textContent)).toEqual(['', '']);
    expect(marks.map((b) => b.getAttribute('aria-label'))).toEqual([
      '1 filtered: Clear the filter',
      '1 layer hidden: Show every layer',
    ]);
    port.dispose?.();
  });

  it('an update between pointerdown and click does not lose the click: the pressed button stays until the press ends', async () => {
    vi.useFakeTimers();
    try {
      const view = ExampleViews.find((v) => v.id === 'headers')!;
      const v = (viewer = createSeqViewer({ ...headersSpec({ parts: 'dom' }), ui: undefined, history: uiHistory() }));
      v.resize(1440);
      await v.load(view.build());
      await v.settle();
      const container = document.createElement('div');
      document.body.append(container);
      const port = partsDomPort()(container, headersOf(v));
      const withActions = () => {
        const state = headersOf(v).partsState();
        return { headers: state.headers.map((h) => ({ ...h, actions: true })) };
      };
      port.update(withActions());
      const hide = () =>
        container.querySelector<HTMLElement>('[role="group"][data-track$="/domains"] [data-part="seq:hide"]')!;
      const button = hide();
      const before = v.view.resolved.value!.sections[0].tracks.length;
      button.dispatchEvent(new PointerEvent('pointerdown', { bubbles: true, pointerId: 1 }));
      // the parts update (the pointer moved onto the part: a hover's geometry): the pressed button is still in place
      port.update(withActions());
      expect(button.isConnected).toBe(true);
      expect(hide()).toBe(button);
      button.dispatchEvent(new PointerEvent('pointerup', { bubbles: true, pointerId: 1 }));
      button.click();
      await v.settle();
      expect(v.view.resolved.value!.sections[0].tracks.length).toBe(before - 1);
      // after the press the held update is drawn
      vi.advanceTimersByTime(1);
      expect(hide()).not.toBe(button);
      port.dispose?.();
    } finally {
      vi.useRealTimers();
    }
  });

  it('the actions of a header have a backdrop over the label, with the ellipsis the canvas draws where the label is cut', async () => {
    const view = ExampleViews.find((v) => v.id === 'headers')!;
    const v = (viewer = createSeqViewer({ ...headersSpec({ parts: 'dom' }), ui: undefined, history: uiHistory() }));
    v.resize(1440);
    await v.load(view.build());
    await v.settle();
    const container = document.createElement('div');
    document.body.append(container);
    const port = partsDomPort()(container, headersOf(v));
    const state = headersOf(v).partsState();
    port.update({ headers: state.headers.map((h) => ({ ...h, actions: true })) });
    const tracks = state.headers.filter((h) => h.geometry.target.track && h.geometry.slots.length);
    expect(tracks.length).toBeGreaterThan(0);
    let cuts = 0;
    for (const h of tracks) {
      const g = h.geometry;
      const group = container.querySelector<HTMLElement>(`[role="group"][data-track="${g.target.track}"]`)!;
      const back = group.querySelector<HTMLElement>('[data-backdrop]')!;
      const x0 = Math.min(...g.slots.map((s) => s.box.x)) - 2;
      const cut = g.labelEnd! > x0;
      expect(back.textContent).toBe(cut ? '…' : '');
      // the backdrop reaches the slots' left edge (and the ellipsis's room further left when it cuts)
      expect(parseFloat(back.style.left)).toBeCloseTo(x0 - (cut ? g.ellipsis! : 0), 5);
      if (cut) cuts++;
    }
    // the story's 80 px column cuts at least one label
    expect(cuts).toBeGreaterThan(0);
    port.dispose?.();
  });
});
