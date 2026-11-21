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
  it('a button a chip with its words (a mark none), a press shows the hidden tracks in one patch', async () => {
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
    const group = (track?: string) =>
      container.querySelector<HTMLElement>(
        track ? `[role="group"][data-track="/views/main/protein/${track}"]` : '[role="group"]:not([data-track])',
      )!;
    // the section's "2 hidden" is a chip in words; Domains' two states are marks (the column holds one slot)
    expect(group()!.textContent).toBe('2 hidden');
    const marks = [...group('domains').querySelectorAll('button')];
    expect(marks.map((b) => b.dataset.part)).toEqual(['seq:filtered', 'seq:hidden-layers']);
    expect(marks.map((b) => b.textContent)).toEqual(['', '']);
    expect(marks.map((b) => b.getAttribute('aria-label'))).toEqual([
      '1 filtered: Clear the filter',
      '1 layer hidden: Show every layer',
    ]);
    // one patch: the press runs the chip's entry
    group()!.querySelector<HTMLElement>('button')!.click();
    await v.settle();
    expect(v.view.resolved.value!.sections[0].tracks).toHaveLength(5);
    port.dispose?.();
  });
});
