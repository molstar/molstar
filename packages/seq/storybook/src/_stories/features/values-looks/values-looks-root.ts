/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { DefaultSeqViewerSpec } from '@molstar/seqx-viewer/default-spec';
import { createSeqViewer } from '@molstar/seqx-viewer/viewer';
import { defineElements, onCleanup } from '@molstar/ui/dom';
import { UiElement } from '@molstar/ui/elements/ui-element';
import { withViewRegistry } from '../../examples/view-registry.js';
import { type CellEntry, Cells, zoomed } from './compose.js';

// The values look story's element: the look in the light theme and the dark one side by
// side, each a viewer of the examples' views per cell, live (hover, drag and the keys work in each). The theme is the
// `data-msui-theme` of the cell's box: the viewer reads its tokens from the nearest one.

const Muted = 'margin: 0; color: var(--msui-text-muted); max-width: 72ch;';

function text(tag: string, content: string, css = ''): HTMLElement {
  const el = document.createElement(tag);
  el.style.cssText = css;
  el.textContent = content;
  return el;
}

/** A viewer of a cell's view in `host`; retired with the story. */
function mount(cell: CellEntry, host: HTMLElement): void {
  const base = DefaultSeqViewerSpec();
  const extra = cell.view.registry?.() ?? [];
  const viewer = createSeqViewer({ ...base, registry: withViewRegistry(base.registry ?? [], extra) }, host);
  onCleanup(() => viewer.retire());
  void viewer
    .load(cell.view.build())
    .then(() => viewer.settle())
    .then(async () => {
      // the range again until the view holds it (a first fit may land after the viewer's first settle)
      for (let i = 0; i < 20 && !zoomed(cell, viewer); i++) await new Promise((r) => setTimeout(r, 100));
    });
}

function theme(name: 'light' | 'dark', cells: readonly CellEntry[]): HTMLElement {
  const box = document.createElement('div');
  box.setAttribute('data-msui-theme', name);
  box.style.cssText =
    'display: grid; gap: var(--msui-space-3); align-content: start; min-width: 0; padding: var(--msui-space-3); background: var(--msui-bg); color: var(--msui-text); border: 1px solid var(--msui-line); border-radius: 8px;';
  box.append(text('h5', name === 'light' ? 'Light' : 'Dark', Muted));
  for (const cell of cells) {
    const host = document.createElement('div');
    host.style.cssText = 'min-width: 0; max-height: 300px; overflow: hidden;';
    host.setAttribute('data-cell', cell.id);
    box.append(text('h6', cell.title, `${Muted} font-size: var(--msui-text-xs);`), host);
    mount(cell, host);
  }
  return box;
}

/** The look, light and dark. */
export class ValuesLooksRoot extends UiElement<readonly CellEntry[]> {
  static readonly tag = 'seqx-values-looks-root';

  protected render(cells: readonly CellEntry[]): HTMLElement {
    const root = document.createElement('div');
    root.style.cssText =
      'display: grid; gap: var(--msui-space-4); align-content: start; padding: var(--msui-space-5); min-width: 0; background: var(--msui-bg); color: var(--msui-text); font: var(--msui-text-sm) / var(--msui-leading-sm) var(--msui-font-sans);';
    const pair = document.createElement('div');
    pair.style.cssText =
      'display: grid; grid-template-columns: repeat(auto-fit, minmax(min(420px, 100%), 1fr)); gap: var(--msui-space-4); min-width: 0;';
    pair.append(theme('light', cells), theme('dark', cells));
    root.append(
      text(
        'p',
        'The values layer: below 4 px a value one stepped profile (the area at 35% under a 0.75 px line in the full color); from 4 px bars with a 1 px gap and rounded tops. The top tenth of the well is headroom, a hairline marks the baseline, a missing value is a gap. Each is live: hover a value to light it.',
        Muted,
      ),
      pair,
    );
    return root;
  }
}

defineElements([ValuesLooksRoot]);

export function valuesLooksStory(): HTMLElement {
  const el = document.createElement(ValuesLooksRoot.tag) as ValuesLooksRoot;
  el.model = Cells;
  return el;
}
