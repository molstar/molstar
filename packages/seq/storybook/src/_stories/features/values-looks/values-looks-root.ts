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
import { type CellEntry, Cells, lookConfig, type LookEntry, Looks } from './compose.js';

// The values looks story's element: each look (today's and the three variants) in the light theme and the dark one side by
// side, each a viewer of the examples' views per cell, live (hover, drag and the keys work in each). The theme is the
// `data-msui-theme` of the cell's box: the viewer reads its tokens from the nearest one.

const Muted = 'margin: 0; color: var(--msui-text-muted); max-width: 72ch;';

function text(tag: string, content: string, css = ''): HTMLElement {
  const el = document.createElement(tag);
  el.style.cssText = css;
  el.textContent = content;
  return el;
}

/** A viewer of a cell's view with the look, in `host`; retired with the story. */
function mount(cell: CellEntry, look: LookEntry, host: HTMLElement): void {
  const base = DefaultSeqViewerSpec();
  const extra = cell.view.registry?.() ?? [];
  const viewer = createSeqViewer(
    {
      ...base,
      registry: withViewRegistry(base.registry ?? [], extra),
      config: [...(base.config ?? []), lookConfig(look.look)],
    },
    host,
  );
  onCleanup(() => viewer.retire());
  void viewer
    .load(cell.view.build())
    .then(() => viewer.settle())
    .then(() => cell.opening?.(viewer));
}

function theme(name: 'light' | 'dark', look: LookEntry): HTMLElement {
  const box = document.createElement('div');
  box.setAttribute('data-msui-theme', name);
  box.style.cssText =
    'display: grid; gap: var(--msui-space-3); align-content: start; min-width: 0; padding: var(--msui-space-3); background: var(--msui-bg); color: var(--msui-text); border: 1px solid var(--msui-line); border-radius: 8px;';
  box.append(text('h5', name === 'light' ? 'Light' : 'Dark', Muted));
  for (const cell of Cells) {
    const host = document.createElement('div');
    host.style.cssText = 'min-width: 0; max-height: 300px; overflow: hidden;';
    host.setAttribute('data-cell', cell.id);
    box.append(text('h6', cell.title, `${Muted} font-size: var(--msui-text-xs);`), host);
    mount(cell, look, host);
  }
  return box;
}

/** Every look beside today's, light and dark. */
export class ValuesLooksRoot extends UiElement<readonly LookEntry[]> {
  static readonly tag = 'seqx-values-looks-root';

  protected render(looks: readonly LookEntry[]): HTMLElement {
    const root = document.createElement('div');
    root.style.cssText =
      'display: grid; gap: var(--msui-space-6); align-content: start; padding: var(--msui-space-5); min-width: 0; background: var(--msui-bg); color: var(--msui-text); font: var(--msui-text-sm) / var(--msui-leading-sm) var(--msui-font-sans);';
    root.append(
      text(
        'p',
        'The values layer in today’s look and three variants, drawn by the one values path: each live, so hover, zoom and pan work in each.',
        Muted,
      ),
    );
    for (const look of looks) {
      const section = document.createElement('section');
      section.style.cssText = 'display: grid; gap: var(--msui-space-3); min-width: 0;';
      section.setAttribute('data-look', look.look);
      const pair = document.createElement('div');
      pair.style.cssText =
        'display: grid; grid-template-columns: repeat(auto-fit, minmax(min(420px, 100%), 1fr)); gap: var(--msui-space-4); min-width: 0;';
      pair.append(theme('light', look), theme('dark', look));
      section.append(text('h4', look.title, 'margin: 0;'), text('p', look.says, Muted), pair);
      root.append(section);
    }
    return root;
  }
}

defineElements([ValuesLooksRoot]);

export function valuesLooksStory(): HTMLElement {
  const el = document.createElement(ValuesLooksRoot.tag) as ValuesLooksRoot;
  el.model = Looks;
  return el;
}
