/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { DefaultSeqViewerSpec } from '@molstar/seqx-viewer/default-spec';
import { createSeqViewer } from '@molstar/seqx-viewer/viewer';
import { defineElements, onCleanup } from '@molstar/ui/dom';
import { UiElement } from '@molstar/ui/elements/ui-element';
import { followSvg } from '../../follow-svg.js';
import { type CapEntry, Caps, capConfig, capSpec } from './compose.js';

// The Layouts story's element: a tracks section past its cap at each cap, in the light theme and the dark one, each a
// live canvas (wheel, rail, Page and Home/End keys scroll the tracks below the fixed ruler) beside its SVG export (every
// track: the figure). The theme is the `data-msui-theme` of the cell's box.

const Muted = 'margin: 0; color: var(--msui-text-muted); max-width: 72ch;';

function text(tag: string, content: string, css = ''): HTMLElement {
  const el = document.createElement(tag);
  el.style.cssText = css;
  el.textContent = content;
  return el;
}

/** A canvas of the cap and its SVG; retired with the story. */
function mount(cap: CapEntry, canvas: HTMLElement, svg: HTMLElement): void {
  const base = DefaultSeqViewerSpec();
  const viewer = createSeqViewer({ ...base, config: [...(base.config ?? []), ...capConfig(cap)] }, canvas);
  const followed = followSvg(viewer, svg, { prefix: `cap-${cap.id}` });
  onCleanup(() => {
    followed.dispose();
    viewer.retire();
  });
  void viewer.load(capSpec(cap)).then(() => viewer.settle());
}

function theme(name: 'light' | 'dark', cap: CapEntry): HTMLElement {
  const box = document.createElement('div');
  box.setAttribute('data-msui-theme', name);
  box.style.cssText =
    'display: grid; gap: var(--msui-space-3); align-content: start; min-width: 0; padding: var(--msui-space-3); background: var(--msui-bg); color: var(--msui-text); border: 1px solid var(--msui-line); border-radius: 8px;';
  box.append(text('h5', name === 'light' ? 'Light' : 'Dark', Muted));
  const canvas = document.createElement('div');
  canvas.style.cssText = 'min-width: 0;';
  canvas.setAttribute('data-cap', cap.id);
  const svg = document.createElement('div');
  svg.style.cssText = 'min-width: 0; max-height: 360px; overflow: auto; border-top: 1px solid var(--msui-line);';
  box.append(text('h6', 'Canvas', Muted), canvas, text('h6', 'SVG: every track', Muted), svg);
  mount(cap, canvas, svg);
  return box;
}

/** Each cap, light and dark, a canvas beside its SVG. */
export class LayoutsRoot extends UiElement<readonly CapEntry[]> {
  static readonly tag = 'seqx-layouts-cap-root';

  protected render(caps: readonly CapEntry[]): HTMLElement {
    const root = document.createElement('div');
    root.style.cssText =
      'display: grid; gap: var(--msui-space-6); align-content: start; padding: var(--msui-space-5); min-width: 0; background: var(--msui-bg); color: var(--msui-text); font: var(--msui-text-sm) / var(--msui-leading-sm) var(--msui-font-sans);';
    root.append(
      text(
        'p',
        'A tracks section taller than its box scrolls inside it: the ruler stays, the tracks below scroll by the wheel (the lanes band first, then the section, then the page), the one rail at the right edge, and Page Up, Page Down, Home and End with the pointer in it. An export shows every track.',
        Muted,
      ),
    );
    for (const cap of caps) {
      const pair = document.createElement('div');
      pair.style.cssText =
        'display: grid; grid-template-columns: repeat(auto-fit, minmax(min(480px, 100%), 1fr)); gap: var(--msui-space-4); min-width: 0;';
      pair.append(theme('light', cap), theme('dark', cap));
      root.append(text('h4', cap.title, 'margin: 0;'), pair);
    }
    return root;
  }
}

defineElements([LayoutsRoot]);

export function tracksPastTheCapStory(): HTMLElement {
  const el = document.createElement(LayoutsRoot.tag) as LayoutsRoot;
  el.model = Caps;
  return el;
}
