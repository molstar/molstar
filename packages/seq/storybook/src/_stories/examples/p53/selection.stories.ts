/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type ExampleView, P53Views } from '@molstar/seqx-examples';
import { defineElements, onCleanup } from '@molstar/ui/dom';
import { UiElement } from '@molstar/ui/elements/ui-element';
import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { composeViewer, describeSelection, selectionSvg } from '../../features/selection/compose.js';

// Selection on the protein on its gene: drag on the gene's letters to select nucleotides; the protein's section shows the
// same selection projected through the codon map, lighter, in the canvas and in the SVG beside it (an export of the
// same spec with the same selection). The hover and the cursor are transient and never in an export.

let instances = 0;

const Box = 'display: grid; gap: var(--msui-space-2); min-width: 0; align-content: start;';
const Title =
  'margin: 0; font-size: var(--msui-text-sm); font-weight: var(--msui-weight-strong); color: var(--msui-text-muted);';

function section(title: string, body: HTMLElement): HTMLElement {
  const el = document.createElement('div');
  el.style.cssText = Box;
  const h = document.createElement('h4');
  h.style.cssText = Title;
  h.textContent = title;
  el.append(h, body);
  return el;
}

class SelectionRoot extends UiElement<ExampleView> {
  static readonly tag = 'seqx-p53-selection-root';

  protected render(view: ExampleView): HTMLElement {
    const root = document.createElement('div');
    root.style.cssText =
      'display: grid; gap: var(--msui-space-4); padding: var(--msui-space-5); min-height: 100vh; background: var(--msui-bg); color: var(--msui-text); font: var(--msui-text-sm) / var(--msui-leading-sm) var(--msui-font-sans);';
    const hint = document.createElement('p');
    hint.style.cssText = 'margin: 0; color: var(--msui-text-muted); max-width: 72ch;';
    hint.textContent =
      'Drag across the gene’s letters (the lower section) to select nucleotides: the protein above shows the codons they fall in. Esc clears.';
    const canvasHost = document.createElement('div');
    canvasHost.style.cssText = 'min-width: 0;';
    const svgHost = document.createElement('div');
    svgHost.style.cssText = 'min-width: 0; overflow-x: auto;';
    const status = document.createElement('p');
    status.style.cssText =
      'margin: 0; color: var(--msui-text-muted); font: var(--msui-text-xs) / var(--msui-leading-xs) var(--msui-font-mono);';
    const pair = document.createElement('div');
    pair.style.cssText =
      'display: grid; grid-template-columns: repeat(auto-fit, minmax(420px, 1fr)); gap: var(--msui-space-5);';
    pair.append(section('Canvas', canvasHost), section('SVG', svgHost));
    root.append(hint, pair, status);

    const viewer = composeViewer(view, { target: canvasHost });
    onCleanup(() => viewer.retire());

    // the SVG follows the selection (once a drag ends, not each move), the highlights, the theme and the view at rest
    const idPrefix = `p53s${++instances}-`;
    let alive = true;
    let latest = 0;
    const drawSvg = async () => {
      const mine = ++latest;
      const r = await selectionSvg(viewer, idPrefix);
      if (!alive || mine !== latest) return;
      // our own output: every string from data or the spec is escaped by the backend
      svgHost.innerHTML = r.svg;
    };
    const showStatus = () => {
      status.textContent = describeSelection(viewer);
    };
    let dragged = false;
    const offs = [
      viewer.selection.value.listen(() => {
        showStatus();
        if (viewer.input.dragging) dragged = true;
        else void drawSvg();
      }),
      viewer.selection.cursor.listen(showStatus),
      viewer.highlight.version.listen(() => void drawSvg()),
      viewer.rested.listen(() => void drawSvg()),
      viewer.tokensVersion.listen(() => void drawSvg()),
    ];
    const onUp = () => {
      if (!dragged) return;
      dragged = false;
      void drawSvg();
    };
    canvasHost.addEventListener('pointerup', onUp);
    canvasHost.addEventListener('lostpointercapture', onUp);
    onCleanup(() => {
      alive = false;
      for (const off of offs) off();
      canvasHost.removeEventListener('pointerup', onUp);
      canvasHost.removeEventListener('lostpointercapture', onUp);
    });
    void viewer
      .load(view.build())
      .then(() => viewer.settle())
      .then(() => {
        showStatus();
        return drawSvg();
      });
    return root;
  }
}

defineElements([SelectionRoot]);

const meta: Meta = { title: 'Examples/p53/Selection' };
export default meta;

/** A drag-select on the gene, projected onto the protein, in the canvas and the SVG. */
export const Selection: StoryObj = {
  render: () => {
    const el = document.createElement(SelectionRoot.tag) as SelectionRoot;
    el.model = P53Views.find((v) => v.id === 'protein-on-gene')!;
    return el;
  },
};
