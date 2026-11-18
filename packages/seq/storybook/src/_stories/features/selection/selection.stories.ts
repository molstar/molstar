/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type ExampleView, P53Views } from '@molstar/seqx-examples';
import { defineElements, onCleanup } from '@molstar/ui/dom';
import { UiElement } from '@molstar/ui/elements/ui-element';
import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { composeViewer, describeSelection, selectionSvg } from './compose.js';

// Selection: the default grammar on the two p53 views (the interaction grammar). Click a residue or a domain to select
// it (again to clear), Shift+click to extend from it, Cmd/Ctrl+click to toggle, drag on the letters to select a range,
// click empty space or press Esc to clear. The arrow keys move the keyboard cursor (Shift extends; Space toggles its
// cell; Enter focuses). On the protein on its gene, what is selected or hovered on one section shows on the other in
// the same look, lighter. The SVG beside the canvas is an export of the same spec with the same selection.

let instances = 0;

const Box = 'display: grid; gap: var(--msui-space-2); min-width: 0; align-content: start;';
const Title =
  'margin: 0; font-size: var(--msui-text-sm); font-weight: var(--msui-weight-strong); color: var(--msui-text-muted);';
const Mono =
  'margin: 0; font: var(--msui-text-xs) / var(--msui-leading-xs) var(--msui-font-mono); color: var(--msui-text-muted);';

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
  static readonly tag = 'seqx-selection-root';

  protected render(view: ExampleView): HTMLElement {
    const root = document.createElement('div');
    root.style.cssText =
      'display: grid; gap: var(--msui-space-4); padding: var(--msui-space-5); min-height: 100vh; background: var(--msui-bg); color: var(--msui-text); font: var(--msui-text-sm) / var(--msui-leading-sm) var(--msui-font-sans);';
    const hint = document.createElement('p');
    hint.style.cssText = 'margin: 0; color: var(--msui-text-muted); max-width: 72ch;';
    hint.textContent =
      'Click selects (again: clears), Shift+click extends, Cmd/Ctrl+click toggles, a drag on the data selects a range, Esc or a click on empty space clears. Arrows move the cursor (Shift extends, Space toggles, Enter focuses).';
    const canvasHost = document.createElement('div');
    canvasHost.style.cssText = 'min-width: 0;';
    const svgHost = document.createElement('div');
    svgHost.style.cssText = 'min-width: 0; overflow-x: auto;';
    const status = document.createElement('p');
    status.style.cssText = Mono;
    const pair = document.createElement('div');
    pair.style.cssText =
      'display: grid; grid-template-columns: repeat(auto-fit, minmax(420px, 1fr)); gap: var(--msui-space-5);';
    pair.append(section('Canvas', canvasHost), section('SVG', svgHost));
    root.append(hint, pair, status);

    const viewer = composeViewer(view, { target: canvasHost });
    onCleanup(() => viewer.retire());

    // the SVG follows the selection, the highlights, the theme and the view at rest; a later draw drops an earlier one
    const idPrefix = `sel${++instances}-`;
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
    const offs = [
      viewer.selection.value.listen(() => {
        showStatus();
        void drawSvg();
      }),
      viewer.selection.cursor.listen(showStatus),
      viewer.highlight.version.listen(() => void drawSvg()),
      viewer.rested.listen(() => void drawSvg()),
      viewer.tokensVersion.listen(() => void drawSvg()),
    ];
    onCleanup(() => {
      alive = false;
      for (const off of offs) off();
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

const meta: Meta = { title: 'Features/Selection' };
export default meta;

function story(id: string): HTMLElement {
  const el = document.createElement(SelectionRoot.tag) as SelectionRoot;
  el.model = P53Views.find((v) => v.id === id)!;
  return el;
}

/** p53's tracks: select residues, domains and ranges; the cursor; the SVG beside it draws the same selection. */
export const Tracks: StoryObj = { render: () => story('tracks') };

/** The protein on its gene: what is selected or hovered on one section shows on the other, lighter. */
export const ProteinOnItsGene: StoryObj = { render: () => story('protein-on-gene') };
