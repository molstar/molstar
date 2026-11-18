/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { ExampleView } from '@molstar/seqx-examples';
import { Button } from '@molstar/ui/components/button/button';
import { defineElements, onCleanup } from '@molstar/ui/dom';
import { UiElement } from '@molstar/ui/elements/ui-element';
import { followSvg } from '../../follow-svg.js';
import { composeViewer, describeSelection } from './compose.js';

// The selection story's root, shared by the Features story and the p53 example (each gives its own hint): the default
// grammar on a view (the interaction grammar). Click a residue or a domain to select it (again to clear), Shift+click to
// extend from it, Cmd/Ctrl+click to toggle, drag on the letters to select a range, click empty space or press Esc to
// clear. The arrow keys move the keyboard cursor (Shift extends; Space toggles its cell; Enter focuses). On the protein
// on its gene, what is selected or hovered on one section shows on the other in the same look, lighter. The SVG beside
// the canvas is an export of the same spec with the same selection and highlights; the hover and the cursor are
// transient and never in an export. The toggle draws the selection as one band a section (`column`, spanning the
// tracks and the gaps between them) or a box a track (`tracks`).

const Box = 'display: grid; gap: var(--msui-space-2); min-width: 0; align-content: start;';
const Title =
  'margin: 0; font-size: var(--msui-text-sm); font-weight: var(--msui-weight-strong); color: var(--msui-text-muted);';
const Mono =
  'margin: 0; font: var(--msui-text-xs) / var(--msui-leading-xs) var(--msui-font-mono); color: var(--msui-text-muted);';

/** The default hint: every gesture. */
export const GestureHint =
  'Click selects (again: clears), Shift+click extends, Cmd/Ctrl+click toggles, a drag on the data selects a range, Esc or a click on empty space clears. Arrows move the cursor (Shift extends, Space toggles, Enter focuses). The hover and the cursor are transient: the SVG export shows the selection and the highlights alone.';

function section(title: string, body: HTMLElement): HTMLElement {
  const el = document.createElement('div');
  el.style.cssText = Box;
  const h = document.createElement('h4');
  h.style.cssText = Title;
  h.textContent = title;
  el.append(h, body);
  return el;
}

/** A story's model: the view, how the selection is drawn, and what to try. */
export interface SelectionStory {
  readonly view: ExampleView;
  readonly shape: 'column' | 'tracks';
  readonly hint?: string;
}

class SelectionRoot extends UiElement<SelectionStory> {
  static readonly tag = 'seqx-selection-root';

  protected render({ view, shape, hint: text }: SelectionStory): HTMLElement {
    const root = document.createElement('div');
    root.style.cssText =
      'display: grid; align-content: start; gap: var(--msui-space-4); padding: var(--msui-space-5); min-height: 100vh; background: var(--msui-bg); color: var(--msui-text); font: var(--msui-text-sm) / var(--msui-leading-sm) var(--msui-font-sans);';
    const hint = document.createElement('p');
    hint.style.cssText = 'margin: 0; color: var(--msui-text-muted); max-width: 72ch;';
    hint.textContent = text ?? GestureHint;
    const tools = document.createElement('div');
    tools.style.cssText = 'display: flex; gap: var(--msui-space-2); align-items: baseline;';
    const toggle = Button({
      children: shape === 'column' ? 'Shape: column (switch to tracks)' : 'Shape: tracks (switch to column)',
      onClick: () => {
        this.model = { view, shape: shape === 'column' ? 'tracks' : 'column', ...(text ? { hint: text } : {}) };
      },
    }) as HTMLElement;
    tools.append(toggle);
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
    root.append(hint, tools, pair, status);

    const viewer = composeViewer(view, { target: canvasHost, look: { shape } });
    onCleanup(() => viewer.retire());

    // the SVG follows the selection (once a drag-select ends, not each move), the highlights, the theme and the view at
    // rest; a later draw drops an earlier one
    const svg = followSvg(viewer, svgHost, { prefix: 'sel', also: [viewer.highlight.version] });
    onCleanup(() => svg.dispose());
    const showStatus = () => {
      status.textContent = describeSelection(viewer);
    };
    let dragged = false;
    const offs = [
      viewer.selection.value.listen(() => {
        showStatus();
        if (viewer.input.dragging) dragged = true;
        else void svg.draw();
      }),
      viewer.selection.cursor.listen(showStatus),
    ];
    // a drag-select's selection is exported once, when the button goes up
    const onUp = () => {
      if (!dragged) return;
      dragged = false;
      void svg.draw();
    };
    canvasHost.addEventListener('pointerup', onUp);
    canvasHost.addEventListener('lostpointercapture', onUp);
    onCleanup(() => {
      for (const off of offs) off();
      canvasHost.removeEventListener('pointerup', onUp);
      canvasHost.removeEventListener('lostpointercapture', onUp);
    });
    void viewer
      .load(view.build())
      .then(() => viewer.settle())
      .then(() => {
        showStatus();
        return svg.draw();
      })
      .catch((e: unknown) => {
        status.textContent = `The view did not load: ${e instanceof Error ? e.message : String(e)}`;
      });
    return root;
  }
}

defineElements([SelectionRoot]);

/** A Storybook `render` result: the view, its selection drawn as a column band, and the hint to show. */
export function selectionStory(view: ExampleView, hint?: string): HTMLElement {
  const el = document.createElement(SelectionRoot.tag) as SelectionRoot;
  el.model = { view, shape: 'column', ...(hint ? { hint } : {}) };
  return el;
}
