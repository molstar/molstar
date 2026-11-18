/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type ExampleView, P53Views } from '@molstar/seqx-examples';
import { Button } from '@molstar/ui/components/button/button';
import { defineElements, onCleanup } from '@molstar/ui/dom';
import { UiElement } from '@molstar/ui/elements/ui-element';
import type { Meta, StoryObj } from '@storybook/web-components-vite';
import {
  cancelDrag,
  composeEditing,
  dragHeight,
  endDrag,
  historyRows,
  lastMove,
  printedSpec,
  setDomainsColor,
  toggleTrack,
} from '../../features/editing/compose.js';
import { addTrack } from './compose.js';

// Edit and undo on p53's tracks: add a track, hide one, set a theme; undo and redo through a stand-in host's stack; a drag
// on the slider previews a track's height (Esc cancels, release commits one step). The SVG follows what is shown, and
// the printed spec beside it. Every string from the spec is put in as text, never markup.

let instances = 0;

const Row = 'display: flex; gap: var(--msui-space-2); align-items: center; flex-wrap: wrap;';
const Mono = 'margin: 0; font: var(--msui-text-xs) / var(--msui-leading-xs) var(--msui-font-mono);';

class EditAndUndoRoot extends UiElement<ExampleView> {
  static readonly tag = 'seqx-p53-edit-and-undo-root';

  protected render(view: ExampleView): HTMLElement {
    const root = document.createElement('div');
    root.style.cssText =
      'display: grid; gap: var(--msui-space-4); padding: var(--msui-space-5); min-height: 100vh; background: var(--msui-bg); color: var(--msui-text); font: var(--msui-text-sm) / var(--msui-leading-sm) var(--msui-font-sans);';
    const hint = document.createElement('p');
    hint.style.cssText = 'margin: 0; color: var(--msui-text-muted); max-width: 72ch;';
    hint.textContent =
      'Add a track, hide the conservation, color the domains, then Undo and Redo; each edit is one step. Drag the slider to preview a height (Esc cancels, release commits one step).';
    const canvasHost = document.createElement('div');
    canvasHost.style.cssText = 'min-width: 0;';
    const svgHost = document.createElement('div');
    svgHost.style.cssText = 'min-width: 0; overflow-x: auto;';
    const status = document.createElement('p');
    status.style.cssText = 'margin: 0; color: var(--msui-text-muted);';
    const say = (p: Promise<string> | string) =>
      void Promise.resolve(p).then((m) => {
        status.textContent = m;
      });

    const s = composeEditing(view, { target: canvasHost });
    onCleanup(() => s.dispose());

    const tools = document.createElement('div');
    tools.style.cssText = Row;
    tools.append(
      Button({ size: 'sm', children: 'Add a track', onClick: () => say(addTrack(s)) }) as HTMLElement,
      Button({
        size: 'sm',
        children: 'Hide/show conservation',
        onClick: () => say(toggleTrack(s, 'conservation')),
      }) as HTMLElement,
      Button({ size: 'sm', children: 'Domains red', onClick: () => say(setDomainsColor(s, '#d62728')) }) as HTMLElement,
      Button({
        size: 'sm',
        children: 'Domains blue',
        onClick: () => say(setDomainsColor(s, '#1f77b4')),
      }) as HTMLElement,
      Button({ size: 'sm', variant: 'primary', children: 'Undo', onClick: () => void s.stack.undo() }) as HTMLElement,
      Button({ size: 'sm', children: 'Redo', onClick: () => void s.stack.redo() }) as HTMLElement,
    );

    const drag = document.createElement('label');
    drag.style.cssText = Row;
    const slider = document.createElement('input');
    slider.type = 'range';
    slider.min = '8';
    slider.max = '80';
    slider.value = '24';
    slider.addEventListener('input', () => dragHeight(s, 'domains', Number(slider.value)));
    slider.addEventListener('change', () => void endDrag(s));
    slider.addEventListener('keydown', (e) => {
      if (e.key === 'Escape') say(cancelDrag(s));
    });
    drag.append('Domains height (drag; Esc cancels)', slider);

    const history = document.createElement('pre');
    history.style.cssText = Mono;
    const printed = document.createElement('pre');
    printed.style.cssText = `${Mono} max-height: 320px; overflow: auto; background: var(--msui-bg-subtle); padding: var(--msui-space-3);`;
    const panes = document.createElement('div');
    panes.style.cssText =
      'display: grid; gap: var(--msui-space-4); grid-template-columns: repeat(auto-fit, minmax(320px, 1fr));';
    panes.append(history, printed);
    root.append(hint, canvasHost, tools, drag, status, panes, svgHost);

    // the SVG is the viewer's export of the spec as shown, drawn when the view rests; a later draw drops an earlier one
    const idPrefix = `p53e${++instances}-`;
    let alive = true;
    let latest = 0;
    const drawSvg = async () => {
      const mine = ++latest;
      const r = await s.viewer.renderSvg({ idPrefix });
      if (!alive || mine !== latest) return;
      // our own output: every string from data or the spec is escaped by the backend
      svgHost.innerHTML = r.svg;
    };
    const draw = () => {
      const rows = historyRows(s);
      history.textContent = [
        'Undo (newest first):',
        ...rows.done.map((l) => `  ${l}`),
        'Redo:',
        ...rows.redo.map((l) => `  ${l}`),
      ].join('\n');
      printed.textContent = printedSpec(s.viewer.spec.shown.value);
    };
    for (const state of [s.stack.steps, s.stack.redoable, s.viewer.spec.shown]) onCleanup(state.listen(draw));
    onCleanup(
      s.stack.lastMove.listen(() => {
        status.textContent = lastMove(s);
      }),
    );
    onCleanup(s.viewer.rested.listen(() => void drawSvg()));
    onCleanup(s.viewer.tokensVersion.listen(() => void drawSvg()));
    onCleanup(() => {
      alive = false;
    });
    void s.viewer
      .load(view.build())
      .then(() => s.viewer.settle())
      .then(() => {
        draw();
        return drawSvg();
      });
    return root;
  }
}

defineElements([EditAndUndoRoot]);

const meta: Meta = { title: 'Examples/p53/Edit and undo' };
export default meta;

/** Add a track, hide one, set a theme; Undo and Redo through a host's stack; drag a height to preview it. */
export const EditAndUndo: StoryObj = {
  render: () => {
    const el = document.createElement(EditAndUndoRoot.tag) as EditAndUndoRoot;
    el.model = P53Views.find((v) => v.id === 'tracks')!;
    return el;
  },
};
