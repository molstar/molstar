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
  bookmarkView,
  cancelDrag,
  composeEditing,
  dragHeight,
  endDrag,
  historyRows,
  lastMove,
  moveFirst,
  printedSpec,
  recolorStructure,
  setDomainsColor,
  toggleTrack,
  Tracks,
} from './compose.js';

// Editing (Streaming and patches): the p53 tracks view edited through commands, each one patch and one step of a
// stand-in host's stack that also holds a stand-in structure's color. Hide or show a track, move one first, color the
// domains, bookmark the view, recolor the structure (the host's own step); Undo and Redo walk both in order. The slider
// is a drag: each move a preview, its release one step, Esc nothing. The printed spec (views and refs) follows what is
// shown. Every string from the spec is put in as text, never markup.

const Row = 'display: flex; gap: var(--msui-space-2); align-items: center; flex-wrap: wrap;';
const Mono = 'margin: 0; font: var(--msui-text-xs) / var(--msui-leading-xs) var(--msui-font-mono);';

class EditingRoot extends UiElement<ExampleView> {
  static readonly tag = 'seqx-editing-root';

  protected render(view: ExampleView): HTMLElement {
    const root = document.createElement('div');
    root.style.cssText =
      'display: grid; gap: var(--msui-space-4); padding: var(--msui-space-5); min-height: 100vh; background: var(--msui-bg); color: var(--msui-text); font: var(--msui-text-sm) / var(--msui-leading-sm) var(--msui-font-sans);';
    const canvasHost = document.createElement('div');
    canvasHost.style.cssText = 'min-width: 0;';
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
    for (const id of Tracks)
      tools.append(
        Button({ size: 'sm', children: `Hide/show ${id}`, onClick: () => say(toggleTrack(s, id)) }) as HTMLElement,
      );
    tools.append(
      Button({
        size: 'sm',
        children: 'Conservation first',
        onClick: () => say(moveFirst(s, 'conservation')),
      }) as HTMLElement,
      Button({ size: 'sm', children: 'Domains red', onClick: () => say(setDomainsColor(s, '#d62728')) }) as HTMLElement,
      Button({
        size: 'sm',
        children: 'Domains blue',
        onClick: () => say(setDomainsColor(s, '#1f77b4')),
      }) as HTMLElement,
      Button({ size: 'sm', children: 'Bookmark the view', onClick: () => say(bookmarkView(s)) }) as HTMLElement,
    );

    const host = document.createElement('div');
    host.style.cssText = Row;
    const swatch = document.createElement('span');
    swatch.style.cssText = 'display: inline-block; width: 12px; height: 12px; border-radius: 2px;';
    const structureLabel = document.createElement('span');
    host.append(
      swatch,
      structureLabel,
      Button({
        size: 'sm',
        children: 'Recolor the structure',
        onClick: () => {
          void recolorStructure(s, s.structure.value === 'grey' ? 'teal' : 'grey');
        },
      }) as HTMLElement,
      Button({ size: 'sm', variant: 'primary', children: 'Undo', onClick: () => void s.stack.undo() }) as HTMLElement,
      Button({ size: 'sm', children: 'Redo', onClick: () => void s.stack.redo() }) as HTMLElement,
    );

    // a drag: the slider's moves preview, its release commits one step, Esc cancels
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
      'display: grid; gap: var(--msui-space-4); grid-template-columns: minmax(0, 1fr) minmax(0, 1fr);';
    panes.append(history, printed);
    root.append(canvasHost, tools, drag, host, status, panes);

    const draw = () => {
      const rows = historyRows(s);
      history.textContent = [
        'Undo (newest first):',
        ...rows.done.map((l) => `  ${l}`),
        'Redo:',
        ...rows.redo.map((l) => `  ${l}`),
      ].join('\n');
      printed.textContent = printedSpec(s.viewer.spec.shown.value);
      swatch.style.background = s.structure.value;
      structureLabel.textContent = `Structure: ${s.structure.value}`;
    };
    for (const state of [s.stack.steps, s.stack.redoable, s.viewer.spec.shown, s.structure])
      onCleanup(state.listen(draw));
    // the status line follows Undo and Redo too
    onCleanup(
      s.stack.lastMove.listen(() => {
        status.textContent = lastMove(s);
      }),
    );
    void s.viewer.load(view.build()).then(draw);
    return root;
  }
}

defineElements([EditingRoot]);

const meta: Meta = { title: 'Features/Editing' };
export default meta;

/** Edits as patches on one host stack with a stand-in structure; a drag as a preview; the printed spec beside. */
export const Editing: StoryObj = {
  render: () => {
    const el = document.createElement(EditingRoot.tag) as EditingRoot;
    el.model = P53Views.find((v) => v.id === 'tracks')!;
    return el;
  },
};
