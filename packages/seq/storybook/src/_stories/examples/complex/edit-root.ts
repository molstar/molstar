/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { ExampleView } from '@molstar/seqx-examples';
import { Button } from '@molstar/ui/components/button/button';
import { defineElements, onCleanup } from '@molstar/ui/dom';
import { UiElement } from '@molstar/ui/elements/ui-element';
import {
  type EditingSession,
  composeEditing,
  historyRows,
  lastMove,
  printedSpec,
} from '../../features/editing/compose.js';
import { followSvg } from '../../follow-svg.js';
import type { SeqViewer } from '@molstar/seqx-viewer/viewer';
import type { SeqViewerSpec } from '@molstar/seqx-viewer/viewer-spec';
import { holdReadout, readout } from './compose.js';

// A Complex view with its edits: the viewer, a button an edit (one
// patch, one step of the host's stack), Undo and Redo, the history and the spec shown printed beside, the SVG export
// following; with `readout`, the host's panel of a grid over the complex (the selected rect and its mean PAE: the
// hovered cell is the viewer's hover card's), or with `lines` another dataset's panel (the Antibody's selection and
// bookmarks). Every string from the spec or data is put in as text, never markup.

/** An edit a button runs: its label and what it does (the status line says what it returns). */
export interface ComplexEdit {
  readonly label: string;
  run(s: EditingSession): Promise<string> | string;
}

export interface ComplexEditModel {
  readonly view: ExampleView;
  readonly edits: readonly ComplexEdit[];
  /** The grid's readout panel (the PAE over the complex). */
  readonly readout?: boolean;
  /** Another dataset's panel instead: its lines, from what the pointer is over and what is selected. */
  readonly lines?: (viewer: SeqViewer) => readonly string[];
  /** Behaviors besides the default spec's (`links.hover`). */
  readonly behaviors?: SeqViewerSpec['behaviors'];
}

const Row = 'display: flex; gap: var(--msui-space-2); align-items: center; flex-wrap: wrap;';
const Mono = 'margin: 0; font: var(--msui-text-xs) / var(--msui-leading-xs) var(--msui-font-mono);';

class ComplexEditRoot extends UiElement<ComplexEditModel> {
  static readonly tag = 'seqx-complex-edit-root';

  protected render(model: ComplexEditModel): HTMLElement {
    const root = document.createElement('div');
    root.style.cssText =
      'display: grid; align-content: start; gap: var(--msui-space-4); padding: var(--msui-space-5); min-height: 100vh; background: var(--msui-bg); color: var(--msui-text); font: var(--msui-text-sm) / var(--msui-leading-sm) var(--msui-font-sans);';
    const summary = document.createElement('p');
    summary.style.cssText = 'margin: 0; color: var(--msui-text-muted); max-width: 72ch;';
    summary.textContent = model.view.summary;
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

    const s = composeEditing(model.view, {
      target: canvasHost,
      ...(model.behaviors ? { behaviors: model.behaviors } : {}),
    });
    onCleanup(() => s.dispose());

    const tools = document.createElement('div');
    tools.style.cssText = Row;
    for (const e of model.edits)
      tools.append(Button({ size: 'sm', children: e.label, onClick: () => say(e.run(s)) }) as HTMLElement);
    const undo = Button({ size: 'sm', children: 'Undo', onClick: () => void s.stack.undo() }) as HTMLButtonElement;
    const redo = Button({ size: 'sm', children: 'Redo', onClick: () => void s.stack.redo() }) as HTMLButtonElement;
    if (model.edits.length) tools.append(undo, redo);

    // the host's panel: what is selected, in words (what the pointer is over is the viewer's hover card)
    const panel = document.createElement('div');
    panel.style.cssText = `display: grid; gap: var(--msui-space-1); padding: var(--msui-space-3); border: 1px solid var(--msui-line); border-radius: var(--msui-radius-md); ${Mono}`;
    const rectLine = document.createElement('div');
    panel.append(rectLine);
    const showReadout = () => {
      if (model.lines) {
        panel.replaceChildren(
          ...model.lines(s.viewer).map((l) => {
            const div = document.createElement('div');
            div.textContent = l;
            return div;
          }),
        );
        return;
      }
      const r = readout(s.viewer);
      rectLine.textContent = r.rect ? r.rect : 'Drag a rectangle on the plot: its residues and its mean PAE';
    };
    const panelShown = model.readout || model.lines !== undefined;
    if (panelShown) {
      // the letters the panel names are held while it is shown; it is written again when they arrive
      if (model.readout) onCleanup(holdReadout(s.viewer));
      onCleanup(s.viewer.selection.value.listen(showReadout));
      onCleanup(s.viewer.store.store.version.listen(showReadout));
      // the bookmarks a panel lists
      if (model.lines) onCleanup(s.viewer.addresses.refs.listen(showReadout));
    }

    const history = document.createElement('pre');
    history.style.cssText = Mono;
    const printed = document.createElement('pre');
    printed.style.cssText = `${Mono} max-height: 320px; overflow: auto; background: var(--msui-bg-subtle); padding: var(--msui-space-3);`;
    const panes = document.createElement('div');
    panes.style.cssText =
      'display: grid; gap: var(--msui-space-4); grid-template-columns: minmax(0, 1fr) minmax(0, 1fr);';
    panes.append(history, printed);
    root.append(summary, canvasHost, ...(panelShown ? [panel] : []), tools, status, panes, svgHost);

    const svg = followSvg(s.viewer, svgHost, { prefix: 'cx' });
    onCleanup(() => svg.dispose());
    const draw = () => {
      const rows = historyRows(s);
      history.textContent = [
        'Undo (newest first):',
        ...rows.done.map((l) => `  ${l}`),
        'Redo:',
        ...rows.redo.map((l) => `  ${l}`),
      ].join('\n');
      printed.textContent = printedSpec(s.viewer.spec.shown.value);
      undo.disabled = rows.done.length === 0;
      redo.disabled = rows.redo.length === 0;
    };
    for (const state of [s.stack.steps, s.stack.redoable, s.viewer.spec.shown]) onCleanup(state.listen(draw));
    onCleanup(
      s.stack.lastMove.listen(() => {
        status.textContent = lastMove(s);
      }),
    );
    void s.viewer
      .load(model.view.build())
      .then(() => s.viewer.settle())
      .then(() => {
        draw();
        if (panelShown) showReadout();
        return svg.draw();
      })
      .catch((e: unknown) => {
        status.textContent = `The view did not load: ${e instanceof Error ? e.message : String(e)}`;
      });
    return root;
  }
}

defineElements([ComplexEditRoot]);

/** A Storybook `render` result: a Complex view with its edits (and the readout panel when asked). */
export function complexEditStory(model: ComplexEditModel): HTMLElement {
  const el = document.createElement(ComplexEditRoot.tag) as ComplexEditRoot;
  el.model = model;
  return el;
}
