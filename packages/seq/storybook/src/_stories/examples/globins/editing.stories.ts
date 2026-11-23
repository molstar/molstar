/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { GlobinsViews } from '@molstar/seqx-examples';
import { Button } from '@molstar/ui/components/button/button';
import { defineElements, onCleanup } from '@molstar/ui/dom';
import { UiElement } from '@molstar/ui/elements/ui-element';
import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { composeEditing, historyRows, lastMove, printedSpec } from '../../features/editing/compose.js';
import { followSvg } from '../../follow-svg.js';
import { keepFolds, pinRows, setScheme, toggleCopies, toggleLogo } from './editing.js';

// The Globins' editing views: the family's residue scheme and its logo, the folded
// clades kept and a row pinned, hemoglobin's copies split and merged; each edit one patch and one step of a host's
// stack, Undo and Redo walking them. The printed spec and the SVG export follow what is shown.

const Row = 'display: flex; gap: var(--msui-space-2); align-items: center; flex-wrap: wrap;';
const Mono = 'margin: 0; font: var(--msui-text-xs) / var(--msui-leading-xs) var(--msui-font-mono);';

type Edit = 'family' | 'folded' | 'hemoglobin';

class GlobinsEditingRoot extends UiElement<Edit> {
  static readonly tag = 'seqx-globins-editing-root';

  protected render(id: Edit): HTMLElement {
    const view = GlobinsViews.find((v) => v.id === id)!;
    const root = document.createElement('div');
    root.style.cssText =
      'display: grid; align-content: start; gap: var(--msui-space-4); padding: var(--msui-space-5); min-height: 100vh; background: var(--msui-bg); color: var(--msui-text); font: var(--msui-text-sm) / var(--msui-leading-sm) var(--msui-font-sans);';
    const canvasHost = document.createElement('div');
    canvasHost.style.cssText = 'min-width: 0;';
    const svgHost = document.createElement('div');
    svgHost.style.cssText = 'min-width: 0; overflow-x: auto;';
    const status = document.createElement('p');
    status.style.cssText = 'margin: 0; color: var(--msui-text-muted);';
    const say = (p: Promise<string>) =>
      void p.then((m) => {
        status.textContent = m;
      });
    const s = composeEditing(view, { target: canvasHost });
    onCleanup(() => s.dispose());

    const tools = document.createElement('div');
    tools.style.cssText = Row;
    const button = (children: string, onClick: () => void) => Button({ size: 'sm', children, onClick }) as HTMLElement;
    if (id === 'family')
      tools.append(
        button('Percent identity', () => say(setScheme(s, 'percent-identity'))),
        button('ESPript', () => say(setScheme(s, 'espript'))),
        button('Hide/show the logo', () => say(toggleLogo(s))),
      );
    if (id === 'folded')
      tools.append(
        button('Keep the folds shown', () => say(keepFolds(s))),
        button('Pin the first row', () => say(pinRows(s, [0]))),
      );
    if (id === 'hemoglobin') tools.append(button('Split/merge the copies', () => say(toggleCopies(s))));
    const undo = button('Undo', () => void s.stack.undo());
    const redo = button('Redo', () => void s.stack.redo());
    tools.append(undo, redo);

    const history = document.createElement('pre');
    history.style.cssText = Mono;
    const printed = document.createElement('pre');
    printed.style.cssText = `${Mono} max-height: 320px; overflow: auto; background: var(--msui-bg-subtle); padding: var(--msui-space-3);`;
    const panes = document.createElement('div');
    panes.style.cssText =
      'display: grid; gap: var(--msui-space-4); grid-template-columns: repeat(auto-fit, minmax(320px, 1fr));';
    panes.append(history, printed);
    root.append(canvasHost, tools, status, panes, svgHost);

    const svg = followSvg(s.viewer, svgHost, { prefix: `globins-${id}` });
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
      (undo as HTMLButtonElement).disabled = rows.done.length === 0;
      (redo as HTMLButtonElement).disabled = rows.redo.length === 0;
    };
    for (const state of [s.stack.steps, s.stack.redoable, s.viewer.spec.shown]) onCleanup(state.listen(draw));
    onCleanup(
      s.stack.lastMove.listen(() => {
        status.textContent = lastMove(s);
      }),
    );
    void s.viewer
      .load(view.build())
      .then(() => s.viewer.settle())
      .then(() => {
        draw();
        return svg.draw();
      })
      .catch((e: unknown) => {
        status.textContent = `The view did not load: ${e instanceof Error ? e.message : String(e)}`;
      });
    return root;
  }
}

defineElements([GlobinsEditingRoot]);

const meta: Meta = { title: 'Examples/Globins/Editing' };
export default meta;

const story = (id: Edit): StoryObj => ({
  render: () => {
    const el = document.createElement(GlobinsEditingRoot.tag) as GlobinsEditingRoot;
    el.model = id;
    return el;
  },
});

/** The family edited: its residue scheme (a theme patch) and its logo hidden or shown; Undo and Redo. */
export const Family = story('family');

/** The folded clades: click a node to fold or unfold (session state), keep the folds shown (one patch), pin a row. */
export const Folded = story('folded');

/** Hemoglobin's copies split into a block a chain and merged again into ×2: one patch each way. */
export const Hemoglobin = story('hemoglobin');
