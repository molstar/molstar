/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type ExampleView, PlasmidViews } from '@molstar/seqx-examples';
import { Button } from '@molstar/ui/components/button/button';
import { defineElements, onCleanup } from '@molstar/ui/dom';
import { UiElement } from '@molstar/ui/elements/ui-element';
import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { composeEditing, historyRows, lastMove, printedSpec } from '../../features/editing/compose.js';
import { followSvg } from '../../follow-svg.js';
import { bookmarkRange, rangeOf, setStart, showEnzymes, toggleStrand } from './editing.js';

// The Plasmid's editing views: the ring's strand tracks hidden, its sites filtered to some enzymes, bla bookmarked from
// its start to its stop across the origin; the rotated linear view's start set. Every edit is one patch and one step of
// the host's stack (Undo and Redo walk them), the printed spec and the SVG follow what is shown.

const Row = 'display: flex; gap: var(--msui-space-2); align-items: center; flex-wrap: wrap;';
const Mono = 'margin: 0; font: var(--msui-text-xs) / var(--msui-leading-xs) var(--msui-font-mono);';
const BlaStart = 2311;
const BlaEnd = 486;

class PlasmidEditingRoot extends UiElement<ExampleView> {
  static readonly tag = 'seqx-plasmid-editing-root';

  protected render(view: ExampleView): HTMLElement {
    const root = document.createElement('div');
    root.style.cssText =
      'display: grid; align-content: start; gap: var(--msui-space-4); padding: var(--msui-space-5); min-height: 100vh; background: var(--msui-bg); color: var(--msui-text); font: var(--msui-text-sm) / var(--msui-leading-sm) var(--msui-font-sans);';
    const canvasHost = document.createElement('div');
    canvasHost.style.cssText = 'min-width: 0;';
    const svgHost = document.createElement('div');
    svgHost.style.cssText = 'min-width: 0; overflow-x: auto;';
    const status = document.createElement('p');
    status.style.cssText = 'margin: 0; color: var(--msui-text-muted);';
    const say = (p: Promise<unknown> | string) =>
      void Promise.resolve(p).then((m) => {
        status.textContent = String(m);
      });
    const s = composeEditing(view, { target: canvasHost });
    onCleanup(() => s.dispose());

    const tools = document.createElement('div');
    tools.style.cssText = Row;
    const button = (children: string, onClick: () => void) => Button({ size: 'sm', children, onClick }) as HTMLElement;
    if (view.id === 'rotated-linear') {
      for (const start of [2000, 2400, 100])
        tools.append(
          button(`Start at ${start}`, () => say(setStart(s, start, 1386).then((ok) => (ok ? 'Set' : 'Refused')))),
        );
    } else {
      tools.append(
        button('Hide/show Features +', () => say(toggleStrand(s, '+'))),
        button('Hide/show Features −', () => say(toggleStrand(s, '-'))),
        button('Only EcoRI, HindIII', () =>
          say(showEnzymes(s, ['EcoRI', 'HindIII']).then((ok) => (ok ? 'Set' : 'Refused'))),
        ),
        button('Every enzyme', () => say(showEnzymes(s).then((ok) => (ok ? 'Set' : 'Refused')))),
        button('Bookmark bla', () => say(bookmarkRange(s, BlaStart, BlaEnd, 'bla-range', 'bla (AmpR), start to stop'))),
      );
    }
    const undo = Button({ size: 'sm', children: 'Undo', onClick: () => void s.stack.undo() }) as HTMLButtonElement;
    const redo = Button({ size: 'sm', children: 'Redo', onClick: () => void s.stack.redo() }) as HTMLButtonElement;
    tools.append(undo, redo);

    const history = document.createElement('pre');
    history.style.cssText = Mono;
    const printed = document.createElement('pre');
    printed.style.cssText = `${Mono} max-height: 320px; overflow: auto; background: var(--msui-bg-subtle); padding: var(--msui-space-3);`;
    const panes = document.createElement('div');
    panes.style.cssText =
      'display: grid; gap: var(--msui-space-4); grid-template-columns: minmax(0, 1fr) minmax(0, 1fr);';
    panes.append(history, printed);
    root.append(canvasHost, tools, status, panes, svgHost);

    const svg = followSvg(s.viewer, svgHost, { prefix: 'plasmid-edit' });
    onCleanup(() => svg.dispose());
    const draw = () => {
      const rows = historyRows(s);
      history.textContent = [
        'Undo (newest first):',
        ...rows.done.map((l) => `  ${l}`),
        'Redo:',
        ...rows.redo.map((l) => `  ${l}`),
        view.id === 'rotated-linear' ? `range: ${JSON.stringify(rangeOf(s))}` : '',
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

defineElements([PlasmidEditingRoot]);

const meta: Meta = { title: 'Examples/Plasmid/Editing' };
export default meta;

const story = (id: string): StoryObj => ({
  render: () => {
    const el = document.createElement(PlasmidEditingRoot.tag) as PlasmidEditingRoot;
    el.model = PlasmidViews.find((v) => v.id === id)!;
    return el;
  },
});

/**
 * The ring edited: hide a strand's track, show the sites of two enzymes only (a filter in the sites layer's data),
 * bookmark bla from its start to its stop across the origin. Each is one patch, one step of Undo and Redo.
 */
export const EditingTheRing: StoryObj = story('ring');

/** The rotated linear view's start set: its range is a section param, a patch, undone. */
export const EditingTheStart: StoryObj = story('rotated-linear');
