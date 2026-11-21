/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { createSeqViewer, type SeqViewer } from '@molstar/seqx-viewer/viewer';
import { Button } from '@molstar/ui/components/button/button';
import { defineElements, onCleanup } from '@molstar/ui/dom';
import { UiElement } from '@molstar/ui/elements/ui-element';
import { DemoExtension, extensibilityComposition, extensibilitySpec, issueLines } from './compose.js';

// The Extensibility story (_spec "Composition", Extensions): a viewer composed without the example extension, a button
// that registers it at run time (`extend`) and takes it away again (its undo), and what the view says under the
// canvas. The element is a UI root, so the theme's tokens apply; the viewer ends when Storybook removes the story.

export interface ExtensibilityStory {
  readonly hint: string;
}

class ExtensibilityRoot extends UiElement<ExtensibilityStory> {
  static readonly tag = 'seqx-extensibility-root';

  protected render(story: ExtensibilityStory): HTMLElement {
    const root = document.createElement('div');
    root.style.cssText =
      'display: grid; gap: var(--msui-space-4); padding: var(--msui-space-5); min-height: 100vh; align-content: start; background: var(--msui-bg); color: var(--msui-text); font: var(--msui-text-sm) / var(--msui-leading-sm) var(--msui-font-sans);';
    const hint = document.createElement('p');
    hint.style.cssText = 'margin: 0; color: var(--msui-text-muted); max-width: 72ch;';
    hint.textContent = story.hint;
    const tools = document.createElement('div');
    tools.style.cssText = 'display: flex; gap: var(--msui-space-3); align-items: baseline;';
    const state = document.createElement('span');
    state.style.cssText = 'color: var(--msui-text-muted);';
    const canvasHost = document.createElement('div');
    canvasHost.style.cssText = 'min-width: 0;';
    const said = document.createElement('pre');
    said.style.cssText =
      'margin: 0; padding: var(--msui-space-3); overflow: auto; white-space: pre-wrap; font: var(--msui-text-xs) / var(--msui-leading-xs) var(--msui-font-mono); border: 1px solid var(--msui-line);';
    root.append(hint, tools, canvasHost, said);

    const viewer: SeqViewer = createSeqViewer(extensibilityComposition(true), canvasHost);
    let undo: (() => void) | undefined;
    let busy = false;
    const show = () => {
      state.textContent = undo ? 'The `demo` extension is registered.' : 'The `demo` extension is not registered.';
      said.textContent = issueLines(viewer).join('\n');
    };
    const offs = [viewer.view.resolved.listen(show), viewer.layoutIssues.listen(show)];
    onCleanup(() => {
      for (const off of offs) off();
      viewer.retire();
    });

    const toggle = Button({
      children: 'Register the extension',
      onClick: async () => {
        if (busy) return;
        busy = true;
        try {
          if (undo) {
            undo();
            undo = undefined;
          } else undo = await viewer.extend(DemoExtension);
        } finally {
          busy = false;
        }
        (toggle as HTMLElement).textContent = undo ? 'Take the extension away' : 'Register the extension';
        show();
      },
    }) as HTMLElement;
    tools.append(toggle, state);
    show();
    void viewer.load(extensibilitySpec()).then(show);
    return root;
  }
}

defineElements([ExtensibilityRoot]);

/** A Storybook `render` result: the viewer, its toggle and what it says. */
export function extensibilityStory(story: ExtensibilityStory): HTMLElement {
  const el = document.createElement(ExtensibilityRoot.tag) as ExtensibilityRoot;
  el.model = story;
  return el;
}
