/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { DefaultSeqViewerSpec } from '@molstar/seqx-viewer/default-spec';
import { createSeqViewer, type SeqViewer } from '@molstar/seqx-viewer/viewer';
import { defineElements, onCleanup } from '@molstar/ui/dom';
import { UiElement } from '@molstar/ui/elements/ui-element';
import { HostSection, hostDrawnSpec } from './compose.js';
import { wireHostNavigator } from './host-navigator-dom.js';

// The host-drawn navigator story (the viewer's _spec "Models", NavigatorModel, the contract): a navigator the host
// draws itself in plain DOM, no canvas (`host-navigator-dom.ts`), above the viewer it moves; the viewer's
// `navigator(def)` model holds every rule. The viewer ends with the story.

export interface HostNavigatorStory {
  readonly hint: string;
}

/** The navigator's height (CSS px): the host's choice. */
const Height = 36;

class HostNavigatorRoot extends UiElement<HostNavigatorStory> {
  static readonly tag = 'seqx-host-navigator-root';

  protected render(story: HostNavigatorStory): HTMLElement {
    const root = document.createElement('div');
    root.style.cssText =
      'display: grid; gap: var(--msui-space-4); padding: var(--msui-space-5); min-height: 100vh; align-content: start; background: var(--msui-bg); color: var(--msui-text); font: var(--msui-text-sm) / var(--msui-leading-sm) var(--msui-font-sans);';
    const hint = document.createElement('p');
    hint.style.cssText = 'margin: 0; color: var(--msui-text-muted); max-width: 72ch;';
    hint.textContent = story.hint;
    const nav = document.createElement('div');
    nav.style.cssText = `height: ${Height}px; overflow: hidden; border-radius: var(--msui-radius-sm, 4px); background: var(--msui-surface-sunken, var(--msui-surface)); outline-offset: 2px;`;
    const readout = document.createElement('div');
    readout.style.cssText = 'min-height: 1.5em; color: var(--msui-text-muted); font-variant-numeric: tabular-nums;';
    const canvasHost = document.createElement('div');
    canvasHost.style.cssText = 'min-width: 0;';
    root.append(hint, nav, readout, canvasHost);

    const viewer: SeqViewer = createSeqViewer(DefaultSeqViewerSpec(), canvasHost);
    const model = viewer.navigator({ id: 'host', of: HostSection, content: 'axis' });
    const unwire = wireHostNavigator(nav, viewer, model, Height, readout);
    onCleanup(() => {
      unwire();
      viewer.retire();
    });
    void viewer.load(hostDrawnSpec());
    return root;
  }
}

defineElements([HostNavigatorRoot]);

/** A Storybook `render` result: the host's navigator above the viewer it moves. */
export function hostNavigatorStory(story: HostNavigatorStory): HTMLElement {
  const el = document.createElement(HostNavigatorRoot.tag) as HostNavigatorRoot;
  el.model = story;
  return el;
}
