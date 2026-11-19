/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { ExampleView } from '@molstar/seqx-examples';
import { DefaultSeqViewerSpec } from '@molstar/seqx-viewer/default-spec';
import { placeCard } from '@molstar/seqx-viewer/hover-card';
import { describeKeys } from '@molstar/seqx-viewer/keys';
import { mountToolbar } from '@molstar/seqx-viewer/mount-toolbar';
import { createSeqViewer } from '@molstar/seqx-viewer/viewer';
import { defineElements, onCleanup } from '@molstar/ui/dom';
import { UiElement } from '@molstar/ui/elements/ui-element';
import { hostCardLines, withCardPort } from './compose.js';

// The hover card's and the keys panel's stories: a view with the default spec (its card `@molstar/ui`'s), with a host's
// own card adapter in its place (`host`: ink on paper, monospace, the same data through the hover card port), or with
// the keys panel open and the host's own list of the same keys beside it (`keys`: `describeKeys`, the panel's data).

export interface HoverCardModel {
  readonly view: ExampleView;
  readonly hint: string;
  /** The host's own card in place of `@molstar/ui`'s. */
  readonly host?: boolean;
  /** The keys panel open at the start, and the keys' data printed beside the view. */
  readonly keys?: boolean;
}

const Mono = 'font: var(--msui-text-xs) / var(--msui-leading-xs) var(--msui-font-mono);';

class HoverCardRoot extends UiElement<HoverCardModel> {
  static readonly tag = 'seqx-hover-card-root';

  protected render(model: HoverCardModel): HTMLElement {
    const root = document.createElement('div');
    root.style.cssText =
      'display: grid; align-content: start; gap: var(--msui-space-4); padding: var(--msui-space-5); min-height: 100vh; background: var(--msui-bg); color: var(--msui-text); font: var(--msui-text-sm) / var(--msui-leading-sm) var(--msui-font-sans);';
    const hint = document.createElement('p');
    hint.style.cssText = 'margin: 0; color: var(--msui-text-muted); max-width: 72ch;';
    hint.textContent = model.hint;
    const canvasHost = document.createElement('div');
    canvasHost.style.cssText = 'min-width: 0;';
    root.append(hint, canvasHost);

    const base = DefaultSeqViewerSpec();
    const spec = { ...base, ui: { ...base.ui, toolbar: mountToolbar } };
    const viewer = createSeqViewer(
      model.host
        ? withCardPort(spec, (container) => {
            // the host's card: its own element and styles, placed by the viewer's `placeCard` at its own size
            const card = document.createElement('pre');
            card.setAttribute('aria-hidden', 'true');
            card.style.cssText = `position: absolute; left: 0; top: 0; margin: 0; display: none; pointer-events: none; padding: 6px 8px; border-radius: 2px; background: var(--msui-text); color: var(--msui-bg); ${Mono}`;
            container?.append(card);
            return {
              show(c, place) {
                card.textContent = hostCardLines(c).join('\n');
                card.style.display = 'block';
                const at = placeCard(place.anchor, place.pointer, place.view, {
                  width: card.offsetWidth,
                  height: card.offsetHeight,
                });
                card.style.transform = `translate(${Math.round(at.x)}px, ${Math.round(at.y)}px)`;
              },
              hide() {
                card.style.display = 'none';
              },
              dispose() {
                card.remove();
              },
            };
          })
        : spec,
      canvasHost,
    );
    onCleanup(() => viewer.retire());

    if (model.keys) {
      // the host's own list from the same data the panel draws
      const list = document.createElement('pre');
      list.style.cssText = `margin: 0; padding: var(--msui-space-3); overflow: auto; max-height: 360px; background: var(--msui-surface); border-radius: var(--msui-radius-md); ${Mono}`;
      list.textContent = describeKeys(viewer)
        .map((g) =>
          [g.title, ...g.entries.map((e) => `  ${e.name}: ${[...e.keys, ...e.gestures].join(', ')}`)].join('\n'),
        )
        .join('\n\n');
      root.append(list);
    }
    void viewer
      .load(model.view.build())
      .then(() => viewer.settle())
      .then(() => {
        if (model.keys) void viewer.runCommand('help.keys', {});
      });
    return root;
  }
}

defineElements([HoverCardRoot]);

export function hoverCardStory(model: HoverCardModel): HTMLElement {
  const el = document.createElement(HoverCardRoot.tag) as HoverCardRoot;
  el.model = model;
  return el;
}
