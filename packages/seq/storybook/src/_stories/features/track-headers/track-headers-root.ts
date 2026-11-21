/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { ExampleView } from '@molstar/seqx-examples';
import type { HeaderTarget } from '@molstar/seqx-viewer/header-entries';
import { mountToolbar } from '@molstar/seqx-viewer/mount-toolbar';
import { createSeqViewer, type SeqViewer } from '@molstar/seqx-viewer/viewer';
import { Button } from '@molstar/ui/components/button/button';
import { defineElements, onCleanup } from '@molstar/ui/dom';
import { UiElement } from '@molstar/ui/elements/ui-element';
import { type HeadersComposition, headerLabel, headersOf, headersSpec, undoKey, widgetPlaces } from './compose.js';

// The track headers' stories (track-header-parts.md): a Seq* view composed by `headersSpec`, a hint, buttons that open a
// header's menu (the only way to a ring track's), and a host's own widgets placed from `viewer.headers`. The element is
// a UI root, so the theme's tokens apply; the viewer ends when Storybook removes the story.

export interface TrackHeadersStory extends HeadersComposition {
  readonly view: ExampleView;
  readonly hint: string;
  /** Buttons that open the menu of a header (its target) at the canvas's top left corner. */
  readonly menus?: readonly { readonly label: string; readonly target: HeaderTarget }[];
  /** The host's own DOM widget on every track: a numbered badge, a click opens that track's menu. */
  readonly widgets?: boolean;
}

const Muted = 'margin: 0; color: var(--msui-text-muted); max-width: 72ch;';

class TrackHeadersRoot extends UiElement<TrackHeadersStory> {
  static readonly tag = 'seqx-track-headers-root';

  protected render(story: TrackHeadersStory): HTMLElement {
    const root = document.createElement('div');
    root.style.cssText =
      'display: grid; gap: var(--msui-space-4); padding: var(--msui-space-5); min-height: 100vh; align-content: start; background: var(--msui-bg); color: var(--msui-text); font: var(--msui-text-sm) / var(--msui-leading-sm) var(--msui-font-sans);';
    const hint = document.createElement('p');
    hint.style.cssText = Muted;
    hint.textContent = story.hint;
    const tools = document.createElement('div');
    tools.style.cssText = 'display: flex; gap: var(--msui-space-2); flex-wrap: wrap; align-items: baseline;';
    // the widgets sit in the canvas's own left margin, so the canvas gets room for them
    const stage = document.createElement('div');
    stage.style.cssText = `position: relative; min-width: 0; ${story.widgets ? 'padding-left: 28px;' : ''}`;
    const canvasHost = document.createElement('div');
    canvasHost.style.cssText = 'min-width: 0; position: relative;';
    const widgetLayer = document.createElement('div');
    widgetLayer.style.cssText = 'position: absolute; left: 28px; top: 0; width: 0; height: 0; pointer-events: none;';
    stage.append(canvasHost, ...(story.widgets ? [widgetLayer] : []));
    const log = document.createElement('pre');
    log.style.cssText =
      'margin: 0; padding: var(--msui-space-3); min-height: 3em; font: var(--msui-text-xs) / var(--msui-leading-xs) var(--msui-font-mono); background: var(--msui-surface, transparent); border: 1px solid var(--msui-line);';
    log.textContent = 'The host’s command log';
    const lines: string[] = [];
    root.append(hint, tools, stage, ...(story.host ? [log] : []));

    const spec = headersSpec({
      ...story,
      registry: story.view.registry?.(),
      say: (line) => {
        lines.push(line);
        log.textContent = lines.slice(-8).join('\n');
      },
    });
    const viewer = createSeqViewer({ ...spec, ui: { ...spec.ui, toolbar: mountToolbar } }, canvasHost);
    let off: (() => void) | undefined;
    // the host's Undo and Redo keys, on the page: a click on a header's part leaves the focus off the canvas
    const onKey = (e: KeyboardEvent) => {
      if (!(e.target instanceof HTMLInputElement) && undoKey(viewer, e)) e.preventDefault();
    };
    document.addEventListener('keydown', onKey);
    onCleanup(() => {
      document.removeEventListener('keydown', onKey);
      off?.();
      viewer.retire();
    });

    // a ring track has no line to point at: the menu opens at the canvas's top left, and this line names its track
    const menuOf = document.createElement('span');
    menuOf.style.cssText = 'color: var(--msui-text-muted);';
    for (const m of story.menus ?? []) {
      tools.append(
        Button({
          children: m.label,
          onClick: () => {
            menuOf.textContent = `The menu open is the track “${headerLabel(viewer, m.target)}”’s (${m.target.track}).`;
            void headersOf(viewer).openMenu(m.target, { x: 120, y: 60 });
          },
        }) as HTMLElement,
      );
    }
    if (story.menus?.length) tools.append(menuOf);
    if (story.widgets) off = placeWidgets(viewer, stage, widgetLayer);
    void viewer.load(story.view.build()).then(() => viewer.settle());
    return root;
  }
}

/** A host's own widgets over the headers' geometry: a numbered badge left of each shown track's label. */
function placeWidgets(viewer: SeqViewer, stage: HTMLElement, layer: HTMLElement): () => void {
  const show = () => {
    // the geometry is in the canvas's own frame: the layer sits at the canvas's corner, whatever is above it (the toolbar)
    const canvas = stage.querySelector('canvas');
    if (canvas) {
      const at = canvas.getBoundingClientRect();
      const from = stage.getBoundingClientRect();
      layer.style.left = `${at.left - from.left}px`;
      layer.style.top = `${at.top - from.top}px`;
    }
    const places = widgetPlaces(viewer.headers?.value ?? []);
    layer.replaceChildren(
      ...places.map((p, i) => {
        const b = document.createElement('button');
        b.type = 'button';
        b.textContent = String(i + 1);
        b.title = `Menu of ${p.target.track}`;
        b.style.cssText = `position: absolute; left: ${p.x}px; top: ${p.y}px; width: ${p.size}px; height: ${p.size}px; padding: 0; pointer-events: auto; font: var(--msui-text-xs) / 1 var(--msui-font-sans); color: var(--msui-text); background: var(--msui-surface, transparent); border: 1px solid var(--msui-line); border-radius: 4px; cursor: pointer;`;
        b.addEventListener('click', () => void headersOf(viewer).openMenu(p.target, { x: 120, y: p.y + p.size }));
        return b;
      }),
    );
  };
  const off = viewer.headers?.listen(show);
  show();
  return () => off?.();
}

defineElements([TrackHeadersRoot]);

/** A Storybook `render` result: the view with its headers composed as the story says. */
export function trackHeadersStory(story: TrackHeadersStory): HTMLElement {
  const el = document.createElement(TrackHeadersRoot.tag) as TrackHeadersRoot;
  el.model = story;
  return el;
}
