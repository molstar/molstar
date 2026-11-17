/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { printSexpr } from '@molstar/core/repr/codecs/sexpr';
import { SeqLanguage } from '@molstar/seqx-core/lang';
import type { ExampleView } from '@molstar/seqx-examples';
import { DefaultSeqViewerSpec } from '@molstar/seqx-viewer/default-spec';
import { createSeqViewer } from '@molstar/seqx-viewer/viewer';
import { defineElements, DomHost, onCleanup } from '@molstar/ui/dom';
import { UiElement } from '@molstar/ui/elements/ui-element';

// A story of an example view (examples-and-stories.md): the view's `build()` shown by a viewer composed from the default
// spec, the canvas and the SVG of the same view side by side, and the spec as it prints, beside them. The element is
// a UI root, so the theme's tokens apply; the viewer is retired when Storybook removes the story.

const Box = 'display: grid; gap: var(--msui-space-2); min-width: 0; align-content: start;';
const Title =
  'margin: 0; font-size: var(--msui-text-sm); font-weight: var(--msui-weight-strong); color: var(--msui-text-muted);';

function section(title: string, body: HTMLElement): HTMLElement {
  const el = document.createElement('div');
  el.style.cssText = Box;
  const h = document.createElement('h4');
  h.style.cssText = Title;
  h.textContent = title;
  el.append(h, body);
  return el;
}

/** The viewer and the SVG of one view, and its printed spec. */
class ExampleRoot extends UiElement<ExampleView> {
  static readonly tag = 'seqx-example-root';

  protected render(view: ExampleView): HTMLElement {
    const root = document.createElement('div');
    root.style.cssText =
      'display: grid; gap: var(--msui-space-5); padding: var(--msui-space-5); min-height: 100vh; background: var(--msui-bg); color: var(--msui-text); font: var(--msui-text-sm) / var(--msui-leading-sm) var(--msui-font-sans);';

    const summary = document.createElement('p');
    summary.style.cssText = 'margin: 0; color: var(--msui-text-muted); max-width: 72ch;';
    summary.textContent = view.summary;

    const canvasHost = document.createElement('div');
    canvasHost.style.cssText = 'min-width: 0;';
    const svgHost = document.createElement('div');
    svgHost.style.cssText = 'min-width: 0; overflow-x: auto;';
    const printed = document.createElement('pre');
    printed.style.cssText =
      'margin: 0; padding: var(--msui-space-3); overflow: auto; max-height: 360px; font: var(--msui-text-xs) / var(--msui-leading-xs) var(--msui-font-mono); background: var(--msui-surface, transparent); border: 1px solid var(--msui-line);';
    const status = document.createElement('p');
    status.style.cssText = 'margin: 0; color: var(--msui-text-muted);';

    const pair = document.createElement('div');
    pair.style.cssText =
      'display: grid; grid-template-columns: repeat(auto-fit, minmax(420px, 1fr)); gap: var(--msui-space-5);';
    pair.append(section('Canvas', canvasHost), section('SVG', svgHost));
    root.append(summary, pair, status, section('Spec', printed));

    // the same spec the view builds, composed from the default spec: the fixtures are read relative to the page
    const viewer = createSeqViewer(DefaultSeqViewerSpec(), canvasHost);
    onCleanup(() => viewer.retire());
    const spec = view.build();
    printed.textContent = printSexpr(spec, { language: SeqLanguage });

    // the SVG follows the canvas's window, drawn once the view rests
    let timer: unknown;
    let alive = true;
    const drawSvg = async () => {
      const r = await viewer.renderSvg();
      if (!alive) return;
      // our own output: every string from data or the spec is escaped by the backend
      svgHost.innerHTML = r.svg;
      status.textContent = `${r.settled ? 'Settled' : 'Not settled'}; ${r.issues.length} issue${r.issues.length === 1 ? '' : 's'}${
        r.issues.length ? `: ${r.issues.map((i) => i.message).join('; ')}` : ''
      }`;
    };
    const later = () => {
      DomHost.clearTimeout(timer);
      timer = DomHost.setTimeout(() => void drawSvg(), 200);
    };
    onCleanup(() => {
      alive = false;
      DomHost.clearTimeout(timer);
    });
    viewer.view.resolved.listen(later);
    viewer.layoutIssues.listen(later);
    viewer.store.store.settled.listen(later);
    void viewer
      .load(spec)
      .then(() => viewer.settle())
      .then(drawSvg);
    return root;
  }
}

defineElements([ExampleRoot]);

/** A Storybook `render` result for an example view. */
export function exampleStory(view: ExampleView): HTMLElement {
  const el = document.createElement(ExampleRoot.tag) as ExampleRoot;
  el.model = view;
  return el;
}
