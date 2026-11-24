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
import { defineElements, onCleanup } from '@molstar/ui/dom';
import { UiElement } from '@molstar/ui/elements/ui-element';
import { followSvg } from '../follow-svg.js';
import { pairGrid } from '../examples/example-layout.js';
import { withViewRegistry } from '../examples/view-registry.js';

// The Docs pages' live examples: one example view in a viewer composed from the default spec, as a host composes it,
// and, where a page asks, the code it was written with, the spec as it prints (S-expression, JSON) and its SVG export
// beside the canvas. A UI root, so the theme's tokens apply; the viewer is retired when Storybook removes the story.

/** What a page shows of a view besides the canvas. */
export interface LiveOptions {
  /** The code the view was written with (a module's source), shown above the view. */
  readonly code?: string;
  /** The spec as it prints, under the view. */
  readonly printed?: readonly ('sexpr' | 'json')[];
  /** The SVG export beside the canvas, drawn again whenever the canvas comes to rest. */
  readonly svg?: boolean;
}

type LiveModel = ExampleView & LiveOptions;

const Muted =
  'margin: 0; font-size: var(--msui-text-sm); font-weight: var(--msui-weight-strong); color: var(--msui-text-muted);';
const Pre =
  'margin: 0; padding: var(--msui-space-3); overflow: auto; max-height: 320px; font: var(--msui-text-xs) / var(--msui-leading-xs) var(--msui-font-mono); background: var(--msui-surface, transparent); color: var(--msui-text); border: 1px solid var(--msui-line);';

function titled(title: string, body: HTMLElement): HTMLElement {
  const box = document.createElement('div');
  box.style.cssText = 'display: grid; grid-template-columns: minmax(0, 1fr); gap: var(--msui-space-2); min-width: 0;';
  const h = document.createElement('p');
  h.style.cssText = Muted;
  h.textContent = title;
  box.append(h, body);
  return box;
}

/** A block of text that scrolls: a tab stop with a name, so a keyboard reaches what overflows (axe's rule). */
function pre(label: string, text: string): HTMLElement {
  const el = document.createElement('pre');
  el.style.cssText = Pre;
  el.tabIndex = 0;
  el.setAttribute('aria-label', label);
  el.textContent = text;
  return el;
}

class DocsLive extends UiElement<LiveModel> {
  static readonly tag = 'seqx-docs-live';

  protected render(view: LiveModel): HTMLElement {
    const root = document.createElement('div');
    root.style.cssText =
      'display: grid; grid-template-columns: minmax(0, 1fr); gap: var(--msui-space-4); padding: var(--msui-space-4); min-width: 0; background: var(--msui-bg); color: var(--msui-text); font: var(--msui-text-sm) / var(--msui-leading-sm) var(--msui-font-sans);';

    const spec = view.build();
    if (view.code) root.append(titled('The code', pre('The code', view.code)));

    const canvasHost = document.createElement('div');
    canvasHost.style.cssText = 'min-width: 0;';
    const svgHost = document.createElement('div');
    svgHost.style.cssText = 'min-width: 0; overflow-x: auto;';
    if (view.svg) {
      const pair = document.createElement('div');
      pair.style.cssText = pairGrid(false);
      pair.append(titled('Canvas', canvasHost), titled('SVG export', svgHost));
      root.append(pair);
    } else root.append(canvasHost);

    const status = document.createElement('p');
    status.style.cssText = 'margin: 0; color: var(--msui-text-muted);';
    status.setAttribute('role', 'status');
    root.append(status);

    for (const p of view.printed ?? []) {
      const text = p === 'json' ? JSON.stringify(spec, null, 2) : printSexpr(spec, { language: SeqLanguage });
      const title = p === 'json' ? 'The spec as JSON' : 'The spec as an S-expression';
      root.append(titled(title, pre(title, text)));
    }

    // the default spec, as a host spreads it: the view's own entries (a format's runtime) added; its fixtures are read
    // relative to the page
    const base = DefaultSeqViewerSpec();
    const viewer = createSeqViewer(
      { ...base, registry: withViewRegistry(base.registry ?? [], view.registry?.() ?? []) },
      canvasHost,
    );
    onCleanup(() => viewer.retire());
    const svg = view.svg ? followSvg(viewer, svgHost, { prefix: 'docs' }) : undefined;
    if (svg) onCleanup(() => svg.dispose());

    const say = () => {
      const issues = viewer.modelIssues.value;
      status.textContent = `${issues.length} issue${issues.length === 1 ? '' : 's'}${issues.length ? `: ${issues.map((i) => i.message).join('; ')}` : ''}`;
    };
    onCleanup(viewer.modelIssues.listen(say));
    status.textContent = 'Opening…';
    void viewer
      .load(spec)
      .then(() => viewer.settle())
      .then(() => {
        say();
        return svg?.draw();
      });
    return root;
  }
}

defineElements([DocsLive]);

/** A Storybook `render` result: the view live, with what `options` asks beside it. */
export function liveStory(view: ExampleView, options: LiveOptions = {}): HTMLElement {
  const el = document.createElement(DocsLive.tag) as DocsLive;
  el.model = { ...view, ...options };
  return el;
}
