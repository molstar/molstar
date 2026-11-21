/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { printSexpr } from '@molstar/core/repr/codecs/sexpr';
import { SeqLanguage } from '@molstar/seqx-core/lang';
import type { ExampleView } from '@molstar/seqx-examples';
import { DefaultSeqViewerSpec } from '@molstar/seqx-viewer/default-spec';
import { toolbarWithHiddenTracks } from '../support/toolbar.js';
import { createSeqViewer, type SeqViewer } from '@molstar/seqx-viewer/viewer';
import type { SeqViewerSpec } from '@molstar/seqx-viewer/viewer-spec';
import { defineElements, onCleanup } from '@molstar/ui/dom';
import { UiElement } from '@molstar/ui/elements/ui-element';
import { followSvg } from '../follow-svg.js';
import { statusIssues } from './status-issues.js';
import { filesOf, fileHead, printSources } from '../features/files/file-head.js';

// A story of an example view (examples-and-stories.md): the view's `build()` shown by a viewer composed from the default
// spec, the canvas and the SVG of the same view side by side, and the spec as it prints, beside them (the SVG with the
// view's navigators where the story asks). The element is a UI root, so the theme's tokens apply; the viewer is retired
// when Storybook removes the story.

/** A story's view, and whether its SVG draws the view's navigators (`renderSvg`'s `navigators`). */
type ExampleModel = ExampleView & {
  readonly navigators?: 'omit' | 'include';
  readonly files?: boolean;
  /** The canvas above the SVG, each at the full width (a figure whose letters need it), not side by side. */
  readonly stacked?: boolean;
  /** What the story does once the view is placed (a move to the region it opens on). */
  readonly opening?: (viewer: SeqViewer) => unknown;
  /** The viewer's config items (`[item, value]` pairs) over the default spec's. */
  readonly config?: SeqViewerSpec['config'];
  /** Behaviors besides the default spec's (a host's choice the story shows on: `links.hover`). */
  readonly behaviors?: SeqViewerSpec['behaviors'];
};

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
class ExampleRoot extends UiElement<ExampleModel> {
  static readonly tag = 'seqx-example-root';

  protected render(view: ExampleModel): HTMLElement {
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
    pair.style.cssText = `display: grid; grid-template-columns: ${view.stacked ? '1fr' : 'repeat(auto-fit, minmax(420px, 1fr))'}; gap: var(--msui-space-5);`;
    pair.append(section('Canvas', canvasHost), section('SVG', svgHost));
    root.append(summary, pair, status);
    if (view.files) root.append(...fileSections(view.build()));
    root.append(section('Spec', printed));

    // the same spec the view builds, composed from the default spec: the fixtures are read relative to the page
    const base = DefaultSeqViewerSpec();
    const extra = view.registry?.() ?? [];
    const viewer = createSeqViewer(
      {
        ...base,
        ui: { ...base.ui, toolbar: toolbarWithHiddenTracks },
        registry: [...(base.registry ?? []), ...extra],
        ...(view.config ? { config: [...(base.config ?? []), ...view.config] } : {}),
        ...(view.behaviors ? { behaviors: [...base.behaviors, ...view.behaviors] } : {}),
      },
      canvasHost,
    );
    onCleanup(() => viewer.retire());
    const spec = view.build();
    printed.textContent = printSexpr(spec, { language: SeqLanguage });

    // the SVG follows the canvas (follow-svg.ts): drawn when the view comes to rest, on an edit or a fold, and when the
    // theme changes (the viewer's tokens are read from this root's CSS), by its own id prefix so SVGs in one page never
    // share a clip or a mask
    const svg = followSvg(viewer, svgHost, {
      prefix: 'ex',
      also: [viewer.layoutIssues, viewer.modelIssues],
      ...(view.navigators ? { svg: { navigators: view.navigators } } : {}),
      onDraw: (r) => {
        // the export's issues and the viewer's model's as it says them (one of a source still opening held back)
        const issues = statusIssues(r.issues, viewer.modelIssues.value);
        status.textContent = `${r.settled ? 'Settled' : 'Not settled'}; ${issues.length} issue${issues.length === 1 ? '' : 's'}${
          issues.length ? `: ${issues.map((i) => i.message).join('; ')}` : ''
        }`;
      },
    });
    onCleanup(() => svg.dispose());
    void viewer
      .load(spec)
      .then(() => viewer.settle())
      .then(() => view.opening?.(viewer))
      .then(() => svg.draw());
    return root;
  }
}

/** What a Files story prints beside the view: the source nodes (the transform) and the head of each file they name. */
function fileSections(spec: Parameters<typeof printSources>[0]): HTMLElement[] {
  const pre = (text: string) => {
    const el = document.createElement('pre');
    el.style.cssText =
      'margin: 0; padding: var(--msui-space-3); overflow: auto; max-height: 280px; font: var(--msui-text-xs) / var(--msui-leading-xs) var(--msui-font-mono); background: var(--msui-surface, transparent); border: 1px solid var(--msui-line);';
    el.textContent = text;
    return el;
  };
  const out = [section('Source node', pre(printSources(spec)))];
  for (const file of filesOf(spec)) {
    const head = pre('Reading…');
    out.push(section(`File: ${file.uri}`, head));
    // read the way the viewer reads it: relative to the page, which serves the examples' fixtures
    void fetch(file.uri)
      .then((r) => (r.ok ? r.arrayBuffer() : Promise.reject(new Error(`${r.status}`))))
      .then((b) => {
        head.textContent = fileHead(new Uint8Array(b), file.format);
      })
      .catch((e) => {
        head.textContent = `The file could not be read: ${e instanceof Error ? e.message : e}`;
      });
  }
  return out;
}

defineElements([ExampleRoot]);

/**
 * A Storybook `render` result for an example view (`navigators`: its SVG draws the view's navigators; `files`: the view's
 * source nodes and the head of each file they name are printed beside it; `config`: the viewer's items).
 */
export function exampleStory(
  view: ExampleView,
  options: {
    readonly navigators?: 'omit' | 'include';
    readonly files?: boolean;
    readonly stacked?: boolean;
    readonly opening?: (viewer: SeqViewer) => unknown;
    readonly config?: SeqViewerSpec['config'];
    readonly behaviors?: SeqViewerSpec['behaviors'];
  } = {},
): HTMLElement {
  const el = document.createElement(ExampleRoot.tag) as ExampleRoot;
  el.model =
    options.navigators || options.files || options.stacked || options.opening || options.config || options.behaviors
      ? { ...view, ...options }
      : view;
  return el;
}
