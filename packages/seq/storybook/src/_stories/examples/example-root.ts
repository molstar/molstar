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
import { Box, PageGrid, pairGrid } from './example-layout.js';
import { statusIssues } from './status-issues.js';
import { withViewRegistry } from './view-registry.js';
import { filesOf, fileText, printSources } from '../features/files/file-head.js';

// A story of an example view: the view's `build()` shown by a viewer composed from the default
// spec, the canvas and the SVG of the same view side by side, and the spec as it prints, beside them (the SVG with the
// view's navigators where the story asks). The element is a UI root, so the theme's tokens apply; the viewer is retired
// when Storybook removes the story.

/** A story's view, and whether its SVG draws the view's navigators (`renderSvg`'s `navigators`). */
type ExampleModel = ExampleView & {
  readonly navigators?: 'omit' | 'include';
  readonly files?: boolean;
  /** The status line says what the view's sources said when they opened (a misread file's lines left out). */
  readonly sourceIssues?: boolean;
  /** The canvas above the SVG, each at the full width (a figure whose letters need it), not side by side. */
  readonly stacked?: boolean;
  /** What the story does once the view is placed (a move to the region it opens on). */
  readonly opening?: (viewer: SeqViewer) => unknown;
  /**
   * A find box above the view (a name typed, a Find button): `run` finds the text and says what it found and where the
   * view went, printed under the box; `initial` is found once the view is placed.
   */
  readonly find?: {
    readonly title: string;
    readonly initial: string;
    readonly run: (viewer: SeqViewer, text: string) => Promise<string>;
  };
  /** The viewer's config items (`[item, value]` pairs) over the default spec's. */
  readonly config?: SeqViewerSpec['config'];
  /** Behaviors besides the default spec's (a host's choice the story shows on: `links.hover`). */
  readonly behaviors?: SeqViewerSpec['behaviors'];
  /** The extensions the viewer is composed with (an extension's showcase: the genome ring's `rondoViewer`). */
  readonly extensions?: SeqViewerSpec['extensions'];
};

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
    root.style.cssText = `${PageGrid} background: var(--msui-bg); color: var(--msui-text); font: var(--msui-text-sm) / var(--msui-leading-sm) var(--msui-font-sans);`;

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
    pair.style.cssText = pairGrid(!!view.stacked);
    pair.append(section('Canvas', canvasHost), section('SVG', svgHost));
    const found = document.createElement('pre');
    const input = document.createElement('input');
    const button = document.createElement('button');
    if (view.find) {
      const box = document.createElement('div');
      box.style.cssText = 'display: flex; gap: var(--msui-space-2); align-items: center;';
      input.type = 'search';
      input.value = view.find.initial;
      input.setAttribute('aria-label', view.find.title);
      input.style.cssText =
        'font: inherit; padding: var(--msui-space-1) var(--msui-space-2); min-width: 16ch; color: var(--msui-text); background: var(--msui-bg); border: 1px solid var(--msui-line); border-radius: 4px;';
      button.textContent = 'Find';
      button.style.cssText =
        'font: inherit; padding: var(--msui-space-1) var(--msui-space-3); color: var(--msui-text); background: var(--msui-surface, transparent); border: 1px solid var(--msui-line); border-radius: 4px; cursor: pointer;';
      found.style.cssText =
        'margin: 0; padding: var(--msui-space-3); overflow: auto; font: var(--msui-text-xs) / var(--msui-leading-xs) var(--msui-font-mono); background: var(--msui-surface, transparent); border: 1px solid var(--msui-line);';
      found.textContent = 'Opening the file…';
      box.append(input, button);
      const wrap = document.createElement('div');
      wrap.style.cssText = Box;
      wrap.append(box, found);
      root.append(summary, section(view.find.title, wrap), pair, status);
    } else root.append(summary, pair, status);
    if (view.files) root.append(...fileSections(view.build()));
    root.append(section('Spec', printed));

    // the same spec the view builds, composed from the default spec: the fixtures are read relative to the page
    const base = DefaultSeqViewerSpec();
    const extra = view.registry?.() ?? [];
    const viewer = createSeqViewer(
      {
        ...base,
        ui: { ...base.ui, toolbar: toolbarWithHiddenTracks },
        registry: withViewRegistry(base.registry ?? [], extra),
        ...(view.config ? { config: [...(base.config ?? []), ...view.config] } : {}),
        ...(view.behaviors ? { behaviors: [...base.behaviors, ...view.behaviors] } : {}),
        ...(view.extensions ? { extensions: [...(base.extensions ?? []), ...view.extensions] } : {}),
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
      also: [viewer.layoutIssues, viewer.modelIssues, ...(view.sourceIssues ? [viewer.store.issues] : [])],
      ...(view.navigators ? { svg: { navigators: view.navigators } } : {}),
      onDraw: (r) => {
        // the export's issues and the viewer's model's as it says them (one of a source still opening held back); with
        // `sourceIssues`, what its sources said when they opened too (a file's lines left out)
        const issues = statusIssues(
          r.issues,
          view.sourceIssues ? [...viewer.modelIssues.value, ...viewer.store.issues.value] : viewer.modelIssues.value,
        );
        status.textContent = `${r.settled ? 'Settled' : 'Not settled'}; ${issues.length} issue${issues.length === 1 ? '' : 's'}${
          issues.length ? `: ${issues.map((i) => i.message).join('; ')}` : ''
        }`;
      },
    });
    onCleanup(() => svg.dispose());
    const find = async () => {
      if (!view.find) return;
      found.textContent = `Finding “${input.value}”…`;
      found.textContent = await view.find.run(viewer, input.value).catch((e) => `${e?.message ?? e}`);
    };
    if (view.find) {
      button.onclick = () => void find();
      input.onkeydown = (e) => {
        if (e.key === 'Enter') void find();
      };
    }
    void viewer
      .load(spec)
      .then(() => viewer.settle())
      .then(() => view.opening?.(viewer))
      .then(find)
      .then(() => svg.draw());
    return root;
  }
}

/** The page's own files by `fetch`: what the Files stories print a head of. */
const fetched = (uri: string) =>
  fetch(uri).then((r) => (r.ok ? r.arrayBuffer() : Promise.reject(new Error(`${r.status}`))));
const pageReader = {
  readText: (uri: string) => fetched(uri).then((b) => new TextDecoder().decode(b)),
  readBytes: (uri: string) => fetched(uri).then((b) => new Uint8Array(b)),
};

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
    out.push(section(file.format === 'vcf' ? `Header: ${file.uri}` : `File: ${file.uri}`, head));
    // read the way the viewer reads it: relative to the page, which serves the examples' fixtures (a recorded file from
    // its pages, a gzip file inflated; a VCF's header by its index)
    void fileText(pageReader, file, new AbortController().signal)
      .then((text) => {
        head.textContent = text;
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
 * source nodes and the head of each file they name are printed beside it; `config`: the viewer's items; `extensions`:
 * the viewer's extensions).
 */
export function exampleStory(
  view: ExampleView,
  options: {
    readonly navigators?: 'omit' | 'include';
    readonly files?: boolean;
    readonly sourceIssues?: boolean;
    readonly stacked?: boolean;
    readonly opening?: (viewer: SeqViewer) => unknown;
    readonly find?: ExampleModel['find'];
    readonly summary?: string;
    readonly config?: SeqViewerSpec['config'];
    readonly behaviors?: SeqViewerSpec['behaviors'];
    readonly extensions?: SeqViewerSpec['extensions'];
  } = {},
): HTMLElement {
  const el = document.createElement(ExampleRoot.tag) as ExampleRoot;
  el.model =
    options.navigators ||
    options.files ||
    options.sourceIssues ||
    options.stacked ||
    options.opening ||
    options.find ||
    options.summary ||
    options.config ||
    options.behaviors ||
    options.extensions
      ? { ...view, ...options }
      : view;
  return el;
}
