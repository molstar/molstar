/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { archiveReader } from '@molstar/seqx-core/data';
import type { ExampleView } from '@molstar/seqx-examples';
import { DefaultSeqViewerSpec } from '@molstar/seqx-viewer/default-spec';
import { createSeqViewer, type SeqViewer } from '@molstar/seqx-viewer/viewer';
import { defineElements, onCleanup } from '@molstar/ui/dom';
import { UiElement } from '@molstar/ui/elements/ui-element';
import { type Hostile, pageAssets, printRefusal, saveView } from './archive-model.js';

// The Archives stories (examples-and-stories.md): a view saved as an archive and drawn from it alone, and archives the
// reader refuses. The element is a UI root, so the theme's tokens apply; the viewer is retired when Storybook removes
// the story.

type ArchiveModel =
  | { readonly kind: 'saved'; readonly view: ExampleView; readonly summary: string }
  | { readonly kind: 'refused'; readonly archives: readonly Hostile[] };

const Pre =
  'margin: 0; padding: var(--msui-space-3); overflow: auto; max-height: 360px; font: var(--msui-text-xs) / var(--msui-leading-xs) var(--msui-font-mono); background: var(--msui-surface, transparent); border: 1px solid var(--msui-line);';

function titled(title: string, body: HTMLElement): HTMLElement {
  const el = document.createElement('div');
  el.style.cssText = 'display: grid; gap: var(--msui-space-2); min-width: 0; align-content: start;';
  const h = document.createElement('h4');
  h.style.cssText =
    'margin: 0; font-size: var(--msui-text-sm); font-weight: var(--msui-weight-strong); color: var(--msui-text-muted);';
  h.textContent = title;
  el.append(h, body);
  return el;
}

const pre = (text: string) => {
  const el = document.createElement('pre');
  el.style.cssText = Pre;
  el.textContent = text;
  return el;
};

class ArchiveRoot extends UiElement<ArchiveModel> {
  static readonly tag = 'seqx-archive-root';

  protected render(model: ArchiveModel): HTMLElement {
    const root = document.createElement('div');
    root.style.cssText =
      'display: grid; gap: var(--msui-space-5); padding: var(--msui-space-5); min-height: 100vh; background: var(--msui-bg); color: var(--msui-text); font: var(--msui-text-sm) / var(--msui-leading-sm) var(--msui-font-sans);';
    if (model.kind === 'refused') {
      for (const h of model.archives) {
        const out = pre('Reading…');
        root.append(titled(h.id, out));
        void printRefusal(h).then((r) => {
          out.textContent = r.text;
        });
      }
      return root;
    }
    const summary = document.createElement('p');
    summary.style.cssText = 'margin: 0; color: var(--msui-text-muted); max-width: 72ch;';
    summary.textContent = model.summary;
    const listing = pre('Saving…');
    const canvasHost = document.createElement('div');
    canvasHost.style.cssText = 'min-width: 0;';
    root.append(summary, titled('The archive', listing), titled('Drawn from the archive alone', canvasHost));
    let retired = false;
    let viewer: SeqViewer | undefined;
    onCleanup(() => {
      retired = true;
      viewer?.retire();
    });
    void (async () => {
      try {
        const saved = await saveView(model.view.build(), pageAssets, { keep: (u) => u.startsWith('generated:') });
        if (retired) return;
        listing.textContent = saved.listing;
        const base = DefaultSeqViewerSpec();
        // the archive is the only place the files are read from
        viewer = createSeqViewer(
          {
            ...base,
            registry: [...(base.registry ?? []), ...(model.view.registry?.() ?? [])],
            ports: { assets: () => ({ reader: archiveReader(saved.archive) }) },
          },
          canvasHost,
        );
        await viewer.load(saved.archive.spec);
        await viewer.settle();
      } catch (e) {
        listing.textContent = `The archive could not be saved: ${e instanceof Error ? e.message : e}`;
      }
    })();
    return root;
  }
}

defineElements([ArchiveRoot]);

/** A Storybook `render` result: the view saved as an archive and drawn from it. */
export function archiveStory(view: ExampleView, summary: string): HTMLElement {
  const el = document.createElement(ArchiveRoot.tag) as ArchiveRoot;
  el.model = { kind: 'saved', view, summary };
  return el;
}

/** A Storybook `render` result: archives the reader refuses, each with what it says. */
export function refusalStory(archives: readonly Hostile[]): HTMLElement {
  const el = document.createElement(ArchiveRoot.tag) as ArchiveRoot;
  el.model = { kind: 'refused', archives };
  return el;
}
