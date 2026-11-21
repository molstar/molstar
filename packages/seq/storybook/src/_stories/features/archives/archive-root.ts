/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { archiveReader, readArchive } from '@molstar/seqx-core/data';
import type { ExampleView } from '@molstar/seqx-examples';
import { DefaultSeqViewerSpec } from '@molstar/seqx-viewer/default-spec';
import { createSeqViewer, type SeqViewer } from '@molstar/seqx-viewer/viewer';
import { Button } from '@molstar/ui/components/button/button';
import { defineElements, onCleanup } from '@molstar/ui/dom';
import { UiElement } from '@molstar/ui/elements/ui-element';
import { type Hostile, listArchive, pageAssets, printRefusal, saveView } from './archive-model.js';

// The Archives stories (examples-and-stories.md): a view saved as an archive and drawn from it alone, and archives the
// reader refuses. The element is a UI root, so the theme's tokens apply; the viewer is retired when Storybook removes
// the story.

type ArchiveModel =
  | { readonly kind: 'saved'; readonly view: ExampleView; readonly summary: string }
  | { readonly kind: 'refused'; readonly archives: readonly Hostile[]; readonly summary: string };

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
      const lead = document.createElement('p');
      lead.style.cssText = 'margin: 0; color: var(--msui-text-muted); max-width: 72ch;';
      lead.textContent = model.summary;
      root.append(lead);
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
    const tools = document.createElement('div');
    tools.style.cssText = 'display: flex; gap: var(--msui-space-2); align-items: baseline; flex-wrap: wrap;';
    const status = document.createElement('span');
    status.style.cssText = 'color: var(--msui-text-muted);';
    root.append(
      summary,
      tools,
      titled('The archive', listing),
      titled('Drawn from the archive alone', canvasHost),
    );
    let retired = false;
    let viewer: SeqViewer | undefined;
    let bytes: Uint8Array | undefined;
    onCleanup(() => {
      retired = true;
      viewer?.retire();
    });
    /** Draws the archive's spec with the archive as the only place files are read from. */
    const draw = async (archive: NonNullable<Awaited<ReturnType<typeof readArchive>>['archive']>) => {
      viewer?.retire();
      canvasHost.replaceChildren();
      const base = DefaultSeqViewerSpec();
      viewer = createSeqViewer(
        {
          ...base,
          registry: [...(base.registry ?? []), ...(model.view.registry?.() ?? [])],
          ports: { assets: () => ({ reader: archiveReader(archive) }) },
        },
        canvasHost,
      );
      await viewer.load(archive.spec);
      await viewer.settle();
    };
    // Save: the archive's bytes as a file the browser downloads; Open: a zip the reader picks, read and drawn instead
    const save = Button({
      children: 'Save the archive',
      onClick: () => {
        if (!bytes) return;
        const url = URL.createObjectURL(new Blob([bytes as BlobPart], { type: 'application/zip' }));
        const a = document.createElement('a');
        a.href = url;
        a.download = `${model.view.id}.seqx.zip`;
        a.click();
        URL.revokeObjectURL(url);
        status.textContent = `Saved ${a.download} (${bytes.byteLength} bytes).`;
      },
    }) as HTMLElement;
    const input = document.createElement('input');
    input.type = 'file';
    input.accept = '.zip,application/zip';
    input.style.display = 'none';
    input.addEventListener('change', () => {
      const file = input.files?.[0];
      if (!file) return;
      void (async () => {
        const picked = new Uint8Array(await file.arrayBuffer());
        const read = await readArchive(picked);
        if (retired) return;
        if (!read.archive) {
          status.textContent = `${file.name} was refused: ${read.issues.map((i) => i.message).join('; ')}`;
          return;
        }
        listing.textContent = `${file.name}\n${listArchive(picked, read.archive)}`;
        await draw(read.archive);
        status.textContent = `Opened ${file.name}: the view below is drawn from it alone.`;
      })();
    });
    const open = Button({ children: 'Open an archive…', onClick: () => input.click() }) as HTMLElement;
    tools.append(save, open, input, status);
    void (async () => {
      try {
        const saved = await saveView(model.view.build(), pageAssets, { keep: (u) => u.startsWith('generated:') });
        if (retired) return;
        bytes = saved.bytes;
        listing.textContent = saved.listing;
        await draw(saved.archive);
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
export function refusalStory(archives: readonly Hostile[], summary: string): HTMLElement {
  const el = document.createElement(ArchiveRoot.tag) as ArchiveRoot;
  el.model = { kind: 'refused', archives, summary };
  return el;
}
