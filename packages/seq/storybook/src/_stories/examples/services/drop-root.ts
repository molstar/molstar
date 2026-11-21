/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { DefaultSeqViewerSpec } from '@molstar/seqx-viewer/default-spec';
import { mountToolbar } from '@molstar/seqx-viewer/mount-toolbar';
import { createSeqViewer, type SeqViewer } from '@molstar/seqx-viewer/viewer';
import { Button } from '@molstar/ui/components/button/button';
import { defineElements, onCleanup } from '@molstar/ui/dom';
import { UiElement } from '@molstar/ui/elements/ui-element';
import { DropRuntimes, droppedUri, specOfDropped } from '../../support/dropped-file.js';

// Your own file (examples.md, Showcases): drop a file onto the view, or pick one, and it is shown: `.seq`, FASTA, Stockholm,
// aligned FASTA, TSV, CSV, JSON, MessagePack, GFF3 and mmCIF files; GenBank says it is not read yet. The page's assets port
// answers the dropped bytes at one location and nothing else, so a file never leaves the page. The buttons load a sample
// of each format from the examples' fixtures, as a dropped file.

export interface DropStory {
  readonly hint: string;
}

const Samples: ReadonlyArray<readonly [name: string, file: string]> = [
  ['2OCJ.seq', 'fixtures/p53/2OCJ.seq'],
  ['1N8Z.fasta', 'fixtures/antibody/1N8Z.fasta'],
  ['PF00042-seed.sto', 'fixtures/globins/PF00042-seed.sto'],
  ['domains.tsv', 'fixtures/files/domains.tsv'],
  ['2ocj.cif.gz', 'fixtures/structures/2OCJ.cif.gz'],
  ['pUC19.gb', 'fixtures/structures/sources.json'],
];

class DropRoot extends UiElement<DropStory> {
  static readonly tag = 'seqx-drop-root';

  protected render(story: DropStory): HTMLElement {
    const root = document.createElement('div');
    root.style.cssText =
      'display: grid; gap: var(--msui-space-4); padding: var(--msui-space-5); min-height: 100vh; align-content: start; background: var(--msui-bg); color: var(--msui-text); font: var(--msui-text-sm) / var(--msui-leading-sm) var(--msui-font-sans);';
    const hint = document.createElement('p');
    hint.style.cssText = 'margin: 0; color: var(--msui-text-muted); max-width: 72ch;';
    hint.textContent = story.hint;
    const zone = document.createElement('div');
    zone.style.cssText =
      'display: grid; gap: var(--msui-space-3); padding: var(--msui-space-5); border: 2px dashed var(--msui-line); border-radius: var(--msui-radius-md, 6px); min-height: 12em; align-content: start;';
    const message = document.createElement('p');
    message.style.cssText = 'margin: 0; color: var(--msui-text-muted);';
    message.textContent = 'Drop a file here.';
    const canvasHost = document.createElement('div');
    canvasHost.style.cssText = 'min-width: 0;';
    const tools = document.createElement('div');
    tools.style.cssText = 'display: flex; gap: var(--msui-space-2); align-items: baseline; flex-wrap: wrap;';
    const input = document.createElement('input');
    input.type = 'file';
    input.setAttribute('aria-label', 'Pick a file to show');
    zone.append(message, tools, canvasHost);
    root.append(hint, zone);

    let viewer: SeqViewer | undefined;
    let gone = false;
    onCleanup(() => {
      gone = true;
      viewer?.retire();
    });

    /** Shows the bytes as the file `name`: a new viewer whose assets port answers them. */
    const show = async (name: string, bytes: Uint8Array) => {
      message.textContent = `Reading ${name}…`;
      viewer?.retire();
      viewer = undefined;
      canvasHost.replaceChildren();
      const dropped = await specOfDropped(name, bytes);
      if (gone) return;
      if ('refused' in dropped) {
        message.textContent = dropped.refused;
        return;
      }
      const base = DefaultSeqViewerSpec();
      const location = droppedUri(name);
      const answer = (uri: string) => {
        if (uri !== location) throw new Error(`Nothing is read but the dropped file (${uri})`);
        return bytes;
      };
      const v = (viewer = createSeqViewer(
        {
          ...base,
          ui: { ...base.ui, toolbar: mountToolbar },
          registry: [
            ...(base.registry ?? []),
            ...DropRuntimes.map((runtime) => ({ kind: 'format-runtime' as const, runtime })),
          ],
          ports: {
            ...base.ports,
            assets: () => ({
              reader: {
                readText: async (uri) => new TextDecoder().decode(answer(uri)),
                readBytes: async (uri) => answer(uri),
              },
            }),
          },
        },
        canvasHost,
      ));
      await v.load(dropped.spec);
      message.textContent = `${name} — ${dropped.summary}`;
    };

    const showFile = async (file: File) => show(file.name, new Uint8Array(await file.arrayBuffer()));
    zone.addEventListener('dragover', (e) => {
      e.preventDefault();
      zone.style.borderColor = 'var(--msui-accent)';
    });
    zone.addEventListener('dragleave', () => {
      zone.style.borderColor = 'var(--msui-line)';
    });
    zone.addEventListener('drop', (e) => {
      e.preventDefault();
      zone.style.borderColor = 'var(--msui-line)';
      const file = e.dataTransfer?.files?.[0];
      if (file) void showFile(file);
    });
    input.addEventListener('change', () => {
      const file = input.files?.[0];
      if (file) void showFile(file);
    });
    tools.append(input);
    for (const [name, file] of Samples) {
      const b = Button({
        children: `Try ${name}`,
        onClick: () => {
          void fetch(file)
            .then((r) => (r.ok ? r.arrayBuffer() : Promise.reject(new Error(`${file}: ${r.status}`))))
            .then((buf) => show(name, new Uint8Array(buf)))
            .catch((e) => {
              message.textContent = `Could not read the sample: ${e instanceof Error ? e.message : e}`;
            });
        },
      }) as HTMLElement;
      tools.append(b);
    }
    return root;
  }
}

defineElements([DropRoot]);

/** A Storybook `render` result: a view that takes a dropped file. */
export function dropStory(story: DropStory): HTMLElement {
  const el = document.createElement(DropRoot.tag) as DropRoot;
  el.model = story;
  return el;
}
