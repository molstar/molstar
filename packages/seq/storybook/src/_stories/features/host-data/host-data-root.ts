/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { AssetsPort } from '@molstar/seqx-core/context';
import { DefaultSeqViewerSpec } from '@molstar/seqx-viewer/default-spec';
import { mountToolbar } from '@molstar/seqx-viewer/mount-toolbar';
import { createSeqViewer, type SeqViewer } from '@molstar/seqx-viewer/viewer';
import { Button } from '@molstar/ui/components/button/button';
import { defineElements, onCleanup } from '@molstar/ui/dom';
import { UiElement } from '@molstar/ui/elements/ui-element';
import { fakeStructureElement } from '../../support/fake-structure/fake-view.js';
import { FakeStructureViewer } from '../../support/fake-structure/fake-viewer.js';
import {
  heldLines,
  BFactorLegend,
  issueLines,
  hostDataSpec,
  hostTableOf,
  HostFile,
  HostTableFile,
  withDefaultMmcif,
  withFakeHost,
} from './compose.js';

// The Host data story (seq-data-ports, The fake host's adapter): one spec naming PDB 2OCJ by its file, a Seq* view of two
// of its chains' B-factors above the fake structure viewer. The fake holds 2OCJ (its recorded chain table) and Seq* is
// given its adapter: the nodes and the B-factors come from the fake's memory, no file is read, and "Drop chain" in the
// fake is followed by Seq* a frame later (the source reopened in place, the nodes of the other chains kept). The button
// beside it composes the same spec with the default adapter, which reads the file through the assets port: the same
// nodes. The element is a UI root, so the theme's tokens apply; the viewer ends when Storybook removes the story.

export interface HostDataStory {
  readonly hint: string;
}

const Box = 'display: grid; gap: var(--msui-space-2); min-width: 0; align-content: start;';
const Title =
  'margin: 0; font-size: var(--msui-text-sm); font-weight: var(--msui-weight-strong); color: var(--msui-text-muted);';

function titled(title: string, ...body: HTMLElement[]): HTMLElement {
  const el = document.createElement('div');
  el.style.cssText = Box;
  const h = document.createElement('h4');
  h.style.cssText = Title;
  h.textContent = title;
  el.append(h, ...body);
  return el;
}

class HostDataRoot extends UiElement<HostDataStory> {
  static readonly tag = 'seqx-host-data-root';

  protected render(story: HostDataStory): HTMLElement {
    const root = document.createElement('div');
    root.style.cssText =
      'display: grid; gap: var(--msui-space-5); padding: var(--msui-space-5); min-height: 100vh; align-content: start; background: var(--msui-bg); color: var(--msui-text); font: var(--msui-text-sm) / var(--msui-leading-sm) var(--msui-font-sans);';
    const hint = document.createElement('p');
    hint.style.cssText = 'margin: 0; color: var(--msui-text-muted); max-width: 72ch;';
    hint.textContent = story.hint;
    const canvasHost = document.createElement('div');
    canvasHost.style.cssText = 'min-width: 0;';
    const legend = document.createElement('p');
    legend.style.cssText = 'margin: 0; color: var(--msui-text-muted); font-size: var(--msui-text-xs); max-width: 72ch;';
    legend.textContent = BFactorLegend;
    const fakeHost = document.createElement('div');
    const held = document.createElement('pre');
    held.style.cssText =
      'margin: 0; padding: var(--msui-space-3); min-height: 8em; overflow: auto; font: var(--msui-text-xs) / var(--msui-leading-xs) var(--msui-font-mono); background: var(--msui-surface, transparent); border: 1px solid var(--msui-line);';
    const tools = document.createElement('div');
    tools.style.cssText = 'display: flex; gap: var(--msui-space-2); align-items: baseline;';
    root.append(
      hint,
      tools,
      titled('Seq*', canvasHost, legend),
      titled('The structure viewer (a fake)', fakeHost),
      titled('What Seq* holds of PDB 2OCJ', held),
    );

    let gone = false;
    let viewer: SeqViewer | undefined;
    let fake: FakeStructureViewer | undefined;
    let by: 'fake' | 'default' = 'fake';
    let reads: string[] = [];
    // what ends with a mounting (the viewer's listeners, the fake's edit listener)
    const offs: (() => void)[] = [];
    let table: ReturnType<typeof hostTableOf> | undefined;
    let edits = 0;

    const show = () => {
      if (!viewer) return;
      held.textContent = [
        by === 'fake'
          ? 'Read by: the fake structure viewer’s adapter, from its chain table'
          : 'Read by: the default adapter, the mmCIF runtime, from the file',
        `Assets read: ${reads.length}${reads.length ? ` (${reads.join(', ')})` : ''}`,
        `Edits the fake told: ${edits}`,
        '',
        ...heldLines(viewer),
        '',
        ...issueLines(viewer),
      ].join('\n');
    };

    // the assets port: the page's own reads, counted (what the fake's adapter must never make)
    const assets = (): AssetsPort => {
      const read = async (uri: string) => {
        reads.push(uri);
        show();
        const r = await fetch(new URL(uri, location.href));
        if (!r.ok) throw new Error(`${r.status}`);
        return r;
      };
      return {
        reader: {
          readText: async (uri) => (await read(uri)).text(),
          readBytes: async (uri) => new Uint8Array(await (await read(uri)).arrayBuffer()),
        },
      };
    };

    // a mounting starts from the table as recorded: a new fake, so an edit made under one adapter never meets the other
    const mount = async (next: 'fake' | 'default') => {
      for (const off of offs.splice(0)) off();
      viewer?.retire();
      fake?.dispose();
      by = next;
      reads = [];
      edits = 0;
      if (!table) return;
      const f = (fake = new FakeStructureViewer(table));
      fakeHost.replaceChildren(fakeStructureElement(f));
      offs.push(
        f.onEdit(() => {
          edits++;
          show();
        }),
      );
      const base = DefaultSeqViewerSpec();
      const v = (viewer = createSeqViewer(
        {
          ...base,
          ui: { ...base.ui, toolbar: mountToolbar },
          registry: by === 'fake' ? withFakeHost(base.registry ?? [], f) : withDefaultMmcif(base.registry ?? []),
          ports: { ...base.ports, assets },
        },
        canvasHost,
      ));
      offs.push(v.store.store.nodeVersion.listen(show), v.modelIssues.listen(show), v.layoutIssues.listen(show));
      show();
      await v.load(hostDataSpec());
      show();
    };

    onCleanup(() => {
      gone = true;
      for (const off of offs.splice(0)) off();
      viewer?.retire();
      fake?.dispose();
    });

    const toggle = Button({
      children: 'Read the file instead (the default adapter; resets the fake)',
      onClick: () => {
        const next = by === 'fake' ? 'default' : 'fake';
        (toggle as HTMLElement).textContent =
          next === 'fake'
            ? 'Read the file instead (the default adapter; resets the fake)'
            : 'Ask the fake again (its adapter; resets the fake)';
        void mount(next);
      },
    }) as HTMLElement;
    tools.append(toggle);

    void (async () => {
      const text = await fetch(HostTableFile).then((r) => (r.ok ? r.text() : Promise.reject(new Error(`${r.status}`))));
      if (gone) return;
      table = hostTableOf(text);
      await mount('fake');
    })().catch((e) => {
      held.textContent = `The story could not start (${HostFile}): ${e instanceof Error ? e.message : e}`;
    });
    return root;
  }
}

defineElements([HostDataRoot]);

/** A Storybook `render` result: the view above the fake structure viewer, served by its adapter. */
export function hostDataStory(story: HostDataStory): HTMLElement {
  const el = document.createElement(HostDataRoot.tag) as HostDataRoot;
  el.model = story;
  return el;
}
