/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { StructureAddresses } from '@molstar/seqx-examples/p53/structure';
import type { ExampleView } from '@molstar/seqx-examples';
import { DefaultSeqViewerSpec } from '@molstar/seqx-viewer/default-spec';
import type { StructureBinding } from '@molstar/seqx-viewer/link-structure';
import { mountToolbar } from '@molstar/seqx-viewer/mount-toolbar';
import { createSeqViewer } from '@molstar/seqx-viewer/viewer';
import { Button } from '@molstar/ui/components/button/button';
import { defineElements, onCleanup } from '@molstar/ui/dom';
import { UiElement } from '@molstar/ui/elements/ui-element';
import { fakeStructureElement } from '../../support/fake-structure/fake-view.js';
import { ColorByEntry, colorByBehavior, colorByLayer } from '../../support/host-adapter.js';
import { chainTableOf, Datasets, type DatasetName, type FakeLink, linkFake } from './compose.js';

// The structure link's stories (structure-syncing): a Seq* view above the fake structure viewer, linked by
// `linkStructure` through the host's adapter. The fake is a small viewer of its own (chains as rows of residue cells, a
// grammar of its own, MolViewSpec states, an optional latency); the log under it lists the messages that cross. The
// element is a UI root, so the theme's tokens apply; the viewer and the link end when Storybook removes the story.

export interface StructureLinkStory {
  readonly view: ExampleView;
  readonly dataset: DatasetName;
  readonly hint: string;
  /** Milliseconds the fake delays what it tells (0). */
  readonly latency?: number;
  /** A toggle for the map the link goes through: the ranked one (SIFTS's) or the alignment's (`via`). */
  readonly chooseMap?: boolean;
  /** The host's header entry "Color the structure by this track" on every track. */
  readonly colorBy?: boolean;
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

/** The theme the root shows (the viewer's tokens read it the same way). */
const schemeOf = (el: Element): 'light' | 'dark' =>
  (el.closest('[data-msui-theme]')?.getAttribute('data-msui-theme') ||
    el.ownerDocument.documentElement.getAttribute('data-msui-theme')) === 'dark'
    ? 'dark'
    : 'light';

class StructureLinkRoot extends UiElement<StructureLinkStory> {
  static readonly tag = 'seqx-structure-link-root';

  protected render(story: StructureLinkStory): HTMLElement {
    const root = document.createElement('div');
    root.style.cssText =
      'display: grid; gap: var(--msui-space-5); padding: var(--msui-space-5); min-height: 100vh; align-content: start; background: var(--msui-bg); color: var(--msui-text); font: var(--msui-text-sm) / var(--msui-leading-sm) var(--msui-font-sans);';
    const hint = document.createElement('p');
    hint.style.cssText = 'margin: 0; color: var(--msui-text-muted); max-width: 72ch;';
    hint.textContent = story.hint;
    const canvasHost = document.createElement('div');
    canvasHost.style.cssText = 'min-width: 0;';
    const fakeHost = document.createElement('div');
    const status = document.createElement('p');
    status.style.cssText = 'margin: 0; color: var(--msui-text-muted);';
    const logEl = document.createElement('pre');
    logEl.style.cssText =
      'margin: 0; padding: var(--msui-space-3); min-height: 12em; overflow: auto; font: var(--msui-text-xs) / var(--msui-leading-xs) var(--msui-font-mono); background: var(--msui-surface, transparent); border: 1px solid var(--msui-line);';
    const tools = document.createElement('div');
    tools.style.cssText = 'display: flex; gap: var(--msui-space-2); align-items: baseline;';
    root.append(
      hint,
      tools,
      titled('Seq*', canvasHost),
      titled('The structure viewer (a fake)', fakeHost),
      titled('Messages', logEl),
      status,
    );

    const dataset = Datasets[story.dataset];
    const base = DefaultSeqViewerSpec();
    let linked: FakeLink | undefined;
    let through: 'ranked' | 'alignment' = 'ranked';
    const bindings = (): StructureBinding[] | undefined =>
      through === 'alignment' ? [{ complex: dataset.complex, via: [StructureAddresses.aln] }] : undefined;
    const viewer = createSeqViewer(
      {
        ...base,
        ui: { ...base.ui, toolbar: mountToolbar },
        registry: [
          ...(base.registry ?? []),
          ...(story.view.registry?.() ?? []),
          ...(story.colorBy ? [ColorByEntry] : []),
        ],
        behaviors: [
          ...base.behaviors,
          ...(story.colorBy
            ? [
                colorByBehavior((layer) =>
                  linked
                    ? colorByLayer(viewer, linked.fake, linked.adapter, layer, schemeOf(root), bindings())
                    : 'The structure is not linked yet',
                ),
              ]
            : []),
        ],
      },
      canvasHost,
    );
    onCleanup(() => {
      linked?.dispose();
      viewer.retire();
    });

    let unlog: (() => void) | undefined;
    const link = (table: Parameters<typeof linkFake>[1]) => {
      linked?.dispose();
      unlog?.();
      const l = (linked = linkFake(viewer, table, dataset.complex, {
        latency: story.latency ?? 0,
        ...(bindings() ? { link: { complexes: bindings()! } } : {}),
      }));
      const show = () => {
        logEl.textContent = l.log.value.join('\n');
      };
      show();
      unlog = l.log.listen(show);
      fakeHost.replaceChildren(fakeStructureElement(l.fake));
      status.textContent = `${table.entry}: ${l.adapter.facts.units.length} units in ${dataset.complex}; messages delayed ${story.latency ?? 0} ms; the link goes through ${through === 'ranked' ? 'the ranked map' : 'the alignment (via)'}.`;
    };

    if (story.chooseMap) {
      const toggle = Button({
        children: 'Link through: the ranked map (switch to the alignment)',
        onClick: () => {
          through = through === 'ranked' ? 'alignment' : 'ranked';
          (toggle as HTMLElement).textContent =
            through === 'ranked'
              ? 'Link through: the ranked map (switch to the alignment)'
              : 'Link through: the alignment (switch to the ranked map)';
          if (table) link(table);
        },
      }) as HTMLElement;
      tools.append(toggle);
    }

    let table: Parameters<typeof linkFake>[1] | undefined;
    void (async () => {
      const [text] = await Promise.all([
        fetch(dataset.tableFile).then((r) => (r.ok ? r.text() : Promise.reject(new Error(`${r.status}`)))),
        viewer.load(story.view.build()).then(() => viewer.settle()),
      ]);
      table = chainTableOf(text);
      link(table);
    })().catch((e) => {
      status.textContent = `The structure could not be linked: ${e instanceof Error ? e.message : e}`;
    });
    return root;
  }
}

defineElements([StructureLinkRoot]);

/** A Storybook `render` result: the view and the fake structure viewer linked. */
export function structureLinkStory(story: StructureLinkStory): HTMLElement {
  const el = document.createElement(StructureLinkRoot.tag) as StructureLinkRoot;
  el.model = story;
  return el;
}
