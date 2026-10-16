/** @vitest-environment jsdom */
/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { scriptedModel } from '@molstar/ai/testing';
import { Scope } from '@molstar/prism-ui/dom/scope';
import type { ViewRenderer } from '@molstar/prism-ui/dom/view';
import { ShellModel, type ShellViewContext } from '@molstar/prism-ui/shell/shell';
import { createViewerAgent } from '@molstar/prism-viewer/agent/session';
import { DefaultPrismViewerSpec } from '@molstar/prism-viewer/default-spec';
import { createPrismViewerHeadless } from '@molstar/prism-viewer/headless';
import type { PrismViewer } from '@molstar/prism-viewer/viewer';
import { afterEach, describe, expect, it } from 'vitest';
import { InspectorPageId } from '../../inspector-access.js';
import { createAssistant } from '../ui.js';

// The assistant's Inspect rows in its tab open the Inspector, a page over the app (agent.md §1), at the view each names;
// its commands do the same.

const created: { dispose(): unknown }[] = [];
afterEach(async () => {
  for (const d of created.splice(0).reverse()) await d.dispose();
  document.body.innerHTML = '';
});

async function mount() {
  const viewer: PrismViewer = await createPrismViewerHeadless({ spec: DefaultPrismViewerSpec() });
  created.push(viewer);
  const model = scriptedModel([{ text: 'Hello.' }]);
  const assistant = await createAssistant({
    viewer,
    createAgent: (models) =>
      createViewerAgent({
        viewer,
        models: { ...models, adapterFor: (p) => (p === 'scripted' ? model.adapter : models.adapterFor(p)) },
      }),
  });
  created.push(assistant);
  const shell = new ShellModel({ navigation: viewer.managers.navigation });
  created.push({ dispose: shell.registerPages(assistant.ui.pages ?? []) });
  const ctx = { shell } as unknown as ShellViewContext;
  const at = (renderer: ViewRenderer<ShellViewContext>) => {
    const host = document.createElement('div');
    document.body.append(host);
    created.push({ dispose: renderer.mount(host, ctx, new Scope()) });
    return host;
  };
  return { assistant, shell, at, ctx };
}

describe('the Inspect rows in the Assistant tab', () => {
  it('list Trace, Tools and Limits and open the Inspector page there', async () => {
    const { assistant, shell, at } = await mount();
    const tab = assistant.ui.tabs!.find((t) => t.id === 'assistant.tab.assistant')!;
    const host = at(tab.render);
    const section = host.querySelector('section[aria-label=Inspect]')!;
    expect(section).not.toBeNull();
    const rows = [...section.querySelectorAll<HTMLButtonElement>('button')];
    expect(rows.map((r) => r.firstElementChild?.textContent)).toEqual(['Trace', 'Tools', 'Limits']);
    expect(rows[0].textContent).toContain('nothing yet');
    expect(rows[1].textContent).toMatch(/\d+ tools, \d+ skills/);

    rows[1].click();
    expect(shell.page.value?.page.id).toBe(InspectorPageId);
    expect(shell.page.value?.focus).toBe('tools');

    // The page at that view: the toolkit's groups; then Limits; the Trace's empty state before a turn.
    const page = shell.page.value!.page;
    const header = at(page.header!);
    const body = at(page.render);
    expect(header.textContent).toContain('Tools');
    expect(body.textContent).toContain('Core tools');
    rows[2].click();
    expect(shell.page.value?.focus).toBe('limits');
    expect(body.textContent).toContain('Turn budget');
    expect(body.textContent).toContain('Context window');
    rows[0].click();
    expect(body.textContent).toContain('The agent has not run yet');

    // The commands open it too.
    shell.closePage();
    expect(shell.page.value).toBeUndefined();
    const ids = (assistant.commands ?? []).map((c) => c.id);
    expect(ids).toEqual(['assistant.trace', 'assistant.tools', 'assistant.limits']);
    void assistant.commands![2].run({} as never, {} as never, new AbortController().signal);
    expect(shell.page.value).toMatchObject({ page: { id: InspectorPageId }, focus: 'limits' });
  });
});
