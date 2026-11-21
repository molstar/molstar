/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Passages } from '@molstar/seqx-examples/services/addressing';
import { Button } from '@molstar/ui/components/button/button';
import { defineElements, onCleanup } from '@molstar/ui/dom';
import { UiElement } from '@molstar/ui/elements/ui-element';
import { hoverLink, leaveLink, LinkStatus } from '../../features/addressing/compose.js';
import { chainTableOf, litLine } from '../../features/structure-link/compose.js';
import { fakeStructureElement } from '../../support/fake-structure/fake-view.js';
import { type Scenario, Scenarios, ServicesApp } from '../../support/services-app.js';
import type { ServiceMode } from '../../support/services-port.js';

// The services stories' page (examples.md, Showcases): an application composing UniProt, PDBe (SIFTS) and AlphaFold DB,
// the Seq* view above the fake structure viewer, linked once SIFTS and the structure are both here. The switch changes
// where the answers come from, never the spec: the recorded files (after a latency, so the order of arrival is seen), or
// the services themselves through the asset policy's hosts. The page lists what was asked and from where. The element is a
// UI root, so the theme's tokens apply; the viewer, the link and the fake end when Storybook removes the story.

export interface ServicesStory {
  readonly scenario: Scenario['id'];
  readonly hint: string;
  /** A recorded answer's latency, ms: long enough to see the view fill in (600). */
  readonly delayMs?: number;
  /** The addressing text beside the view. */
  readonly text?: boolean;
}

const Box = 'display: grid; gap: var(--msui-space-2); min-width: 0; align-content: start;';
const Title =
  'margin: 0; font-size: var(--msui-text-sm); font-weight: var(--msui-weight-strong); color: var(--msui-text-muted);';
const Mono = 'font: var(--msui-text-xs) / var(--msui-leading-xs) var(--msui-font-mono);';
const Link = 'color: var(--msui-accent); cursor: pointer; text-decoration: underline;';

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

const fetchBytes = async (file: string) => {
  const r = await fetch(file);
  if (!r.ok) throw new Error(`${file}: ${r.status}`);
  return new Uint8Array(await r.arrayBuffer());
};

class ServicesRoot extends UiElement<ServicesStory> {
  static readonly tag = 'seqx-services-root';

  protected render(story: ServicesStory): HTMLElement {
    const scenario = Scenarios[story.scenario];
    const root = document.createElement('div');
    root.style.cssText =
      'display: grid; gap: var(--msui-space-4); padding: var(--msui-space-5); min-height: 100vh; align-content: start; background: var(--msui-bg); color: var(--msui-text); font: var(--msui-text-sm) / var(--msui-leading-sm) var(--msui-font-sans);';
    const hint = document.createElement('p');
    hint.style.cssText = 'margin: 0; color: var(--msui-text-muted); max-width: 72ch;';
    hint.textContent = story.hint;
    const tools = document.createElement('div');
    tools.style.cssText = 'display: flex; gap: var(--msui-space-3); align-items: baseline; flex-wrap: wrap;';
    const status = document.createElement('span');
    status.style.cssText = 'color: var(--msui-text-muted);';
    const text = document.createElement('div');
    text.style.cssText = 'display: grid; gap: var(--msui-space-2); max-width: 72ch;';
    const copied = document.createElement('code');
    copied.style.cssText = `${Mono} color: var(--msui-text-muted);`;
    const canvasHost = document.createElement('div');
    canvasHost.style.cssText = 'min-width: 0;';
    const fakeHost = document.createElement('div');
    const lit = document.createElement('p');
    lit.style.cssText = `margin: 0; ${Mono}`;
    const pre = (min: string) => {
      const el = document.createElement('pre');
      el.style.cssText = `margin: 0; padding: var(--msui-space-3); min-height: ${min}; overflow: auto; ${Mono} background: var(--msui-surface, transparent); border: 1px solid var(--msui-line);`;
      return el;
    };
    const readsEl = pre('6em');
    const logEl = pre('8em');
    const linkStatus = document.createElement('p');
    linkStatus.style.cssText = 'margin: 0; color: var(--msui-text-muted);';
    root.append(
      hint,
      tools,
      ...(story.text ? [text] : []),
      titled('Seq*', canvasHost),
      titled('The structure viewer (a fake)', lit, fakeHost),
      titled('What the services were asked', readsEl),
      titled('Messages between the view and the structure', logEl),
      linkStatus,
    );

    let gone = false;
    let app: ServicesApp | undefined;
    let links: LinkStatus | undefined;
    const offs: (() => void)[] = [];
    const rebind: (() => void)[] = [];

    const showStatus = () => {
      if (!app) return;
      const sifts = app.sifts.value.length ? `; SIFTS: ${app.sifts.value.length} chain mappings` : '';
      status.textContent = `${app.mode.value === 'fixtures' ? 'Recorded answers' : 'The live services'} · ${app.phase.value}${app.message.value ? `: ${app.message.value}` : ''}${sifts}`;
    };
    const showReads = () => {
      if (!app) return;
      readsEl.textContent =
        app.reads.value
          .map((r) =>
            r.mode === 'fixtures'
              ? `${r.service ?? 'file'}: ${r.uri}\n   recorded: ${r.file ?? '—'}`
              : `${r.service ?? '?'}: ${r.uri}\n   live`,
          )
          .join('\n') || 'nothing asked yet';
    };

    // what belongs to one viewer and one link: bound again each time the application makes another
    const bind = () => {
      for (const off of rebind.splice(0)) off();
      links?.dispose();
      links = undefined;
      const a = app;
      if (!a?.viewer) {
        fakeHost.replaceChildren();
        logEl.textContent = '';
        return;
      }
      const viewer = a.viewer;
      links = new LinkStatus(viewer, (m) => {
        linkStatus.textContent = m;
      });
      const link = a.link;
      if (link) {
        fakeHost.replaceChildren(fakeStructureElement(link.fake));
        const showLog = () => {
          logEl.textContent = link.log.value.join('\n');
        };
        const showLit = () => {
          lit.textContent = litLine(link.fake);
        };
        showLog();
        showLit();
        rebind.push(link.log.listen(showLog), link.fake.highlight.listen(showLit));
      } else {
        fakeHost.replaceChildren();
        lit.textContent = 'The structure is not linked yet';
        logEl.textContent = '';
      }
    };

    onCleanup(() => {
      gone = true;
      for (const off of offs.splice(0)) off();
      for (const off of rebind.splice(0)) off();
      links?.dispose();
      app?.dispose();
    });

    const toggle = Button({
      children: 'Switch to the live services',
      onClick: () => {
        if (!app) return;
        const next: ServiceMode = app.mode.value === 'fixtures' ? 'live' : 'fixtures';
        (toggle as HTMLElement).textContent =
          next === 'fixtures' ? 'Switch to the live services' : 'Back to the recorded answers';
        canvasHost.replaceChildren();
        void app.setMode(next);
      },
    }) as HTMLElement;
    tools.append(toggle, status);

    if (story.text) {
      const copy = Button({
        children: 'Copy link',
        onClick: () => {
          const link = app?.viewer?.addresses.copyLink();
          copied.textContent = link ?? 'Nothing selected or under the pointer';
          if (link) void navigator.clipboard?.writeText(link).catch(() => undefined);
        },
      }) as HTMLElement;
      tools.append(copy, copied);
      for (const passage of Passages) {
        const p = document.createElement('p');
        p.style.cssText = 'margin: 0;';
        for (const run of passage) {
          if (typeof run === 'string') {
            p.append(run);
            continue;
          }
          const a = document.createElement('a');
          a.style.cssText = Link;
          a.textContent = run.text;
          a.setAttribute('href', run.href);
          a.title = run.href;
          a.addEventListener('click', (e) => {
            e.preventDefault();
            if (links) linkStatus.textContent = links.follow(run.href);
          });
          a.addEventListener('pointerenter', () => app?.viewer && hoverLink(app.viewer, run.href));
          a.addEventListener('pointerleave', () => app?.viewer && leaveLink(app.viewer));
          p.append(a);
        }
        text.append(p);
      }
    }

    void (async () => {
      const table = chainTableOf(new TextDecoder().decode(await fetchBytes(scenario.tableFile)));
      if (gone) return;
      const a = (app = new ServicesApp({
        scenario,
        table,
        fixture: fetchBytes,
        delayMs: story.delayMs ?? 600,
        target: canvasHost,
        scheme: () => schemeOf(root),
      }));
      offs.push(
        a.phase.listen(showStatus),
        a.message.listen(showStatus),
        a.sifts.listen(showStatus),
        a.mode.listen(showStatus),
        a.reads.listen(showReads),
        a.objects.listen(bind),
      );
      showStatus();
      showReads();
      await a.start();
    })().catch((e) => {
      status.textContent = `The story could not start: ${e instanceof Error ? e.message : e}`;
    });
    return root;
  }
}

defineElements([ServicesRoot]);

/** A Storybook `render` result: the application composing the services, a scenario of it. */
export function servicesStory(story: ServicesStory): HTMLElement {
  const el = document.createElement(ServicesRoot.tag) as ServicesRoot;
  el.model = story;
  return el;
}
