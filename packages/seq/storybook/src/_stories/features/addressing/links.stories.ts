/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type ExampleView, FeatureViews } from '@molstar/seqx-examples';
import { Passages } from '@molstar/seqx-examples/p53/addressing';
import { Button } from '@molstar/ui/components/button/button';
import { defineElements, onCleanup } from '@molstar/ui/dom';
import { UiElement } from '@molstar/ui/elements/ui-element';
import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { composeViewer, hoverLink, LinkStatus, leaveLink } from './compose.js';

// Links (addressing): a short text beside the p53 view whose links act on addresses in Mol*'s link grammar
// (`!seq-focus=/refs/dbd`): focus, highlight on a channel, select; one residue linked in two numberings. Hovering a link
// highlights what it names on the links' own channel. The highlight channels are drawn on the view, each a tint of its
// color, and listed below it with the selection and the focus. "Copy link" gives the link of the selection (else of what
// the pointer is over) as the shortest stable address. A link's text and href are put in as text and attributes, never
// markup; a click runs the link's actions instead of navigating.

const Link = 'color: var(--msui-accent); cursor: pointer; text-decoration: underline;';
const Mono = 'font: var(--msui-text-xs) / var(--msui-leading-xs) var(--msui-font-mono); color: var(--msui-text-muted);';

class LinksRoot extends UiElement<ExampleView> {
  static readonly tag = 'seqx-links-root';

  protected render(view: ExampleView): HTMLElement {
    const root = document.createElement('div');
    root.style.cssText =
      'display: grid; gap: var(--msui-space-4); padding: var(--msui-space-5); min-height: 100vh; background: var(--msui-bg); color: var(--msui-text); font: var(--msui-text-sm) / var(--msui-leading-sm) var(--msui-font-sans);';
    const text = document.createElement('div');
    text.style.cssText = 'display: grid; gap: var(--msui-space-2); max-width: 72ch;';
    const canvasHost = document.createElement('div');
    canvasHost.style.cssText = 'min-width: 0;';
    const tools = document.createElement('div');
    tools.style.cssText = 'display: flex; gap: var(--msui-space-2); align-items: baseline; flex-wrap: wrap;';
    const copy = Button({
      children: 'Copy link',
      onClick: () => {
        const link = viewer.addresses.copyLink();
        copied.textContent = link ?? 'Nothing selected or under the pointer';
        if (link) void navigator.clipboard?.writeText(link).catch(() => undefined);
      },
    }) as HTMLElement;
    const copied = document.createElement('code');
    copied.style.cssText = Mono;
    tools.append(copy, copied);
    const session = document.createElement('p');
    session.style.cssText = `margin: 0; ${Mono}`;
    const note = document.createElement('p');
    note.style.cssText = 'margin: 0; color: var(--msui-text-muted);';
    note.textContent =
      'Highlights (a link’s channel, and hovering a link) are drawn on the view as tints of their channel’s color, and listed here.';
    const status = document.createElement('p');
    status.style.cssText = 'margin: 0; color: var(--msui-text-muted);';
    root.append(text, canvasHost, tools, session, note, status);

    const viewer = composeViewer(view, { target: canvasHost });
    onCleanup(() => viewer.retire());
    // one message a link: a wait says it once, then `Focused` when the view has moved
    const links = new LinkStatus(viewer, (m) => {
      status.textContent = m;
    });
    onCleanup(() => links.dispose());

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
          status.textContent = links.follow(run.href);
        });
        a.addEventListener('pointerenter', () => hoverLink(viewer, run.href));
        a.addEventListener('pointerleave', () => leaveLink(viewer));
        p.append(a);
      }
      text.append(p);
    }

    // the session as addresses: what is selected, the highlight channels, and the focus
    const showSession = () => {
      const channels = viewer.highlight.names.value.map((n) => `${n}: ${viewer.highlight.get(n)?.address ?? ''}`);
      session.textContent = `selection: ${viewer.addresses.selectionAddress() ?? '—'} · focus: ${
        viewer.addresses.focusAddress() ?? '—'
      } · highlights: ${channels.length ? channels.join(', ') : '—'}`;
    };
    const offs = [
      viewer.selection.value.listen(showSession),
      viewer.highlight.names.listen(showSession),
      viewer.rested.listen(showSession),
    ];
    onCleanup(() => {
      for (const off of offs) off();
    });
    void viewer.load(view.build()).then(showSession);
    return root;
  }
}

defineElements([LinksRoot]);

const meta: Meta = { title: 'Features/Addressing/Links' };
export default meta;

/** A text whose links focus, highlight and select in the view; one residue in two numberings; copy link. */
export const Links: StoryObj = {
  render: () => {
    const el = document.createElement(LinksRoot.tag) as LinksRoot;
    el.model = FeatureViews.find((v) => v.id === 'addressing')!;
    return el;
  },
};
