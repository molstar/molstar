/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type ExampleView, FeatureViews } from '@molstar/seqx-examples';
import { Button } from '@molstar/ui/components/button/button';
import { defineElements, onCleanup } from '@molstar/ui/dom';
import { UiElement } from '@molstar/ui/elements/ui-element';
import type { Meta, StoryObj } from '@storybook/web-components-vite';
import {
  bookmarkRows,
  bookmarkView,
  composeViewer,
  focusBookmark,
  LinkStatus,
} from '../../features/addressing/compose.js';
import { RegionText } from './compose.js';

// Bookmarks and links on the gene's region, read region by region: each bookmark is a ref; the transcripts' arrive with
// their part of the region, so their bookmarks wait, then resolve. A link in the text moves the view the same way (a
// focus asked meanwhile runs when its node arrives); "Bookmark this view" adds a ref of what is shown. Every string from
// the spec or the text is put in as text, never markup.

const Muted =
  'color: var(--msui-text-muted); font: var(--msui-text-xs) / var(--msui-leading-xs) var(--msui-font-mono);';

class BookmarksAndLinksRoot extends UiElement<ExampleView> {
  static readonly tag = 'seqx-p53-bookmarks-and-links-root';

  protected render(view: ExampleView): HTMLElement {
    const root = document.createElement('div');
    root.style.cssText =
      'display: grid; align-content: start; gap: var(--msui-space-4); padding: var(--msui-space-5); min-height: 100vh; background: var(--msui-bg); color: var(--msui-text); font: var(--msui-text-sm) / var(--msui-leading-sm) var(--msui-font-sans);';
    const hint = document.createElement('p');
    hint.style.cssText = 'margin: 0; color: var(--msui-text-muted); max-width: 72ch;';
    hint.textContent =
      'Click a link in the text or a bookmark right away: the transcripts’ wait for their data, then focus the view. Bookmark the view to add one of your own.';
    const text = document.createElement('p');
    text.style.cssText = 'margin: 0; max-width: 72ch;';
    const canvasHost = document.createElement('div');
    canvasHost.style.cssText = 'min-width: 0;';
    const list = document.createElement('ul');
    list.style.cssText = 'margin: 0; padding: 0; list-style: none; display: grid; gap: var(--msui-space-2);';
    const status = document.createElement('p');
    status.style.cssText = 'margin: 0; color: var(--msui-text-muted);';
    const add = Button({
      variant: 'primary',
      children: 'Bookmark this view',
      onClick: () => {
        void Promise.resolve(bookmarkView(viewer)).then((m) => {
          status.textContent = m;
        });
      },
    }) as HTMLElement;
    const tools = document.createElement('div');
    tools.append(add);
    root.append(hint, text, canvasHost, tools, list, status);

    const viewer = composeViewer(view, { target: canvasHost });
    onCleanup(() => viewer.retire());
    // a link on a node that arrives later says it waits, and `Focused` when the view has moved to it
    const links = new LinkStatus(viewer, (m) => {
      status.textContent = m;
    });
    onCleanup(() => links.dispose());

    for (const run of RegionText) {
      if (typeof run === 'string') {
        text.append(run);
        continue;
      }
      const a = document.createElement('a');
      a.style.cssText = 'color: var(--msui-accent); cursor: pointer; text-decoration: underline;';
      a.textContent = run.text;
      a.setAttribute('href', run.href);
      a.title = run.href;
      a.addEventListener('click', (e) => {
        e.preventDefault();
        status.textContent = links.follow(run.href);
      });
      text.append(a);
    }

    const draw = () => {
      list.replaceChildren(
        ...bookmarkRows(viewer).map((r) => {
          const li = document.createElement('li');
          li.style.cssText = 'display: flex; gap: var(--msui-space-2); align-items: baseline; flex-wrap: wrap;';
          const b = Button({
            size: 'sm',
            children: r.label,
            onClick: () => {
              status.textContent = focusBookmark(viewer, r.address);
            },
          }) as HTMLElement;
          const state = document.createElement('span');
          state.style.cssText = Muted;
          state.textContent = r.status === 'waiting' ? 'waiting for its data' : r.status;
          const to = document.createElement('code');
          to.style.cssText = Muted;
          to.textContent = `${r.address} → ${r.to}`;
          li.append(b, state, to);
          return li;
        }),
      );
    };
    onCleanup(viewer.addresses.refs.listen(draw));
    void viewer
      .load(view.build())
      .then(() => viewer.settle())
      .then(() => {
        draw();
      })
      .catch((e: unknown) => {
        status.textContent = `The view did not load: ${e instanceof Error ? e.message : String(e)}`;
      });
    return root;
  }
}

defineElements([BookmarksAndLinksRoot]);

const meta: Meta = { title: 'Examples/p53/Bookmarks and links' };
export default meta;

/** A bookmark on a transcript that arrives later waits, then resolves; a link in the text moves the view. */
export const BookmarksAndLinks: StoryObj = {
  render: () => {
    const el = document.createElement(BookmarksAndLinksRoot.tag) as BookmarksAndLinksRoot;
    el.model = FeatureViews.find((v) => v.id === 'bookmarks')!;
    return el;
  },
};
