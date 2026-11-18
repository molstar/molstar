/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type ExampleView, FeatureViews } from '@molstar/seqx-examples';
import { defineElements, onCleanup } from '@molstar/ui/dom';
import { UiElement } from '@molstar/ui/elements/ui-element';
import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { bookmarkRows, bookmarkView, composeViewer, focusBookmark } from './compose.js';

// Bookmarks (addressing): the gene's region read region by region, with refs as bookmarks. Each bookmark focuses the
// view; one on a transcript waits until its region's answer brings the transcript, then resolves (a focus asked
// meanwhile runs then). "Bookmark this view" adds a ref of what is shown: one patch. Every string from the spec is put
// in as text, never markup.

const Row = 'display: flex; gap: var(--msui-space-2); align-items: baseline; flex-wrap: wrap;';
const Muted =
  'color: var(--msui-text-muted); font: var(--msui-text-xs) / var(--msui-leading-xs) var(--msui-font-mono);';

class BookmarksRoot extends UiElement<ExampleView> {
  static readonly tag = 'seqx-bookmarks-root';

  protected render(view: ExampleView): HTMLElement {
    const root = document.createElement('div');
    root.style.cssText =
      'display: grid; gap: var(--msui-space-4); padding: var(--msui-space-5); min-height: 100vh; background: var(--msui-bg); color: var(--msui-text); font: var(--msui-text-sm) / var(--msui-leading-sm) var(--msui-font-sans);';
    const summary = document.createElement('p');
    summary.style.cssText = 'margin: 0; color: var(--msui-text-muted); max-width: 72ch;';
    summary.textContent = view.summary;
    const canvasHost = document.createElement('div');
    canvasHost.style.cssText = 'min-width: 0;';
    const list = document.createElement('ul');
    list.style.cssText = 'margin: 0; padding: 0; list-style: none; display: grid; gap: var(--msui-space-2);';
    const add = document.createElement('button');
    add.textContent = 'Bookmark this view';
    const status = document.createElement('p');
    status.style.cssText = 'margin: 0; color: var(--msui-text-muted);';
    root.append(summary, canvasHost, add, list, status);

    const viewer = composeViewer(view, { target: canvasHost });
    onCleanup(() => viewer.retire());

    const draw = () => {
      list.replaceChildren(
        ...bookmarkRows(viewer).map((r) => {
          const li = document.createElement('li');
          li.style.cssText = Row;
          const b = document.createElement('button');
          b.textContent = r.label;
          if (r.color) b.style.borderColor = r.color;
          b.addEventListener('click', () => {
            status.textContent = focusBookmark(viewer, r.address);
          });
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
    add.addEventListener('click', () => {
      status.textContent = bookmarkView(viewer);
    });
    const off = viewer.addresses.refs.listen(draw);
    onCleanup(off);
    void viewer.load(view.build()).then(draw);
    return root;
  }
}

defineElements([BookmarksRoot]);

const meta: Meta = { title: 'Features/Addressing/Bookmarks' };
export default meta;

/** Refs as bookmarks: each focuses the view; one on a transcript waits for its region; add one of the view. */
export const Bookmarks: StoryObj = {
  render: () => {
    const el = document.createElement(BookmarksRoot.tag) as BookmarksRoot;
    el.model = FeatureViews.find((v) => v.id === 'bookmarks')!;
    return el;
  },
};
