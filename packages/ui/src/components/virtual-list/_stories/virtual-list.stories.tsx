/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { Note, page, story } from '../../../_stories/support.js';
import { setText } from '../../../dom/index.js';
import { state } from '../../../state.js';
import { Icon } from '../../icon/icon.js';
import { RecycledList } from '../recycled-list.js';
import { VirtualList } from '../virtual-list.js';

const meta: Meta = { title: 'Components/Virtual list' };
export default meta;

interface Item {
  readonly id: number;
  readonly label: string;
}

const items = state<readonly Item[]>(Array.from({ length: 10_000 }, (_, i) => ({ id: i, label: `Residue ${i + 1}` })));

// A list on a surface (styles §5): rows `.msui-row` (24px), the icon in its slot on the inset, the text after it, the
// count's end on the inner padding. The skeleton behind (a fast scroll) has its bars on the same text column.
const box = 'display: flex; width: 280px; height: 360px;';

/** 10,000 keyed rows (`VirtualList`): only what shows (and a margin) is in the DOM. */
export const Keyed: StoryObj = {
  render: () =>
    story(() => (
      <div style={page}>
        <div class="msui-surface" style={box}>
          <VirtualList items={items.asReadonly} itemKey={(r) => r.id} edge aria-label="Residues">
            {(r) => (
              <div class="msui-row">
                <span class="msui-slot">
                  <Icon name="circle" />
                </span>
                <span>{r.label}</span>
                <span class="msui-num" style="color: var(--msui-text-muted)">
                  {String(r.id % 97)}
                </span>
              </div>
            )}
          </VirtualList>
        </div>
        <Note>Scroll fast: the skeleton's bars sit where the rows' text does.</Note>
      </div>
    )),
};

/** The same with recycled rows (`RecycledList`): a pool of slots re-pointed as the window moves. */
export const Recycled: StoryObj = {
  render: () =>
    story(() => (
      <div style={page}>
        <div class="msui-surface" style={box}>
          <RecycledList
            items={items.asReadonly}
            itemKey={(r) => r.id}
            edge
            aria-label="Residues"
            slot={() => {
              const text = document.createTextNode('');
              const count = document.createTextNode('');
              const el = (
                <div class="msui-row">
                  <span class="msui-slot">
                    <Icon name="circle" />
                  </span>
                  <span>{text}</span>
                  <span class="msui-num" style="color: var(--msui-text-muted)">
                    {count}
                  </span>
                </div>
              ) as HTMLElement;
              return {
                el,
                show(item: Item | undefined) {
                  if (!item) return;
                  setText(text, item.label);
                  setText(count, String(item.id % 97));
                },
              };
            }}
          />
        </div>
      </div>
    )),
};
