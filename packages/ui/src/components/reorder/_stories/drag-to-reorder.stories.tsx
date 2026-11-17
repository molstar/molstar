/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { derive, type State } from '@molstar/core/reactive';
import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { Note, page, story } from '../../../_stories/support.js';
import { For } from '../../../dom/index.js';
import { state } from '../../../state.js';
import { Icon } from '../../icon/icon.js';
import { dragToReorder, reorderRows } from '../drag-to-reorder.js';

const meta: Meta = { title: 'Components/Drag to reorder' };
export default meta;

// The marks the gesture sets, styled by the list: rows on the alignment primitives (styles §5), the grip in the icon
// slot, the row being moved faded, a 2px line where it lands.
const css = `
.story-reorder-list { width: 280px; }
.story-reorder-list > [data-id] { position: relative; }
.story-reorder-list > [data-id]:hover { background: var(--msui-control-hover); }
.story-reorder-list > [data-dragging] { opacity: 0.5; }
.story-reorder-list > [data-drop]::after { content: ""; position: absolute; inset-inline: var(--msui-inner); height: 2px; border-radius: 1px; background: var(--msui-accent); }
.story-reorder-list > [data-drop="before"]::after { top: -1px; }
.story-reorder-list > [data-drop="after"]::after { bottom: -1px; }
.story-reorder-grip { cursor: grab; touch-action: none; color: var(--msui-text-subtle); }
`;

const initial = ['Overview', 'Binding site', 'Ligand', 'Surface', 'Final'];

function moved(rows: readonly string[], id: string, index: number): string[] {
  const next = rows.filter((r) => r !== id);
  next.splice(index, 0, id);
  return next;
}

function Rows(props: { rows: State<readonly string[]>; gesture: 'attached' | 'press' }) {
  const { rows } = props;
  let list: HTMLElement | undefined;
  return (
    <div
      class="msui-surface story-reorder-list"
      ref={(el: HTMLElement) => {
        list = el;
        if (props.gesture === 'attached')
          reorderRows(el, '.story-reorder-grip', (id, i) => rows.set(moved(rows.value, id, i)));
      }}
    >
      <For each={rows.asReadonly} by={(r) => r}>
        {(label) => (
          <div class="msui-row" data-id={label}>
            <span
              class="msui-slot story-reorder-grip"
              aria-hidden="true"
              onPointerDown={
                props.gesture === 'press'
                  ? (e) => dragToReorder(e, list, label, (i) => rows.set(moved(rows.value, label, i)))
                  : undefined
              }
            >
              <Icon name="grip-vertical" />
            </span>
            <span>{label}</span>
          </div>
        )}
      </For>
    </div>
  );
}

/** Rows moved by their grip (`reorderRows`, attached once with the list): the row fades, a line shows where it lands. */
export const List: StoryObj = {
  render: () =>
    story(() => {
      const rows = state<readonly string[]>(initial);
      return (
        <div style={page}>
          <style>{css}</style>
          <Rows rows={rows} gesture="attached" />
          <Note>order: {derive(rows.asReadonly, (r) => r.join(', '))}</Note>
        </div>
      );
    }),
};

/** The same with `dragToReorder` called on each grip's press (outside a component). */
export const OnPress: StoryObj = {
  render: () =>
    story(() => {
      const rows = state<readonly string[]>(initial);
      return (
        <div style={page}>
          <style>{css}</style>
          <Rows rows={rows} gesture="press" />
          <Note>order: {derive(rows.asReadonly, (r) => r.join(', '))}</Note>
        </div>
      );
    }),
};
