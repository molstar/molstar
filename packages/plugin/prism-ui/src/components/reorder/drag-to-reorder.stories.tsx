/** @jsxImportSource @molstar/plugin-prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { map } from 'rxjs';
import { story } from '../../_stories/story-root.js';
import { For } from '../../dom/control-flow.js';
import { observableState } from '../../reactive/observable-state.js';
import { Icon } from '../icon/icon.js';
import { dragToReorder } from './drag-to-reorder.js';

const meta: Meta = { title: 'Components/Drag to reorder' };
export default meta;

const panel =
  'display: grid; gap: var(--msp-space-3); width: 300px; padding: var(--msp-space-4); background: var(--msp-color-panel); color: var(--msp-color-text); font-size: var(--msp-font-size-sm);';
const readout = 'color: var(--msp-color-text-muted);';
// The marks `dragToReorder` sets, styled as the snapshots panel does.
const css = `
.story-reorder-item { display: flex; align-items: center; gap: var(--msp-space-3); height: 30px; padding: 0 var(--msp-space-3); border-radius: var(--msp-radius-md); background: var(--msp-color-control); }
.story-reorder-item[data-dragging] { opacity: 0.5; }
.story-reorder-item[data-drop="before"] { box-shadow: 0 -2px 0 var(--msp-color-text-muted); }
.story-reorder-item[data-drop="after"] { box-shadow: 0 2px 0 var(--msp-color-text-muted); }
.story-reorder-grip { display: grid; place-items: center; width: 14px; height: 22px; color: var(--msp-color-text-subtle); cursor: grab; touch-action: none; }
`;

/** Five rows moved by their grip: the row being moved fades and a line shows where it lands. */
export const Rows: StoryObj = {
  render: () =>
    story(() => {
      const rows = observableState(['Overview', 'Binding site', 'Ligand', 'Surface', 'Final']);
      let list: HTMLElement | undefined;
      const moveTo = (id: string, index: number) => {
        const next = rows.value.filter((r) => r !== id);
        next.splice(index, 0, id);
        rows.set(next);
      };
      return (
        <div style={panel}>
          <style>{css}</style>
          <div
            style="display: grid; gap: var(--msp-space-2)"
            ref={(el: HTMLElement) => {
              list = el;
            }}
          >
            <For each={rows.asReadonly} key={(r) => r}>
              {(label) => (
                <div class="story-reorder-item" data-id={label}>
                  <span
                    class="story-reorder-grip"
                    aria-hidden="true"
                    onPointerDown={(e) => dragToReorder(e, list, label, (index) => moveTo(label, index))}
                  >
                    <Icon name="grip-vertical" size={14} />
                  </span>
                  <span>{label}</span>
                </div>
              )}
            </For>
          </div>
          <span style={readout}>order: {rows.asReadonly.pipe(map((r) => r.join(', ')))}</span>
        </div>
      );
    }),
};
