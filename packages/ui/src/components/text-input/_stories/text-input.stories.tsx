/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { story } from '../../../_stories/support.js';
import { state } from '../../../state.js';
import { TextInput } from '../text-input.js';

const meta: Meta = { title: 'Components/Text input' };
export default meta;

// A surface with a form (styles §5): labels on the inset, fields from the values' column to the inner padding.
const panel = 'width: 320px; margin: var(--msui-space-4);';

/** In a form (`size="sm"`): the field is a control in its row, its text where a plain value's starts; the value follows the typing. */
export const InAForm: StoryObj = {
  render: () =>
    story(() => {
      const id = state('1tqn');
      const key = state('');
      return (
        <div class="msui-surface" style={panel}>
          <div class="msui-form">
            <span class="msui-label">PDB ID</span>
            <TextInput size="sm" aria-label="PDB ID" placeholder="e.g. 1tqn" value={id} onInput={(v) => id.set(v)} />
            <span class="msui-label">Value</span>
            <span class="msui-value">{id}</span>
            <span class="msui-label">API key</span>
            <TextInput size="sm" type="password" aria-label="API key" value={key} onInput={(v) => key.set(v)} />
            <span class="msui-label">Disabled</span>
            <TextInput size="sm" aria-label="Disabled" value="read only" disabled />
          </div>
        </div>
      );
    }),
};

/** On its own (the default size, a row's height: 24px), a search field: Escape empties it; Enter submits (echoed below). */
export const Search: StoryObj = {
  render: () =>
    story(() => {
      const query = state('');
      const submitted = state('');
      return (
        <div class="msui-surface" style={panel}>
          <div style="display: grid; padding: var(--msui-space-1) var(--msui-inner);">
            <TextInput
              search
              placeholder="Filter, then Enter"
              aria-label="Filter"
              value={query}
              onInput={(v) => query.set(v)}
              onSubmit={(v) => submitted.set(v)}
            />
          </div>
          <div class="msui-form">
            <span class="msui-label">Typed</span>
            <span class="msui-value">{query}</span>
            <span class="msui-label">Submitted</span>
            <span class="msui-value">{submitted}</span>
          </div>
        </div>
      );
    }),
};
