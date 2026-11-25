/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { of } from '@molstar/core/reactive';
import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { type Child } from '../../../dom/index.js';
import { Heading, Note, page, story } from '../../../_stories/support.js';
import { state } from '../../../state.js';
import { Icon } from '../../icon/icon.js';
import { ActionFooter } from '../action-footer.js';

const meta: Meta = { title: 'Components/Action footer' };
export default meta;

// A block on a surface (styles §5): a row with an icon and a form row above the footer, so its reset sits in the icon
// slot and its info at the text column; the action's label ends at the inset.
function Block(props: { title: string; children: Child }) {
  return (
    <div style="display: grid; gap: var(--msui-space-1);">
      <Heading>{props.title}</Heading>
      <div class="msui-surface" style="width: 320px;">
        <div class="msui-row">
          <span class="msui-slot">
            <Icon name="camera" />
          </span>
          <span>Image</span>
        </div>
        <div class="msui-form">
          <label class="msui-label" htmlFor="story-af-name">
            Name
          </label>
          <input id="story-af-name" class="msui-control" value="1cbs" />
        </div>
        {props.children}
      </div>
    </div>
  );
}

/** A plain footer, one with reset and an info line, and one with a secondary action; each echoes what it did. */
export const Layouts: StoryObj = {
  render: () =>
    story(() => {
      const last = state('none');
      const did = (what: string) => () => last.set(what);
      return (
        <div style={page}>
          <Block title="Plain">
            <ActionFooter label="Download" onApply={did('download')} />
          </Block>
          <Block title="Reset and info">
            <ActionFooter label="Save" info="1920 × 1080 px" onReset={did('reset')} onApply={did('save')} />
          </Block>
          <Block title="Secondary action">
            <ActionFooter
              label="Download"
              secondary={{ label: 'Copy', onClick: did('copy') }}
              onApply={did('download')}
            />
          </Block>
          <Note>last: {last.asReadonly} (Enter in the field applies too)</Note>
        </div>
      );
    }),
};

/** Busy: applying shows a spinner and disables the action for 1.5 s; disabled: not possible yet, no spinner. */
export const States: StoryObj = {
  render: () =>
    story(() => {
      const busy = state(false);
      const apply = () => {
        busy.set(true);
        setTimeout(() => busy.set(false), 1500);
      };
      return (
        <div style={page}>
          <Block title="Busy (press Apply)">
            <ActionFooter
              label="Apply"
              busy={busy.asReadonly}
              secondary={{ label: 'Copy', onClick: () => {} }}
              onReset={() => {}}
              onApply={apply}
            />
          </Block>
          <Block title="Disabled">
            <ActionFooter label="Apply" disabled={of(true)} onApply={() => {}} />
          </Block>
        </div>
      );
    }),
};
