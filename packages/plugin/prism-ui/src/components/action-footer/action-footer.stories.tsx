/** @jsxImportSource @molstar/plugin-prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { of } from 'rxjs';
import { story } from '../../_stories/story-root.js';
import { observableState } from '../../reactive/observable-state.js';
import { ActionFooter } from './action-footer.js';

const meta: Meta = { title: 'Components/Action footer' };
export default meta;

const panel =
  'display: grid; gap: var(--msp-space-4); width: 300px; padding: var(--msp-space-4); background: var(--msp-color-panel); color: var(--msp-color-text); font-size: var(--msp-font-size-sm);';
// The footer undoes its container's 8px padding at the sides and bottom: the raised block an open action sits in.
const block = 'padding: var(--msp-space-3); border-radius: var(--msp-radius-md); background: var(--msp-color-raised);';
const caption = 'margin: 0 0 var(--msp-space-2); color: var(--msp-color-text-muted);';

function Block(props: { title: string; children: Node }) {
  return (
    <div>
      <p style={caption}>{props.title}</p>
      <div style={block}>{props.children}</div>
    </div>
  );
}

/** A plain footer, one with reset and an info line, and one with a secondary action; each echoes what it did. */
export const Layouts: StoryObj = {
  render: () =>
    story(() => {
      const last = observableState('none');
      const did = (what: string) => () => last.set(what);
      return (
        <div style={panel}>
          <Block title="Plain">
            <ActionFooter icon="download" label="Download" onApply={did('download')} />
          </Block>
          <Block title="Reset and info">
            <ActionFooter
              icon="camera"
              label="Save"
              info="1920 × 1080 px"
              onReset={did('reset')}
              onApply={did('save')}
            />
          </Block>
          <Block title="Secondary action">
            <ActionFooter
              icon="download"
              label="Download"
              secondary={{ icon: 'copy', label: 'Copy', onClick: did('copy') }}
              onApply={did('download')}
            />
          </Block>
          <span style="color: var(--msp-color-text-muted)">last: {last.asReadonly}</span>
        </div>
      );
    }),
};

/** Busy: applying shows a spinner and disables the action for 1.5 s; disabled: not possible yet, no spinner. */
export const States: StoryObj = {
  render: () =>
    story(() => {
      const busy = observableState(false);
      const apply = () => {
        busy.set(true);
        setTimeout(() => busy.set(false), 1500);
      };
      return (
        <div style={panel}>
          <Block title="Busy (press Apply)">
            <ActionFooter
              icon="download"
              label="Apply"
              busy={busy.asReadonly}
              secondary={{ icon: 'copy', label: 'Copy', onClick: () => {} }}
              onReset={() => {}}
              onApply={apply}
            />
          </Block>
          <Block title="Disabled">
            <ActionFooter icon="download" label="Apply" disabled={of(true)} onApply={() => {}} />
          </Block>
        </div>
      );
    }),
};
