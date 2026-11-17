/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { State } from '@molstar/core/reactive';
import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { story } from '../../_stories/support.js';
import { HoverLabel } from '../hover-label.js';
import { structureLines, StagePlaceholder } from '../placeholders.js';
import type { HoverLine } from '../stage-model.js';
import { StageBox } from './demo.js';

const meta: Meta = { title: 'Stage/Hover label' };
export default meta;

/** A label's lines: the bead's, and on every fifth bead a long one (it is cut at its cap and fades). */
function lines(i: number): readonly HoverLine[] {
  if (i % 5) return structureLines(i);
  return [
    ...structureLines(i),
    { role: 'item', runs: [{ text: '• Ligand ' }, { text: 'NAG', code: true }, { text: ' bound at 2.8 Å' }] },
    {
      role: 'text',
      runs: [
        { text: 'C', bold: true },
        { text: '6', sub: true },
        { text: 'H' },
        { text: '12', sub: true },
        { text: 'O' },
        { text: '6', sub: true },
      ],
    },
    {
      role: 'text',
      runs: [
        { text: 'See the entry', link: true },
        { text: ' (a link is shown, not followed)', small: true },
      ],
    },
  ];
}

function Over(props: { placement: 'corner' | 'cursor' }) {
  const hover = new State<readonly HoverLine[]>([]);
  let area!: HTMLElement;
  return (
    <StageBox>
      <div style="position: absolute; inset: 0; background: #f4f5f7" ref={(el: HTMLElement) => (area = el)}>
        <StagePlaceholder onHover={(i) => hover.set(i === undefined ? [] : lines(i))} />
        <HoverLabel model={{ lines: hover }} area={() => area} placement={props.placement} />
      </div>
    </StageBox>
  );
}

/** In the corner (the default): it stays put; hover the beads. */
export const Corner: StoryObj = { render: () => story(() => <Over placement="corner" />) };

/** Beside the pointer: placed once a frame, flipped near the edges; hidden while a button is held. */
export const Cursor: StoryObj = { render: () => story(() => <Over placement="cursor" />) };
