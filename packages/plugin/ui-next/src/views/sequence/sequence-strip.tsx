/** @jsxImportSource @molstar/plugin-ui-next/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Button } from '../../components/button/button.js';
import { Icon } from '../../components/icon/icon.js';
import { classes } from './_generated/sequence-strip.module.css.js';

export interface SequenceChain {
  /** The structure it belongs to ("1TQN"). */
  structure: string;
  /** The chain ("A") and its entity description. */
  name: string;
  info?: string;
  /** One-letter codes. */
  sequence: string;
  /** The number of the first residue. */
  start: number;
}

/**
 * A sequence strip for the top dock (a mock of the sequence view). Residues are rendered once; hover and selection
 * are applied to the existing nodes, so thousands of residues cost no subscriptions or re-rendering.
 */
export function SequenceStrip(props: {
  chain: SequenceChain;
  onHover?: (seqNumber: number | undefined) => void;
  onToggle?: (seqNumber: number, selected: boolean) => void;
  /** Collapses the strip (to the dock's slim bar). */
  onCollapse?: () => void;
}) {
  const { chain } = props;
  const residues = (<div class={classes.residues} />) as HTMLElement;
  const fragment = document.createDocumentFragment();
  for (let i = 0; i < chain.sequence.length; i++) {
    const span = document.createElement('span');
    span.className = classes.residue;
    span.textContent = chain.sequence[i];
    const number = chain.start + i;
    span.dataset.seq = String(number);
    if (number % 10 === 0) span.dataset.number = String(number);
    fragment.appendChild(span);
  }
  residues.appendChild(fragment);
  // One listener per strip (delegation), not per residue.
  residues.addEventListener('pointerover', (e) => {
    const seq = (e.target as HTMLElement).dataset?.seq;
    props.onHover?.(seq === undefined ? undefined : Number(seq));
  });
  residues.addEventListener('pointerleave', () => props.onHover?.(undefined));
  residues.addEventListener('click', (e) => {
    const target = e.target as HTMLElement;
    if (target.dataset?.seq === undefined) return;
    const selected = target.classList.toggle(classes.selected);
    props.onToggle?.(Number(target.dataset.seq), selected);
  });
  // The current structure and chain as selectors on top (as in the classic UI), so the residues get the full width.
  return (
    <div class={classes.strip}>
      <div class={classes.header}>
        <button type="button" class={classes.selector} title="Structure">
          {chain.structure} <Icon name="chevron-down" size={12} />
        </button>
        <button type="button" class={classes.selector} title="Chain">
          {chain.name}
          {chain.info ? <span class={classes.info}>{chain.info}</span> : null}
          <Icon name="chevron-down" size={12} />
        </button>
        <span class={classes.spacer} />
        {props.onCollapse ? (
          <Button
            variant="ghost"
            size="sm"
            icon={<Icon name="chevron-up" />}
            title="Hide the sequence"
            onClick={() => props.onCollapse?.()}
          />
        ) : null}
      </div>
      {residues}
    </div>
  );
}
