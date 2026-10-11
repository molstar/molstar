/** @jsxImportSource @molstar/plugin-ui-next/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { classes } from './_generated/sequence-strip.module.css.js';

export interface SequenceChain {
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
  return (
    <div class={classes.strip}>
      <div class={classes.chain}>
        <span class={classes.chainName}>{chain.name}</span>
        {chain.info ? <span class={classes.chainInfo}>{chain.info}</span> : null}
      </div>
      {residues}
    </div>
  );
}
