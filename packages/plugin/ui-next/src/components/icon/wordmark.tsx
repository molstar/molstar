/** @jsxImportSource @molstar/plugin-ui-next/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { classes } from './_generated/wordmark.module.css.js';

// The Mol* wordmark as outlines (after the logo): the M, the rounded o, and the l whose foot runs into the asterisk's
// left arm (the asterisk is the Mol* mark's polygon, see molecular.ts). Logo units; strokes keep their width at any size.
const M = 'M85 90 L155 90 L237 262 L320 90 L390 90 L390 432 L318 432 L318 232 L237 398 L157 232 L157 432 L85 432Z';
const LStar =
  'M919.6 105.1 L980.4 105.1 L980.4 171.7 L1027.5 124.6 L1070.4 167.5 L1023.3 214.6 L1089.9 214.6 L1089.9 275.4 L1023.3 275.4 L1070.4 322.5 L1027.5 365.4 L980.4 318.3 L980.4 384.9 L919.6 384.9 L919.6 318.3 L872.5 365.4 L829.6 322.5 L876.7 275.4 L722 275.4 L722 432 L652 432 L652 70 L722 70 L722 214.6 L876.7 214.6 L829.6 167.5 L872.5 124.6 L919.6 171.7Z';

/** The Mol* wordmark, outlined in the text color: Home's heading. */
export function MolstarWordmark(props: { class?: string }) {
  return (
    <svg
      class={`${classes.wordmark} ${props.class ?? ''}`}
      viewBox="75 60 1025 382"
      role="img"
      aria-label="Mol*"
      fill="none"
      stroke="currentColor"
      stroke-linejoin="miter"
    >
      <path d={M} vector-effect="non-scaling-stroke" />
      <rect x="410" y="168" width="222" height="264" rx="100" vector-effect="non-scaling-stroke" />
      <rect x="470" y="228" width="102" height="144" rx="48" vector-effect="non-scaling-stroke" />
      <path d={LStar} vector-effect="non-scaling-stroke" />
    </svg>
  );
}
