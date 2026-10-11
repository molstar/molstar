/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Generated from reset.css by scripts/css.mjs; do not edit.

import { defineStyles } from '../registry.js';

defineStyles('styles/reset.css', `@property --msp-fade {
  syntax: "<length>";
  inherits: false;
  initial-value: 0;
}

@keyframes msp-fade-overflow {
  from, to {
    --msp-fade: 20px;
  }
}

@layer msp.reset {
  :where(.msp-root) {
    font-family: var(--msp-font-sans);
    font-size: var(--msp-font-size-md);
    line-height: var(--msp-line-height);
    color: var(--msp-color-text);
    -webkit-font-smoothing: antialiased;
    text-size-adjust: none;
    cursor: default;
  }

  :where(.msp-root) :where(input:not([type="range"], [type="checkbox"], [type="radio"]), textarea, [contenteditable]) {
    cursor: text;
  }

  :where(.msp-root), :where(.msp-root) *, :where(.msp-root) :before, :where(.msp-root) :after {
    box-sizing: border-box;
  }

  :where(.msp-root) :where(h1, h2, h3, h4, h5, h6, p, ul, ol, dl, dd, figure, blockquote, pre) {
    margin: 0;
  }

  :where(.msp-root) :where(ul, ol) {
    padding: 0;
    list-style: none;
  }

  :where(.msp-root) :where(button, input, select, textarea) {
    font: inherit;
    color: inherit;
    letter-spacing: inherit;
  }

  :where(.msp-root) :where(button) {
    cursor: pointer;
    background: none;
    border: 0;
    margin: 0;
    padding: 0;
  }

  :where(.msp-root) :where(button:disabled) {
    cursor: default;
  }

  :where(.msp-root) :where(svg) {
    flex-shrink: 0;
    display: block;
  }

  :where(.msp-root) :where(:focus-visible) {
    outline: 2px solid var(--msp-color-focus-ring);
    outline-offset: 1px;
  }

  :where(.msp-root) :where([hidden]) {
    display: none !important;
  }

  @media (prefers-reduced-motion: reduce) {
    :where(.msp-root), :where(.msp-root) *, :where(.msp-root) :before, :where(.msp-root) :after {
      transition-duration: .01ms !important;
      animation-duration: .01ms !important;
      animation-iteration-count: 1 !important;
    }
  }
}
`);
