/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Generated from reset.css by scripts/css.mjs; do not edit.

import { defineStyles } from '../registry.js';

defineStyles('styles/reset.css', `@property --msui-fade {
  syntax: "<length>";
  inherits: false;
  initial-value: 0;
}

@supports (animation-timeline: scroll()) {
  @keyframes msui-fade-overflow {
    from, to {
      --msui-fade: 20px;
    }
  }
}

@property --msui-fade-top {
  syntax: "<length>";
  inherits: false;
  initial-value: 0;
}

@property --msui-fade-bottom {
  syntax: "<length>";
  inherits: false;
  initial-value: 0;
}

@supports (animation-timeline: scroll()) {
  @keyframes msui-fade-scroll-top {
    from {
      --msui-fade-top: 0px;
    }

    to {
      --msui-fade-top: 16px;
    }
  }

  @keyframes msui-fade-scroll-bottom {
    from {
      --msui-fade-bottom: 16px;
    }

    to {
      --msui-fade-bottom: 0px;
    }
  }
}

@keyframes msui-reveal {
  from {
    opacity: 0;
    height: 0;
    margin-block: 0;
    padding-block: 0;
    overflow: clip;
  }

  to {
    overflow: clip;
  }
}

@keyframes msui-reveal-side {
  from {
    opacity: 0;
    translate: 12px;
  }
}

@keyframes msui-reveal-fade {
  from {
    opacity: 0;
  }
}

@layer msui.reset {
  :where(.msui-root) {
    interpolate-size: allow-keywords;
    font-family: var(--msui-font-sans);
    font-size: var(--msui-text-sm);
    line-height: var(--msui-leading-sm);
    color: var(--msui-text);
    font-style: normal;
    font-weight: var(--msui-weight-regular);
    font-variant: normal;
    text-align: start;
    text-transform: none;
    text-indent: 0;
    letter-spacing: normal;
    word-spacing: normal;
    white-space: normal;
    -webkit-font-smoothing: antialiased;
    text-size-adjust: none;
    cursor: default;
  }

  :where(.msui-root) :where(input:not([type="range"], [type="checkbox"], [type="radio"]), textarea, [contenteditable]) {
    cursor: text;
  }

  :where(.msui-root), :where(.msui-root) *, :where(.msui-root) :before, :where(.msui-root) :after {
    box-sizing: border-box;
  }

  :where(.msui-root) :where(h1, h2, h3, h4, h5, h6, p, ul, ol, dl, dd, figure, blockquote, pre) {
    margin: 0;
  }

  :where(.msui-root) :where(ul, ol) {
    padding: 0;
    list-style: none;
  }

  :where(.msui-root) :where(button, input, select, textarea) {
    font: inherit;
    color: inherit;
    letter-spacing: inherit;
  }

  :where(.msui-root) :where(button) {
    cursor: pointer;
    background: none;
    border: 0;
    margin: 0;
    padding: 0;
  }

  :where(.msui-root) :where(button:disabled) {
    cursor: default;
  }

  :where(.msui-root) :where(svg) {
    flex-shrink: 0;
    display: block;
  }

  :where(.msui-root) :where(:focus-visible) {
    outline: 2px solid var(--msui-focus);
    outline-offset: 1px;
  }

  :where(.msui-root) :where([data-reveal]) {
    animation: msui-reveal var(--msui-motion-base) var(--msui-ease-in-out);
  }

  :where(.msui-root) :where([data-reveal="side"]) {
    animation-name: msui-reveal-side;
  }

  :where(.msui-root) :where([data-reveal="fade"]) {
    animation: msui-reveal-fade var(--msui-motion-slow) var(--msui-ease-in-out);
  }

  @supports (corner-shape: squircle) {
    :where(.msui-root), :where(.msui-root) *, :where(.msui-root) :before, :where(.msui-root) :after {
      corner-shape: var(--msui-corner-shape);
    }
  }

  :where(.msui-root) :where([hidden]) {
    display: none !important;
  }
}
`);
