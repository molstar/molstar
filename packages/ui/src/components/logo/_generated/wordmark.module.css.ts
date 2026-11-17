/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Generated from wordmark.module.css by scripts/css.mjs; do not edit.

import { defineStyles } from '../../../styles/registry.js';

export const classes = {
  'flash': 'msui-logo-wordmark-flash',
  'flash-hover': 'msui-logo-wordmark-flash-hover',
  'flash-often': 'msui-logo-wordmark-flash-often',
  'flash-rare': 'msui-logo-wordmark-flash-rare',
  'flash-reveal': 'msui-logo-wordmark-flash-reveal',
  'live': 'msui-logo-wordmark-live',
  'outline': 'msui-logo-wordmark-outline',
  'wordmark': 'msui-logo-wordmark-wordmark',
} as const;

defineStyles('components/logo/wordmark.module.css', `@layer msui.components {
.msui-logo-wordmark-wordmark {
  width: calc(57px * var(--msui-wordmark-scale, 1));
  height: calc(22px * var(--msui-wordmark-scale, 1));
  flex: none;
  display: block;
  position: relative;

  & svg {
    stroke-width: calc(1.5 / sqrt(var(--msui-wordmark-scale, 1)));
    position: absolute;
    inset: 0;
    overflow: visible;
  }
}

.msui-logo-wordmark-outline {
  color: inherit;
}

.msui-logo-wordmark-flash {
  animation: 14s ease-in-out .3s infinite both msui-logo-wordmark-flash-rare;
  mask-image: linear-gradient(100deg, #0000 0 38%, #000 47% 53%, #0000 62% 100%);
  mask-position: 100% 0;
  mask-size: 300% 100%;
  mask-repeat: no-repeat;

  [data-flash="often"] > & {
    animation: 7s ease-in-out .3s infinite both msui-logo-wordmark-flash-often;
  }

  [data-flash="reveal"] > & {
    animation: 1.4s ease-out both msui-logo-wordmark-flash-reveal;
    mask-image: linear-gradient(100deg, #000 0 42%, #0000 58% 100%);
  }

  .msui-logo-wordmark-wordmark:not([data-flash="reveal"], [data-live]):hover & {
    animation: 2.8s ease-in-out both msui-logo-wordmark-flash-hover;
  }
}

@keyframes msui-logo-wordmark-flash-rare {
  0% {
    mask-position: 100% 0;
  }

  20%, 100% {
    mask-position: 0 0;
  }
}

@keyframes msui-logo-wordmark-flash-often {
  0% {
    mask-position: 100% 0;
  }

  40%, 100% {
    mask-position: 0 0;
  }
}

@keyframes msui-logo-wordmark-flash-reveal {
  from {
    mask-position: 100% 0;
  }

  to {
    mask-position: 0 0;
  }
}

@keyframes msui-logo-wordmark-flash-hover {
  from {
    mask-position: 100% 0;
  }

  to {
    mask-position: 0 0;
  }
}

.msui-logo-wordmark-wordmark[data-live] {
  pointer-events: auto;
}

.msui-logo-wordmark-live {
  opacity: 0;
  --msui-wordmark-strip: calc(14px * var(--msui-wordmark-scale, 1));
  mask-image: linear-gradient(90deg,
    transparent calc(var(--msui-wordmark-x, 50%) - var(--msui-wordmark-strip)),
    #000 calc(var(--msui-wordmark-x, 50%) - var(--msui-wordmark-strip) / 3),
    #000 calc(var(--msui-wordmark-x, 50%) + var(--msui-wordmark-strip) / 3),
    transparent calc(var(--msui-wordmark-x, 50%) + var(--msui-wordmark-strip)));
  transition: opacity .6s var(--msui-ease);

  .msui-logo-wordmark-wordmark:hover > & {
    opacity: 1;
    transition-duration: .2s;
  }
}

@media (prefers-reduced-motion: reduce) {
  .msui-logo-wordmark-flash {
    animation: none;
  }
}
}
`);
