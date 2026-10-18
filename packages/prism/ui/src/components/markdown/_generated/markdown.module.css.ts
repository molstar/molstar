/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Generated from markdown.module.css by scripts/css.mjs; do not edit.

import { defineStyles } from '../../../styles/registry.js';

export const classes = {
  'markdown': 'msp-markdown-markdown',
  'media': 'msp-markdown-media',
  'table': 'msp-markdown-table',
} as const;

defineStyles('components/markdown/markdown.module.css', `@layer msp.components {
.msp-markdown-markdown {
  gap: var(--msp-space-3);
  line-height: 1.5;
  display: grid;

  & h3, & h4, & h5, & h6 {
    font-size: var(--msp-font-size-md);
    font-weight: 600;
  }

  & ul, & ol {
    gap: var(--msp-space-1);
    padding-left: var(--msp-space-5);
    display: grid;
  }

  & ul {
    list-style: outside;
  }

  & ol {
    list-style: decimal;
  }

  & code {
    padding: 0 var(--msp-space-1);
    border-radius: var(--msp-radius-sm);
    background: var(--msp-color-control);
    font-family: var(--msp-font-mono);
    font-size: .92em;
  }

  & pre {
    padding: var(--msp-space-2);
    border-radius: var(--msp-radius-sm);
    background: var(--msp-color-control);
    overflow-x: auto;

    & code {
      background: none;
      padding: 0;
    }
  }

  & hr {
    border: 0;
    border-top: 1px solid var(--msp-color-border);
  }

  & img {
    border-radius: var(--msp-radius-sm);
    max-width: 100%;
  }

  & a {
    color: var(--msp-color-accent);
    text-decoration: none;

    &:hover {
      text-decoration: underline;
    }
  }
}

.msp-markdown-media {
  display: grid;

  & audio {
    width: 100%;
    height: 32px;
  }
}

.msp-markdown-table {
  overflow-x: auto;

  & table {
    border-collapse: collapse;
    font-variant-numeric: tabular-nums;
  }

  & th, & td {
    padding: var(--msp-space-1) var(--msp-space-2);
    border-bottom: 1px solid var(--msp-color-border);
    text-align: left;
  }

  & th {
    font-weight: 600;
  }
}
}
`);
