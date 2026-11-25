/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Generated from markdown.module.css by scripts/css.mjs; do not edit.

import { defineStyles } from '../../../styles/registry.js';

export const classes = {
  'markdown': 'msui-markdown-markdown',
  'media': 'msui-markdown-media',
  'table': 'msui-markdown-table',
} as const;

defineStyles('components/markdown/markdown.module.css', `@layer msui.components {
.msui-markdown-markdown {
  --md-text: var(--msui-text-sm);
  --md-leading: var(--msui-leading-sm);
  --md-heading: var(--msui-text-md);
  --md-heading-leading: var(--msui-leading-md);
  --md-code: var(--msui-text-xs);
  --md-gap: var(--msui-space-2);
  gap: var(--md-gap);
  min-width: 0;
  font-size: var(--md-text);
  line-height: var(--md-leading);
  display: grid;

  &[data-size="md"] {
    --md-text: var(--msui-text-md);
    --md-leading: var(--msui-leading-md);
    --md-heading: var(--msui-text-lg);
    --md-heading-leading: var(--msui-leading-lg);
    --md-code: var(--msui-text-sm);
    --md-gap: var(--msui-space-3);
  }

  & h3, & h4, & h5, & h6 {
    font-size: var(--md-text);
    line-height: var(--md-leading);
    font-weight: var(--msui-weight-strong);
    margin: 0;
  }

  & h3, & h4 {
    font-size: var(--md-heading);
    line-height: var(--md-heading-leading);
  }

  & p {
    margin: 0;
  }

  & ul, & ol {
    gap: var(--msui-space-1);
    margin: 0;
    padding-inline-start: var(--msui-indent);
    display: grid;
  }

  & ul {
    list-style: outside;
  }

  & ol {
    list-style: decimal;
  }

  & code {
    padding: 0 var(--msui-space-1);
    border-radius: var(--msui-radius-sm);
    background: var(--msui-control);
    font-family: var(--msui-font-mono);
    font-size: var(--md-code);
    line-height: 1;
  }

  & pre {
    padding: var(--msui-space-2);
    border-radius: var(--msui-radius-md);
    background: var(--msui-control);
    font-family: var(--msui-font-mono);
    font-size: var(--md-code);
    line-height: var(--msui-leading-sm);
    margin: 0;
    overflow-x: auto;

    & code {
      line-height: inherit;
      background: none;
      padding: 0;
    }
  }

  & hr {
    border: 0;
    border-top: 1px solid var(--msui-line);
    margin: 0 0 -1px;
  }

  & img {
    border-radius: var(--msui-radius-sm);
    max-width: 100%;
  }

  & a {
    color: var(--msui-accent);
    text-decoration: none;

    &:hover {
      text-decoration: underline;
    }

    &:focus-visible {
      border-radius: var(--msui-radius-sm);
      outline: 2px solid var(--msui-focus);
    }
  }
}

.msui-markdown-media {
  display: grid;

  & audio {
    width: 100%;
    height: var(--msui-header);
  }
}

.msui-markdown-table {
  overflow-x: auto;

  & table {
    border-collapse: collapse;
    font-variant-numeric: tabular-nums;
  }

  & th, & td {
    padding: var(--msui-space-1) var(--msui-space-2);
    box-shadow: inset 0 -1px var(--msui-line);
    text-align: left;
  }

  & th:first-child, & td:first-child {
    padding-inline-start: 0;
  }

  & th:last-child, & td:last-child {
    padding-inline-end: 0;
  }

  & th {
    font-weight: var(--msui-weight-strong);
  }
}
}
`);
