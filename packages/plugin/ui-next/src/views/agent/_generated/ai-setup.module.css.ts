/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Generated from ai-setup.module.css by scripts/css.mjs; do not edit.

import { defineStyles } from '../../../styles/registry.js';

export const classes = {
  'action': 'msp-agent-ai-setup-action',
  'caption': 'msp-agent-ai-setup-caption',
  'configured': 'msp-agent-ai-setup-configured',
  'copy': 'msp-agent-ai-setup-copy',
  'fold': 'msp-agent-ai-setup-fold',
  'grow': 'msp-agent-ai-setup-grow',
  'keyRow': 'msp-agent-ai-setup-keyRow',
  'label': 'msp-agent-ai-setup-label',
  'lead': 'msp-agent-ai-setup-lead',
  'link': 'msp-agent-ai-setup-link',
  'needed': 'msp-agent-ai-setup-needed',
  'note': 'msp-agent-ai-setup-note',
  'own': 'msp-agent-ai-setup-own',
  'ownFoot': 'msp-agent-ai-setup-ownFoot',
  'provider': 'msp-agent-ai-setup-provider',
  'providerHead': 'msp-agent-ai-setup-providerHead',
  'providerName': 'msp-agent-ai-setup-providerName',
  'providerState': 'msp-agent-ai-setup-providerState',
  'providers': 'msp-agent-ai-setup-providers',
  'row': 'msp-agent-ai-setup-row',
  'setup': 'msp-agent-ai-setup-setup',
} as const;

defineStyles('views/agent/ai-setup.module.css', `@layer msp.components {
.msp-agent-ai-setup-setup {
  font-size: var(--msp-font-size-sm);
}

.msp-agent-ai-setup-needed, .msp-agent-ai-setup-configured, .msp-agent-ai-setup-own {
  gap: var(--msp-space-2);
  display: grid;
}

.msp-agent-ai-setup-lead {
  padding: 0 var(--msp-space-3);
  margin: 0;
  font-weight: 600;
}

.msp-agent-ai-setup-copy {
  padding: 0 var(--msp-space-3);
  color: var(--msp-color-text-muted);
  margin: 0;
  line-height: 1.45;
}

.msp-agent-ai-setup-note {
  padding: 0 var(--msp-space-3);
  color: var(--msp-color-text-subtle);
  font-size: var(--msp-font-size-xs);
  margin: 0;
}

.msp-agent-ai-setup-caption {
  margin: var(--msp-space-1) 0 0;
  padding: 0 var(--msp-space-3);
  color: var(--msp-color-text-subtle);
  font-size: var(--msp-font-size-xs);
  font-weight: 500;
}

.msp-agent-ai-setup-row {
  grid-template-columns: var(--msp-label-width, 96px) minmax(0, 1fr);
  align-items: start;
  gap: var(--msp-space-3);
  padding-left: var(--msp-space-3);
  display: grid;
}

.msp-agent-ai-setup-label {
  color: var(--msp-color-text-muted);
  line-height: 24px;
}

.msp-agent-ai-setup-grow {
  flex: 1;
  min-width: 0;
}

.msp-agent-ai-setup-keyRow {
  gap: var(--msp-space-2);
  display: flex;
}

.msp-agent-ai-setup-action {
  height: 24px;
  padding: 0 var(--msp-space-3);
  border-radius: var(--msp-radius-md);
  background: color-mix(in oklab, var(--msp-color-accent), transparent 84%);
  color: var(--msp-color-accent-strong);
  flex: none;
  font-weight: 600;

  &:hover {
    background: color-mix(in oklab, var(--msp-color-accent), transparent 76%);
  }
}

.msp-agent-ai-setup-providers {
  display: grid;
}

.msp-agent-ai-setup-provider {
  gap: var(--msp-space-1);
  padding: 2px 0;
  display: grid;
}

.msp-agent-ai-setup-providerHead {
  align-items: center;
  gap: var(--msp-space-2);
  min-height: 24px;
  padding-left: var(--msp-space-3);
  display: flex;
}

.msp-agent-ai-setup-providerName {
  min-width: 0;
}

.msp-agent-ai-setup-providerState {
  color: var(--msp-color-text-subtle);
  font-size: var(--msp-font-size-xs);
  flex: 1;

  &[data-has] {
    color: var(--msp-color-text-muted);
  }
}

.msp-agent-ai-setup-link {
  height: 22px;
  padding: 0 var(--msp-space-2);
  border-radius: var(--msp-radius-md);
  color: var(--msp-color-text-muted);
  font-size: var(--msp-font-size-xs);

  &:hover {
    background: var(--msp-color-hover);
    color: var(--msp-color-text);
  }
}

.msp-agent-ai-setup-fold {
  height: 24px;
  padding: 0 var(--msp-space-3);
  border-radius: var(--msp-radius-md);
  color: var(--msp-color-text-muted);
  justify-self: start;

  &:hover {
    background: var(--msp-color-hover);
    color: var(--msp-color-text);
  }

  &:before {
    content: "+ ";
  }

  &[aria-expanded="true"]:before {
    content: "− ";
  }
}

.msp-agent-ai-setup-own {
  padding: var(--msp-space-2) 0;
  border-radius: var(--msp-radius-md);
  background: var(--msp-color-well);
}

.msp-agent-ai-setup-ownFoot {
  padding: 0 var(--msp-space-2);
  justify-content: flex-end;
  display: flex;
}
}
`);
