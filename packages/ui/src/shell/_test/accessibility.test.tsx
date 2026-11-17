/** @vitest-environment jsdom */
/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import axe from 'axe-core';
import { afterEach, describe, expect, it } from 'vitest';
import { flush } from '../../dom/index.js';
import type { ShellModel } from '../shell-model.js';
import { microtasks, mountShell, probeView } from './fixture.js';

// axe over the shell's states (shell §4): WCAG 2.2 A and AA and best practices, less what jsdom cannot evaluate (no
// layout, no painted colors) and the page-level rules (a shell is mounted in a bare container).
const Off = [
  'color-contrast',
  'color-contrast-enhanced',
  'target-size',
  'scrollable-region-focusable',
  'document-title',
  'html-has-lang',
  'html-lang-valid',
  'landmark-one-main',
  'page-has-heading-one',
  'region',
  'bypass',
];

async function violations(root: Element): Promise<string> {
  const rules: axe.RuleObject = {};
  for (const id of Off) rules[id] = { enabled: false };
  const r = await axe.run(root, {
    runOnly: { type: 'tag', values: ['wcag2a', 'wcag2aa', 'wcag21a', 'wcag21aa', 'wcag22aa', 'best-practice'] },
    rules,
    iframes: false,
    resultTypes: ['violations'],
  });
  return r.violations
    .map((v) => `${v.id}: ${v.help}\n${v.nodes.map((n) => `  ${n.html.slice(0, 140)}`).join('\n')}`)
    .join('\n');
}

let current: { dispose(): void } | undefined;
afterEach(() => {
  current?.dispose();
  current = undefined;
});

function shellWith(setup: (s: ShellModel) => void, options = {}) {
  const views = [
    probeView('home', { place: 'head', order: 0, compactBar: true }),
    probeView('layers', { order: 1, compactBar: true }),
    probeView('inspect', { order: 2 }),
    probeView('settings', { place: 'foot', order: 10 }),
    probeView('export', { home: 'card' }),
    probeView('report', { home: 'page' }),
    probeView('confirm', { home: 'modal' }),
    probeView('plot', { home: 'float', instances: 'many' }),
  ];
  const m = mountShell(views, options);
  m.shell.registerRailItems([{ id: 'theme', icon: 'sun-moon', title: 'Theme', run: () => {} }]);
  setup(m.shell);
  flush();
  current = m;
  return m;
}

describe('the shell is accessible (axe)', () => {
  const states: [string, (s: ShellModel) => void, object?][] = [
    ['folded to the rail', () => {}],
    ['a view open', (s) => s.show('home')],
    [
      'the column split, tabs below',
      (s) => {
        s.show('layers');
        s.splitBelow('inspect');
      },
    ],
    [
      'windows, one of two views',
      (s) => {
        s.show('plot', { key: 'a' });
        s.show('plot', { key: 'b' });
        s.float('layers');
        s.move('settings', { kind: 'group', id: s.layout.value.windows.at(-1)!.group.id });
      },
    ],
    ['a card', (s) => s.show('export')],
    ['a page', (s) => s.show('report')],
    ['a modal', (s) => s.show('confirm')],
    ['the rail hidden', (s) => s.toggleRail(false)],
    ['compact', (s) => s.show('layers'), { mode: 'compact' }],
    ['embedded', () => {}, { mode: 'embedded' }],
  ];
  for (const [name, setup, options] of states) {
    it(name, async () => {
      const m = shellWith(setup, options);
      await microtasks();
      flush();
      expect(await violations(m.root)).toBe('');
    });
  }
});
