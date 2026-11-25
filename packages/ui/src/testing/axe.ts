/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import axe from 'axe-core';

// axe for the views' tests: WCAG 2.2 A and AA and best practices, less what jsdom cannot evaluate (no layout, no painted
// colors) and the page-level rules (a view is mounted in a bare container).

const Tags = ['wcag2a', 'wcag2aa', 'wcag21a', 'wcag21aa', 'wcag22aa', 'best-practice'];
const Off = ['color-contrast', 'region', 'target-size'];

/** The violations under `root`, one per line with its nodes; '' when there are none. */
export async function axeViolations(root: Element, off: readonly string[] = []): Promise<string> {
  const rules: Record<string, { enabled: boolean }> = {};
  for (const id of [...Off, ...off]) rules[id] = { enabled: false };
  const r = await axe.run(root, { runOnly: { type: 'tag', values: Tags }, rules, resultTypes: ['violations'] });
  return r.violations.map((v) => `${v.id}: ${v.help}\n${v.nodes.map((n) => `  ${n.html}`).join('\n')}`).join('\n');
}
