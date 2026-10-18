/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import axe from 'axe-core';

// Automated accessibility checks of a mounted view in a unit test (`ui/_spec/accessibility.md`, Automated checks):
// axe-core over a DOM subtree in jsdom, WCAG 2.2 A and AA plus its best practices. Its own export
// (`@molstar/prism-ui/testing/axe`), so that stories, which import `testing`, do not load it. Only tests import it.

/**
 * Rules jsdom cannot evaluate, always off. jsdom has no layout and does not load the CSS modules, so a rule that
 * reads painted colors or boxes would only guess. Their gaps are in the audit table of `accessibility.md`.
 */
export const AxeUnavailableInJsdom: readonly string[] = [
  // Needs painted colors; the tokens are checked by the contrast of the themes, not here.
  'color-contrast',
  'color-contrast-enhanced',
  // Needs the rendered size of the target; `--msp-touch-target` and the 24px floor are CSS (`styles:check`).
  'target-size',
  // Needs scroll boxes (layout).
  'scrollable-region-focusable',
  // Whole-page rules: a view is mounted in a bare container, not in a page with a title, language and landmarks.
  'document-title',
  'html-has-lang',
  'html-lang-valid',
  'landmark-one-main',
  'page-has-heading-one',
  'region',
  'bypass',
];

export interface AxeOptions {
  /** More rules to turn off for this view, each with the reason in a comment where it is passed. */
  readonly disable?: readonly string[];
}

export type AxeViolation = axe.Result;

/** The violations axe finds in `root` (an element, its shadow tree included). */
export async function axeViolations(root: Element, options: AxeOptions = {}): Promise<AxeViolation[]> {
  const rules: axe.RuleObject = {};
  for (const id of [...AxeUnavailableInJsdom, ...(options.disable ?? [])]) rules[id] = { enabled: false };
  const result = await axe.run(root, {
    runOnly: { type: 'tag', values: ['wcag2a', 'wcag2aa', 'wcag21a', 'wcag21aa', 'wcag22aa', 'best-practice'] },
    rules,
    // The default is the main document's frames; a view has none.
    iframes: false,
    resultTypes: ['violations'],
  });
  return result.violations;
}

/** A readable report: the rule, its help, and the elements with the failing part. */
export function describeViolations(violations: readonly AxeViolation[]): string {
  return violations
    .map((v) => {
      const nodes = v.nodes
        .slice(0, 6)
        .map(
          (n) =>
            `    ${n.target.join(' ')}\n      ${n.html.slice(0, 160)}\n      ${n.failureSummary?.replace(/\n/g, '\n      ')}`,
        )
        .join('\n');
      return `${v.id} (${v.impact}): ${v.help} [${v.nodes.length}]\n${nodes}`;
    })
    .join('\n\n');
}

/** Fails, with the report, if axe finds a violation in `root`. */
export async function expectNoAxeViolations(root: Element, options: AxeOptions = {}): Promise<void> {
  const violations = await axeViolations(root, options);
  if (violations.length)
    throw new Error(`${violations.length} accessibility violation(s):\n\n${describeViolations(violations)}`);
}
