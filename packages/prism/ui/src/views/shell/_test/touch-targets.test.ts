/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { readFileSync } from 'node:fs';
import { resolve } from 'node:path';
import { describe, expect, it } from 'vitest';

// Touch targets (accessibility.md, Targets): in compact mode the shell sets `--msp-touch-target` to 44px, and every
// fixed-size control of the rail, the tab strip and the canvas chrome grows to it (`max(<size>, var(--msp-touch-target,
// 0px))`). jsdom has no layout, so the guard reads the rules: a control that loses the clause fails here.

const src = resolve(import.meta.dirname, '../../..');

const rule = (css: string, name: string) => new RegExp(`(?:^|\\n)\\.${name}\\s*\\{([^}]*)\\}`).exec(css)?.[1];

const Targets: [file: string, classes: string[]][] = [
  ['views/shell/dock.module.css', ['tab']],
  ['views/shell/canvas.module.css', ['tool', 'pill', 'pillButton', 'railMenuButton', 'stackToggle']],
  ['components/button/button.module.css', ['small', 'tiny']],
];

describe('touch targets', () => {
  it('the compact shell sets the 44px target', () => {
    const shell = readFileSync(resolve(src, 'views/shell/shell.module.css'), 'utf8');
    expect(shell).toMatch(/\[data-mode="compact"\][^}]*--msp-touch-target:\s*44px/);
  });

  for (const [file, classes] of Targets) {
    it(`${file}: ${classes.join(', ')} grow to it`, () => {
      const css = readFileSync(resolve(src, file), 'utf8');
      for (const name of classes) {
        const body = rule(css, name);
        expect(body, `.${name} is in ${file}`).toBeDefined();
        expect(body, `.${name} grows to the touch target`).toContain('--msp-touch-target');
      }
    });
  }
});
