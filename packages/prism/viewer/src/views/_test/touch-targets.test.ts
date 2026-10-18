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

// Touch targets of the viewer's own controls (accessibility.md, Targets; the shell's: prism-ui's `touch-targets.test.ts`):
// the scene tree's rows and buttons and the canvas's selection bar grow to `--msp-touch-target` in compact mode. jsdom
// has no layout, so the guard reads the rules.

const views = resolve(import.meta.dirname, '..');

const rule = (css: string, name: string) => new RegExp(`(?:^|\\n)\\.${name}\\s*\\{([^}]*)\\}`).exec(css)?.[1];

const Targets: [file: string, classes: string[]][] = [
  ['structure/scene-tree.module.css', ['row', 'expander', 'action', 'menuItem']],
  ['canvas/selection-bar.module.css', ['bar', 'action']],
  ['state/snapshots-panel.module.css', ['item', 'action']],
];

describe('touch targets', () => {
  for (const [file, classes] of Targets) {
    it(`${file}: ${classes.join(', ')} grow to the touch target`, () => {
      const css = readFileSync(resolve(views, file), 'utf8');
      for (const name of classes) {
        const body = rule(css, name);
        expect(body, `.${name} is in ${file}`).toBeDefined();
        expect(body, `.${name} grows to the touch target`).toContain('--msp-touch-target');
      }
    });
  }
});
