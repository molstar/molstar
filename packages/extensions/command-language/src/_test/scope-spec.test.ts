/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as fs from 'node:fs';
import { describe, expect, it } from 'vitest';
import { PymolReps, PymolRepsUnsupported } from '../language/pymol/args.js';
import { PymolCommands } from '../language/pymol/commands.js';
import { PymolSettings } from '../language/pymol/settings.js';
import { PymolUnsupported } from '../language/pymol/unsupported.js';
import type { RepType } from '../language/ir.js';
import type { Unsupported } from '../language/table.js';
import { VmdCommands } from '../language/vmd/commands.js';
import { VmdColorMethodNames, VmdStyleNames } from '../language/vmd/styles.js';
import { VmdUnsupported } from '../language/vmd/unsupported.js';
import { MolstarReps } from '../prism-viewer/reps.js';

// The scope tables (_spec/pymol.md, _spec/vmd.md) are written from the command tables, so they say what the code
// does. This test fails when a file is not what the code says; `UPDATE_SPEC=1 pnpm test:unit
// packages/extensions/command-language/src/_test/scope-spec.test.ts` writes them.

const cell = (s: string) => s.replace(/\|/g, '\\|').replace(/\n/g, ' ');
const code = (s: string) => `\`${cell(s)}\``;

const molstar = (reps: readonly RepType[] | 'all') =>
  reps === 'all'
    ? 'every type'
    : reps
        .map((r) => {
          const m = MolstarReps[r];
          const extra = [
            m.typeParams
              ? Object.entries(m.typeParams)
                  .map(([k, v]) => `${k}: ${JSON.stringify(v)}`)
                  .join(', ')
              : '',
            m.sizeTheme ? `size ${m.sizeTheme.name}` : '',
          ].filter(Boolean);
          return `\`${m.type}\`${extra.length ? ` (${extra.join('; ')})` : ''}`;
        })
        .join(', ');

const unsupportedTable = (list: readonly Unsupported[]) => [
  '| Command | Why | Instead |',
  '| --- | --- | --- |',
  ...list.map((u) => `| ${code(u.name)} | ${cell(u.why)} | ${cell(u.hint ?? '')} |`),
];

function pymol(): string {
  return [
    '# PyMOL: scope',
    '',
    'What the PyMOL dialect does ([README](README.md)), written from the command tables by',
    '`src/_test/scope-spec.test.ts` (`UPDATE_SPEC=1` writes it). An omitted selection is `all`, as in PyMOL.',
    '',
    '## Commands',
    '',
    '| Command | Usage | What it does in Mol\\* |',
    '| --- | --- | --- |',
    ...PymolCommands.map(
      (c) =>
        `| ${code(c.name)}${c.aliases?.length ? ` (${c.aliases.map(code).join(', ')})` : ''} | ${code(c.usage)} | ${cell(c.summary)} |`,
    ),
    '',
    'Also `help [command]`, and `@file.pml` or `load file.pml`: a script given with this one runs in its place.',
    '',
    '## Representations',
    '',
    '| PyMOL | Mol\\* |',
    '| --- | --- |',
    ...Object.entries(PymolReps).map(
      ([k, v]) =>
        `| ${code(k)} | ${molstar(v)}${k === 'nb_spheres' || k === 'nonbonded' ? ', atoms without bonds' : ''} |`,
    ),
    ...Object.entries(PymolRepsUnsupported).map(([k, why]) => `| ${code(k)} | not supported: ${cell(why)} |`),
    '',
    '## Settings',
    '',
    '| Setting | What | Default |',
    '| --- | --- | --- |',
    ...Object.entries(PymolSettings).map(([k, s]) => `| ${code(k)} | ${cell(s.summary)} | ${code(s.default)} |`),
    '',
    'Any other setting fails, naming these.',
    '',
    '## Not supported',
    '',
    ...unsupportedTable(PymolUnsupported),
    '',
  ].join('\n');
}

function vmd(): string {
  return [
    '# VMD: scope',
    '',
    'What the VMD dialect does ([README](README.md)), written from the command tables by',
    '`src/_test/scope-spec.test.ts` (`UPDATE_SPEC=1` writes it). `molid` is a molecule id or `top`; a selection applies',
    'to the command’s molecule. Tcl: `set`, `unset`, `$name`, `[atomselect …]`, `[molinfo top]`, `puts`; atomselect',
    'handles’ `num`, `text`, `molid`, `get`, `delete`.',
    '',
    '## Commands',
    '',
    '| Command | Usage | What it does in Mol\\* |',
    '| --- | --- | --- |',
    ...VmdCommands.map(
      (c) =>
        `| ${code(c.name)}${c.aliases?.length ? ` (${c.aliases.map(code).join(', ')})` : ''} | ${code(c.usage)} | ${cell(c.summary)} |`,
    ),
    '',
    'Also `help [command]`, and `source file.tcl` or `play file.vmd`: a script given with this one runs in its place.',
    '',
    '## Styles and coloring methods',
    '',
    `Styles: ${VmdStyleNames.map(code).join(', ')}. Coloring methods: ${VmdColorMethodNames.map(code).join(', ')}.`,
    'Materials: `Opaque`, `Transparent` (others draw as Opaque, with a warning).',
    '',
    '## Not supported',
    '',
    ...unsupportedTable(VmdUnsupported),
    '',
  ].join('\n');
}

describe('the scope spec', () => {
  for (const [name, text] of [
    ['pymol.md', pymol()],
    ['vmd.md', vmd()],
  ] as const) {
    it(`${name} is what the code says`, () => {
      const file = new URL(`../../_spec/${name}`, import.meta.url);
      if (process.env.UPDATE_SPEC) fs.writeFileSync(file, text);
      // a Windows checkout has the file with CRLF line endings
      expect(fs.existsSync(file) ? fs.readFileSync(file, 'utf8').replace(/\r\n/g, '\n') : '').toBe(text);
    });
  }
});
