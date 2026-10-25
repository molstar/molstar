/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as fs from 'node:fs';
import * as path from 'node:path';
import { describe, expect, it } from 'vitest';
import type { OperationDefinition } from '../api/area.js';
import { DefaultPrismViewerSpec } from '../default-spec.js';
import { isLazyExtension } from '../spec.js';
import { headlessViewer } from './harness.js';

// `_spec/operations.md` (the package's spec, README "Operations"): every operation of the default spec, written from the
// code. This test fails when the file is not what the code says; `UPDATE_SPEC=1 pnpm test:unit` writes it.

const file = path.resolve(import.meta.dirname, '../../_spec/operations.md');

function flags(d: OperationDefinition<any, any>): string {
  const out: string[] = [];
  if (d.read) out.push('read');
  if (d.history === false && !d.read) out.push('view');
  if (d.undoable === false) out.push('no undo');
  if (d.approval) out.push(d.approval + (d.approval === 'network' ? (d.addresses ? ' (addresses)' : ' (asks)') : ''));
  if (d.modules?.length) out.push(`needs ${d.modules.join(', ')}`);
  if (d.surface) out.push(d.surface.kind);
  return out.join('; ');
}

async function generate(): Promise<string> {
  const v = await headlessViewer({ spec: { options: { ui: { showRemoteStates: true } } } });
  for (const e of DefaultPrismViewerSpec().extensions ?? []) if (isLazyExtension(e)) await v.extensions.load(e.id);
  const rows = [...v.api.definitions.value.values()]
    .sort((a, b) => a.id.localeCompare(b.id))
    .map((d) => {
      const params = Object.keys(v.api.describe(d.id)?.params ?? {}).join(', ');
      return `| \`${d.id}\` | ${d.title} | ${flags(d)} | ${params} |`;
    });
  return [
    '# Operations',
    '',
    'Every operation of `DefaultPrismViewerSpec()` (the core and the built-in features, lazy ones loaded), written from',
    'the code by `src/_test/operations-spec.test.ts` (`UPDATE_SPEC=1 pnpm test:unit packages/prism/viewer/src/_test/operations-spec.test.ts`).',
    "Flags: `read` (no step), `view` (`history: false`), `no undo` (`undoable: false`), the agent's approval (`network`",
    "with `addresses` or always asking), the modules it needs, a `surface` limit. Mol\\*'s extensions list theirs in",
    'their own packages (`@molstar/<name>-extension/prism-viewer`).',
    '',
    '| Operation | Title | Flags | Params |',
    '| --- | --- | --- | --- |',
    ...rows,
    '',
  ].join('\n');
}

describe('the operations spec', () => {
  it('is what the code says', async () => {
    const text = await generate();
    if (process.env.UPDATE_SPEC) fs.writeFileSync(file, text);
    expect(fs.existsSync(file) ? fs.readFileSync(file, 'utf8').replace(/\r\n/g, '\n') : '').toBe(text);
  });
});
