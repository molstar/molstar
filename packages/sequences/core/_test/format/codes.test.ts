/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

// V6's drift guards: every `svs.*` literal in the sources is registered, and every `model.*` literal of the model lands
// on a registered code through `svsIssue`. The codes' table is the generated reference (`_spec/reference/format.md`).

import { readdirSync, readFileSync } from 'node:fs';
import { join } from 'node:path';
import { fileURLToPath } from 'node:url';
import { Issue } from '@molstar/seq-core/model';
import { describe, expect, it } from 'vitest';
import { issueOf } from '../../src/format/index.js';
import { svsIssue } from '../../src/format/load/issues.js';
import { kernel } from './fixtures.js';

const here = (p: string) => fileURLToPath(new URL(p, import.meta.url));
const sources = (dir: string): [string, string][] =>
  readdirSync(dir, { recursive: true, encoding: 'utf8' })
    .filter((f) => f.endsWith('.ts'))
    .map((f) => [f, readFileSync(join(dir, f), 'utf8').replace(/\r\n/g, '\n')]);

describe('the issue codes', () => {
  const reg = kernel();

  it('every svs.* literal in the sources is registered', () => {
    const used = new Set<string>();
    for (const [, text] of sources(here('../../src/format'))) {
      for (const [, code] of text.matchAll(/["'](svs\.[a-z_]+(?:\.[a-z_]+)+)["']/g)) used.add(code);
    }
    expect(used.size).toBeGreaterThan(40);
    expect([...used].filter((c) => !reg.get('issueCodes', c))).toEqual([]);
  });

  it('svs.provider.unmapped is a registered warning', () => {
    expect(reg.get('issueCodes', 'svs.provider.unmapped')?.severity).toBe('warning');
  });

  it('issueOf takes the registered severity, an error when unregistered', () => {
    expect(issueOf(reg, 'svs.provider.unavailable', '/a', 'm')).toEqual(
      Issue.warning('svs.provider.unavailable', '/a', 'm'),
    );
    expect(issueOf(reg, 'svs.params.type', '', 'm').severity).toBe('error');
    expect(issueOf(reg, 'svs.nope.never', '', 'm').severity).toBe('error');
  });

  it('every model.* literal, pushed through svsIssue, lands on a registered code of the same severity', () => {
    const seen = new Map<string, Set<string>>();
    for (const [, raw] of sources(here('../../src/model'))) {
      // `const FORMAT = "model.source.format"` stands for its literal where it is used.
      let text = raw;
      for (const [, name, code] of raw.matchAll(/const ([A-Z_]+) = (["']model\.[a-z_.]+["']);/g)) {
        text = text.replaceAll(new RegExp(`\\b${name}\\b(?! =)`, 'g'), code);
      }
      for (const m of text.matchAll(/["'](model\.[a-z_]+(?:\.[a-z_]+)+)["']/g)) {
        const before = text.slice(Math.max(0, m.index - 80), m.index);
        if (/const [A-Z_]+ = $/.test(before)) continue; // the constant's own declaration
        const rest = text.slice(m.index);
        // A warning or a note where `Issue.warning(` or `Issue.info(` holds it, a warning for `report(…, true)`
        // (annotation.ts), or, in a table of codes (`["model.x", …]`), the next `Issue.` call's after it.
        const severity =
          /Issue\.(warning|info)\($/.exec(before)?.[1] ??
          (/^[^\n)]*, true\);?\n/.test(rest) && !/Issue\.error/.test(before) ? 'warning' : undefined) ??
          (/\[$/.test(before) ? /^[\s\S]{0,300}?Issue\.(warning|error|info)/.exec(rest)?.[1] : undefined) ??
          'error';
        seen.set(m[1], (seen.get(m[1]) ?? new Set()).add(severity));
      }
    }
    expect(seen.size).toBeGreaterThan(40);
    const wrong: string[] = [];
    for (const [code, severities] of seen) {
      // One code, one severity (seqspec:format#validation): every site's must be the registered one.
      const out = [...severities].map((s) =>
        svsIssue({ code, severity: s as Issue['severity'], path: '', message: '' }),
      );
      const registered = reg.get('issueCodes', out[0].code)?.severity;
      if (!registered) wrong.push(`${code} -> ${out[0].code}: not registered`);
      else if (!out.every((i) => i.severity === registered)) wrong.push(`${code} -> ${out[0].code}: vs ${registered}`);
    }
    expect(wrong).toEqual([]);
    // The table's codes for these, not near ones: a provider that cannot run is a warning (seqspec:format#layer-types).
    const as = (code: string) => svsIssue(Issue.error(code, '', ''));
    expect(as('model.numbering.invalid').code).toBe('svs.numbering.invalid');
    // ST1b-a: a removal of a key the set does not hold is the stream phase's (seqspec:format#validation), a warning
    // like the model's.
    expect(svsIssue(Issue.warning('model.key.unknown', '', ''))).toMatchObject({ code: 'svs.stream.ref' });
    expect(as('model.alignment.parent').code).toBe('svs.alignment.parent');
    expect([as('model.provider.lazy'), as('model.provider.alphabet')]).toMatchObject([
      { code: 'svs.provider.lazy', severity: 'warning' },
      { code: 'svs.provider.alphabet', severity: 'warning' },
    ]);
  });
});
