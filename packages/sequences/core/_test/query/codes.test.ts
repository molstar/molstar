/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

// SeqQL's code registry (seqspec:seqql#errors): one severity per code, and the runtime's warnings and failures each
// under a code of their severity.
import { readdirSync, readFileSync } from 'node:fs';
import { join } from 'node:path';
import { fileURLToPath } from 'node:url';
import { describe, expect, it } from 'vitest';
import { queryIssue, SeqqlIssueCodes, seqqlSeverity } from '../../src/query/index.js';

const here = (p: string) => fileURLToPath(new URL(p, import.meta.url));
const sources = (dir: string): string[] =>
  readdirSync(dir, { recursive: true, encoding: 'utf8' })
    .filter((f) => f.endsWith('.ts'))
    .map((f) => readFileSync(join(dir, f), 'utf8'));

describe("SeqQL's issue codes", () => {
  it('lists each code once, snake_case, and an issue takes the registered severity', () => {
    const codes = SeqqlIssueCodes.map((c) => c.code);
    expect(new Set(codes).size).toBe(codes.length);
    expect(codes.filter((c) => !/^seqql(\.[a-z_]+)+$/.test(c))).toEqual([]);
    // An empty result explained is a note (seqspec:seqql#empty-results); a result that may be incomplete, a warning.
    expect(queryIssue('seqql.empty', 'm').severity).toBe('info');
    expect(queryIssue('seqql.lazy', 'm').severity).toBe('warning');
    expect(queryIssue('seqql.type', 'm').severity).toBe('error');
    // An extension's own code, unlisted: an error.
    expect(seqqlSeverity('acme.own')).toBe('error');
  });

  it('every code the sources name is registered; the runtime warns under warnings or notes and fails under errors', () => {
    const text = [...sources(here('../../src/query')), ...sources(here('../../src/model/query'))].join('\n');
    const named = new Set([...text.matchAll(/["'](seqql\.[a-z_]+(?:\.[a-z_]+)*)["']/g)].map((m) => m[1]));
    expect([...named].filter((c) => !SeqqlIssueCodes.some((e) => e.code === c))).toEqual([]);
    const wrong: string[] = [];
    for (const [, how, code = ''] of text.matchAll(/\.(warn|fail)\(\s*["'](seqql\.[a-z_.]+)["']/g)) {
      const want = how === 'warn' ? ['warning', 'info'] : ['error'];
      if (!want.includes(seqqlSeverity(code))) wrong.push(`${how} ${code}: registered ${seqqlSeverity(code)}`);
    }
    expect(wrong).toEqual([]);
  });
});
