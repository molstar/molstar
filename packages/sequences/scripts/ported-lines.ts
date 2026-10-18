/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// New production lines per phase (the size budgets of packages/sequences/AGENTS.md). Production lines are as in the
// repository's line counts: not blank, not comments. Counted over the files under packages/sequences/*/src (any depth) that
// are added or changed since --base (default HEAD: staged, unstaged and untracked).
//   - A file that existed at the base counts its added lines as new (a longest-common-subsequence diff of lines).
//   - Any other file is all new.
// (The name stays from when a rebuilt file was diffed against its old source and counted as ported.)
// Usage: node packages/sequences/scripts/ported-lines.ts [--base <rev>]
import { execFileSync } from 'node:child_process';
import { existsSync, readFileSync } from 'node:fs';
import { join, relative } from 'node:path';

const SEQSTAR = join(import.meta.dirname, '..');
const root = execFileSync('git', ['rev-parse', '--show-toplevel'], { cwd: SEQSTAR, encoding: 'utf8' }).trim();
/** Seq*'s folder as git paths start: "packages/sequences/". */
const PREFIX = relative(root, SEQSTAR) ? `${relative(root, SEQSTAR)}/` : '';
const WATCHED = new RegExp(`^${PREFIX.replace(/[/.]/g, '\\$&')}[^/]+/src/.*\\.(tsx?|css)$`);

/** The production lines of a text, trimmed: not blank, not `//` comments, not `/* … *​/` blocks. */
export function productionLines(text: string): string[] {
  const out: string[] = [];
  let block = false;
  for (const raw of text.split('\n')) {
    const line = raw.trim();
    if (block) {
      block = !line.includes('*/');
      continue;
    }
    if (!line || line.startsWith('//')) continue;
    if (line.startsWith('/*')) {
      block = !line.includes('*/');
      continue;
    }
    out.push(line);
  }
  return out;
}

/** The length of the longest common subsequence of two line lists. */
export function common(a: readonly string[], b: readonly string[]): number {
  let prev = new Int32Array(b.length + 1);
  let cur = new Int32Array(b.length + 1);
  for (let i = 1; i <= a.length; i++) {
    for (let j = 1; j <= b.length; j++) {
      cur[j] = a[i - 1] === b[j - 1] ? prev[j - 1] + 1 : Math.max(prev[j], cur[j - 1]);
    }
    [prev, cur] = [cur, prev];
  }
  return prev[b.length];
}

/** The production lines of a file that its base lacks: the new ones. */
export function added(file: string, base: string): number {
  const lines = productionLines(file);
  return lines.length - common(productionLines(base), lines);
}

const git = (...args: string[]) =>
  execFileSync('git', args, { cwd: root, maxBuffer: 1 << 28, stdio: ['ignore', 'pipe', 'ignore'] }).toString();

if (import.meta.main) {
  const i = process.argv.indexOf('--base');
  const base = i > 0 ? process.argv[i + 1] : 'HEAD';
  const files = [
    ...git('diff', '--name-only', base).split('\n'),
    ...git('ls-files', '-o', '--exclude-standard').split('\n'),
  ].filter((f, k, all) => WATCHED.test(f) && all.indexOf(f) === k && existsSync(join(root, f)));
  let total = 0;
  for (const f of files.sort()) {
    const text = readFileSync(join(root, f), 'utf8');
    let at = '';
    try {
      at = git('show', `${base}:${f}`);
    } catch {}
    const n = added(text, at);
    total += n;
    console.log(`${String(n).padStart(6)} new  ${relative('.', f)}`);
  }
  console.log(`${String(total).padStart(6)} new  total, ${files.length} files`);
}
