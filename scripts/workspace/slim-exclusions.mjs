/**
 * The excluded-module list of the slim-plugin acceptance target (.v6/designs/plugin-composition.md §12), shared by the
 * import-graph check (`import-graph.mjs`) and the bundle check (`slim-bundle.mjs`).
 *
 * `slim-exclusions.json` pins the list as source-path globs relative to the repository root:
 *   - `entry`: the example entry whose value-import graph and bundles must not contain an excluded module;
 *   - `excluded`: named groups of path globs;
 *   - `knownLeaks`: excluded modules that the entry is known to reach and that need a design decision to remove. Each
 *     has a `path` (a glob), a `reason`, and a `decision` naming what has to be decided. They are the only exceptions,
 *     and an entry that no longer matches a reached module is an error (remove it).
 *
 * Globs support `**` (any number of directories), `*` and `?` (within a path segment) and `{a,b}` alternatives.
 */
import fs from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';

const repoRoot = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '../..');

export function loadExclusions(root = repoRoot) {
  return JSON.parse(fs.readFileSync(path.join(root, 'scripts/workspace/slim-exclusions.json'), 'utf8'));
}

/** Expands the first `{a,b}` group recursively. */
function expandBraces(glob) {
  const open = glob.indexOf('{');
  if (open < 0) return [glob];
  const close = glob.indexOf('}', open);
  if (close < 0) throw new Error(`Unbalanced '{' in glob: ${glob}`);
  const head = glob.slice(0, open);
  const tail = glob.slice(close + 1);
  return glob
    .slice(open + 1, close)
    .split(',')
    .flatMap((option) => expandBraces(head + option + tail));
}

function segmentToSource(segment) {
  let out = '';
  for (const ch of segment) {
    if (ch === '*') out += '[^/]*';
    else if (ch === '?') out += '[^/]';
    else out += ch.replace(/[.+^${}()|[\]\\]/g, '\\$&');
  }
  return out;
}

/** Compiles a glob to a RegExp matching whole posix paths. */
export function globToRegExp(glob) {
  const alternatives = expandBraces(glob).map((g) => {
    const segments = g.split('/');
    let source = '';
    segments.forEach((segment, i) => {
      const last = i === segments.length - 1;
      if (segment === '**') source += last ? '.*' : '(?:[^/]+/)*';
      else source += segmentToSource(segment) + (last ? '' : '/');
    });
    return source;
  });
  return new RegExp(`^(?:${alternatives.join('|')})$`);
}

function compile(patterns) {
  return patterns.map((glob) => ({ glob, regExp: globToRegExp(glob) }));
}

/** Validates the shape of the exclusions file; returns error messages. */
export function exclusionErrors(exclusions) {
  const errors = [];
  if (typeof exclusions.entry !== 'string') errors.push('slim-exclusions.json: entry must be a string');
  if (!Array.isArray(exclusions.excluded) || !exclusions.excluded.length)
    errors.push('slim-exclusions.json: excluded must be a non-empty array');
  for (const group of exclusions.excluded ?? []) {
    if (!group.name || !Array.isArray(group.paths) || !group.paths.length)
      errors.push(`slim-exclusions.json: invalid excluded group: ${JSON.stringify(group)}`);
  }
  for (const leak of exclusions.knownLeaks ?? []) {
    if (!leak.path || !leak.reason || !leak.decision)
      errors.push(`slim-exclusions.json: a known leak needs a path, a reason and a decision: ${JSON.stringify(leak)}`);
  }
  return errors;
}

/**
 * Splits the modules that match the exclusion list into real violations and known leaks.
 * `files` are repository-relative posix paths. Returns `{ violations: Map(file -> group name), known: Map(file -> leak),
 * unusedLeaks }`; `unusedLeaks` lists the known-leak entries that match none of `files`.
 */
export function findExcluded(files, exclusions) {
  const groups = exclusions.excluded.map((group) => ({ name: group.name, patterns: compile(group.paths) }));
  const leaks = (exclusions.knownLeaks ?? []).map((leak) => ({ leak, regExp: globToRegExp(leak.path) }));
  const violations = new Map();
  const known = new Map();
  const usedLeaks = new Set();
  for (const file of files) {
    const group = groups.find((g) => g.patterns.some((p) => p.regExp.test(file)));
    if (!group) continue;
    const leak = leaks.find((l) => l.regExp.test(file));
    if (leak) {
      known.set(file, leak.leak);
      usedLeaks.add(leak.leak);
    } else violations.set(file, group.name);
  }
  return { violations, known, unusedLeaks: leaks.map((l) => l.leak).filter((leak) => !usedLeaks.has(leak)) };
}
