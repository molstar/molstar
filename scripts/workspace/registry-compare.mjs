#!/usr/bin/env node
/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Acceptance comparison for the plugin-composition refactor (.v6/plans/plugin-composition.md, section 6): dumps what
 * the default spec and the Viewer register (`registry-dump.mjs`) and compares it with the step-0 baseline in
 * `.v6/baselines`. Providers must match by name and order, apart from the intentional differences of section 6.1, and
 * the baseline transformer ids must all still be registered.
 *
 * Usage: node scripts/workspace/registry-compare.mjs
 *
 * Requires compiled library output (`pnpm build:lib`).
 */

import { execFileSync } from 'node:child_process';
import { mkdtempSync, readFileSync, rmSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { dirname, join, resolve } from 'node:path';
import { fileURLToPath } from 'node:url';

const root = resolve(dirname(fileURLToPath(import.meta.url)), '../..');
const baselineDir = join(root, '.v6/baselines');

/** Section 6.1: the differences the comparison allows. They are removed from the current dump before comparing. */
const Allowed = {
  /** `external-structure` and `external-volume` are registered in all three scopes (as in 5.x); the prototype dropped them. */
  externalThemes: ['external-structure', 'external-volume'],
  /** The open-anything handler is the registered `DefaultDragAndDrop` fallback instead of being hard-coded. */
  dragAndDrop: ['open-files'],
};

/**
 * Particles are opt-in (`ParticlesRegistry`, breaking-v6-changes.md): the default spec registers none of them, so the
 * baseline's particle entries are taken out before comparing it; the Viewer composes `ParticlesRegistry` (as in 5.12)
 * and is compared with them.
 */
const ParticleTransformer =
  /^ms-plugin\.(particle|parse-artiatomi-em|parse-cryoet-data-portal-ndjson|parse-dynamo-tbl|parse-simularium)/;
const isParticleTheme = (name) => name.startsWith('particle-');

function withoutParticles(target) {
  const copy = structuredClone(target);
  copy.representations.particles = [];
  copy.themes.particles = { color: [], size: [] };
  for (const scope of Object.values(copy.themes)) {
    scope.color = scope.color.filter((n) => !isParticleTheme(n));
    scope.size = scope.size.filter((n) => !isParticleTheme(n));
  }
  copy.formats = copy.formats.filter((n) => !n.endsWith('_particles'));
  copy.animations = copy.animations.filter((n) => n !== 'built-in.animate-particle-trajectory');
  const data = copy.actions.data;
  data.all = data.all.filter((id) => !ParticleTransformer.test(id));
  for (const [from, ids] of Object.entries(data.byFromType)) {
    const kept = ids.filter((id) => !ParticleTransformer.test(id));
    if (kept.length) data.byFromType[from] = kept;
    else delete data.byFromType[from];
  }
  return copy;
}

const read = (dir, file) => JSON.parse(readFileSync(join(dir, file), 'utf8'));

function withoutAllowed(target) {
  const copy = structuredClone(target);
  for (const scope of Object.values(copy.themes)) {
    scope.color = scope.color.filter((n) => !Allowed.externalThemes.includes(n));
  }
  copy.dragAndDrop = copy.dragAndDrop.filter((n) => !Allowed.dragAndDrop.includes(n));
  return copy;
}

function diff(path, a, b, out) {
  if (JSON.stringify(a) === JSON.stringify(b)) return;
  if (a && b && typeof a === 'object' && typeof b === 'object' && !Array.isArray(a) && !Array.isArray(b)) {
    for (const k of new Set([...Object.keys(a), ...Object.keys(b)])) diff(`${path}.${k}`, a[k], b[k], out);
    return;
  }
  if (Array.isArray(a) && Array.isArray(b)) {
    const missing = a.filter((x) => !b.includes(x));
    const extra = b.filter((x) => !a.includes(x));
    const order = missing.length === 0 && extra.length === 0 ? ' (order differs)' : '';
    out.push(`${path}: missing [${missing}] extra [${extra}]${order}`);
    return;
  }
  out.push(`${path}: baseline ${JSON.stringify(a)}, current ${JSON.stringify(b)}`);
}

const out = mkdtempSync(join(tmpdir(), 'molstar-registry-'));
const problems = [];
try {
  execFileSync(process.execPath, [join(root, 'scripts/workspace/registry-dump.mjs'), '--out', out], {
    cwd: root,
    stdio: ['ignore', 'ignore', 'inherit'],
  });
  const baseline = read(baselineDir, 'registry.json');
  const current = read(out, 'registry.json');
  const baselineIds = read(baselineDir, 'transformer-ids.json');
  const currentIds = read(out, 'transformer-ids.json');

  for (const target of Object.keys(baseline.targets)) {
    const c = current.targets[target];
    if (!c) {
      problems.push(`${target}: target missing from the current dump`);
      continue;
    }
    const expected = target === 'viewer' ? baseline.targets[target] : withoutParticles(baseline.targets[target]);
    diff(target, expected, withoutAllowed(c), problems);

    // The allowed differences must actually be present, so an obsolete allowance is noticed.
    for (const [scope, t] of Object.entries(c.themes)) {
      // opt-in particles: no particles scope to check
      if (scope === 'particles') continue;
      for (const name of Allowed.externalThemes) {
        if (!t.color.includes(name)) problems.push(`${target}.themes.${scope}.color: expected '${name}'`);
      }
    }
    for (const name of Allowed.dragAndDrop) {
      if (!c.dragAndDrop.includes(name)) problems.push(`${target}.dragAndDrop: expected '${name}'`);
    }

    const have = new Set(currentIds.targets[target]);
    const lost = baselineIds.targets[target].filter((id) => !have.has(id) && !ParticleTransformer.test(id));
    if (lost.length > 0) problems.push(`${target}: transformer ids no longer registered: ${lost.join(', ')}`);
  }
} finally {
  rmSync(out, { recursive: true, force: true });
}

if (problems.length > 0) {
  console.error('The registries differ from the step-0 baseline beyond the plan section 6.1 differences:');
  for (const p of problems) console.error(`  ${p}`);
  process.exit(1);
}
console.log('Registries match the step-0 baseline apart from the plan section 6.1 differences.');
