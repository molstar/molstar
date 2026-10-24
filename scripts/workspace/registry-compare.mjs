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
 *
 * @author David Sehnal <david.sehnal@gmail.com>
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
  /** Formats new in v6: NAMD extended system files (periodic boxes for coordinates). */
  formats: ['xsc'],
  /**
   * Parsers take String, Binary and Response Data (`Download` and `ReadFile` give Response Data, `.v6/data-response`):
   * the data actions that list `Response Data` are taken out of the String and Binary Data lists of both dumps.
   */
  responseData: 'Response Data',
};

/** Section 6.1, the response data row: the parsers the current dump lists for Response Data, out of both dumps. */
function withoutResponseData(expected, current) {
  const e = structuredClone(expected);
  const c = structuredClone(current);
  const readers = new Set(c.actions.data.byFromType[Allowed.responseData] ?? []);
  delete c.actions.data.byFromType[Allowed.responseData];
  for (const t of [e, c]) {
    for (const from of ['String Data', 'Binary Data']) {
      const ids = t.actions.data.byFromType[from];
      if (!ids) continue;
      const kept = ids.filter((id) => !readers.has(id));
      if (kept.length) t.actions.data.byFromType[from] = kept;
      else delete t.actions.data.byFromType[from];
    }
  }
  return [e, c];
}

/**
 * Particles are opt-in (`ParticlesRegistry`, breaking-v6-changes.md): the default spec registers none of them, so the
 * baseline's particle entries are taken out before comparing it. The Viewer composes `ParticlesRegistry` (as in 5.12)
 * and is compared with them, as that registry adds them (`asComposed`): the particle themes in the particles scope only
 * (not the structure and volume scopes, where they never apply), and the particle formats and animation as a block of
 * their own, so only the order within the particle entries and within the others is compared (a format is chosen by
 * its priority: registration order only breaks ties).
 */
const ParticleTransformer =
  /^ms-plugin\.(particle|parse-artiatomi-em|parse-cryoet-data-portal-ndjson|parse-dynamo-tbl|parse-simularium)/;
const isParticleTheme = (name) => name.startsWith('particle-');
const isParticleFormat = (name) => name.endsWith('_particles');
const isParticleAnimation = (name) => name === 'built-in.animate-particle-trajectory';

/** The others in their order, then the particle entries in theirs */
const particlesLast = (xs, isParticle) => [...xs.filter((x) => !isParticle(x)), ...xs.filter(isParticle)];

/** A Viewer registry with its particle formats and animation as one block (the current one, or the baseline's) */
function particleBlocks(target) {
  const copy = structuredClone(target);
  copy.formats = particlesLast(copy.formats, isParticleFormat);
  copy.animations = particlesLast(copy.animations, isParticleAnimation);
  return copy;
}

/** The baseline as `[...DefaultRegistry, ...ParticlesRegistry]` registers it */
function asComposed(target) {
  const copy = particleBlocks(target);
  for (const [name, scope] of Object.entries(copy.themes)) {
    if (name === 'particles') continue;
    scope.color = scope.color.filter((n) => !isParticleTheme(n));
    scope.size = scope.size.filter((n) => !isParticleTheme(n));
  }
  return copy;
}

function withoutParticles(target) {
  const copy = structuredClone(target);
  copy.representations.particles = [];
  copy.themes.particles = { color: [], size: [] };
  for (const scope of Object.values(copy.themes)) {
    scope.color = scope.color.filter((n) => !isParticleTheme(n));
    scope.size = scope.size.filter((n) => !isParticleTheme(n));
  }
  copy.formats = copy.formats.filter((n) => !isParticleFormat(n));
  copy.animations = copy.animations.filter((n) => !isParticleAnimation(n));
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
  copy.formats = copy.formats.filter((n) => !Allowed.formats.includes(n));
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
    const viewer = target === 'viewer';
    const [expected, compared] = withoutResponseData(
      viewer ? asComposed(baseline.targets[target]) : withoutParticles(baseline.targets[target]),
      withoutAllowed(viewer ? particleBlocks(c) : c),
    );
    diff(target, expected, compared, problems);

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
    for (const name of Allowed.formats) {
      if (!c.formats.includes(name)) problems.push(`${target}.formats: expected '${name}'`);
    }
    if (!c.actions.data.byFromType[Allowed.responseData]?.includes('ms-plugin.parse-cif')) {
      problems.push(`${target}.actions.data.byFromType: expected '${Allowed.responseData}' with 'ms-plugin.parse-cif'`);
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
