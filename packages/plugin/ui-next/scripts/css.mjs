/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Compiles every `.css` file under `src/` into a module in the `_generated/` folder next to it (`x.css` →
// `_generated/x.css.ts`) that exports its class names (for `*.module.css`) and registers its CSS with
// `styles/registry.ts`. The generated modules are committed, so tsc, Jest, bundlers and browsers need no CSS support.
//
//   node scripts/css.mjs            regenerate everything
//   node scripts/css.mjs --check    fail if a generated module is missing, stale or orphaned
//   node scripts/css.mjs --watch    regenerate on change

import fs from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';
import { transform } from 'lightningcss';

const packageRoot = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '..');
const srcRoot = path.join(packageRoot, 'src');
const registry = path.join(srcRoot, 'styles', 'registry.js');

const HEADER = `/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */`;

const GENERATED = '_generated';

/** Where a CSS file's module goes. */
function outputOf(file) {
  return path.join(path.dirname(file), GENERATED, `${path.basename(file)}.ts`);
}

function listCss(dir = srcRoot, out = []) {
  for (const entry of fs.readdirSync(dir, { withFileTypes: true })) {
    const full = path.join(dir, entry.name);
    if (entry.isDirectory()) {
      if (entry.name !== GENERATED) listCss(full, out);
    } else if (entry.name.endsWith('.css')) out.push(full);
  }
  return out;
}

function listGenerated(dir = srcRoot, out = []) {
  for (const entry of fs.readdirSync(dir, { withFileTypes: true })) {
    const full = path.join(dir, entry.name);
    if (entry.isDirectory()) listGenerated(full, out);
    else if (path.basename(dir) === GENERATED && entry.name.endsWith('.css.ts')) out.push(full);
  }
  return out;
}

/** The class-name prefix of a module: its file name, with its folder when they differ (`tree`, `scene-row`). */
function modulePrefix(file) {
  const base = path.basename(file, '.module.css');
  const dir = path.basename(path.dirname(file));
  return base === dir ? base : `${dir}-${base}`;
}

function escapeTemplate(text) {
  return text.replace(/\\/g, '\\\\').replace(/`/g, '\\`').replace(/\$\{/g, '\\${');
}

/** The generated module for a CSS file. */
export function compile(file) {
  const isModule = file.endsWith('.module.css');
  const id = path.relative(srcRoot, file).split(path.sep).join('/');
  const result = transform({
    filename: id,
    code: fs.readFileSync(file),
    // Grid area and container names stay global: they are layout contracts between modules (shell grid, docks,
    // the canvas container).
    cssModules: isModule ? { pattern: `msp-${modulePrefix(file)}-[local]`, grid: false, container: false } : false,
    minify: false,
  });
  let css = result.code.toString().trim();
  if (isModule) css = `@layer msp.components {\n${css}\n}`;
  const out = outputOf(file);
  let importPath = path.relative(path.dirname(out), registry).split(path.sep).join('/');
  if (!importPath.startsWith('.')) importPath = `./${importPath}`;
  const lines = [HEADER, '', `// Generated from ${path.basename(file)} by scripts/css.mjs; do not edit.`, ''];
  lines.push(`import { defineStyles } from '${importPath}';`, '');
  if (isModule) {
    const classes = Object.entries(result.exports ?? {})
      .sort(([a], [b]) => (a < b ? -1 : a > b ? 1 : 0))
      .map(([local, e]) => `  '${local}': '${e.name}',`);
    lines.push('export const classes = {', ...classes, '} as const;', '');
  }
  lines.push(`defineStyles('${id}', \`${escapeTemplate(css)}\n\`);`, '');
  return { out, content: lines.join('\n'), prefix: isModule ? modulePrefix(file) : undefined };
}

/** Regenerates (or, with `check`, verifies) all generated modules; returns the problems found. */
export function generateAll({ check = false } = {}) {
  const problems = [];
  const prefixes = new Map();
  const expected = new Set();
  for (const file of listCss()) {
    const { out, content, prefix } = compile(file);
    expected.add(out);
    if (prefix) {
      const other = prefixes.get(prefix);
      if (other) problems.push(`${rel(file)}: class prefix "${prefix}" is also used by ${rel(other)}`);
      prefixes.set(prefix, file);
    }
    const current = fs.existsSync(out) ? fs.readFileSync(out, 'utf8') : undefined;
    if (current === content) continue;
    if (check) problems.push(`${rel(out)} is ${current === undefined ? 'missing' : 'stale'}`);
    else write(out, content);
  }
  for (const out of listGenerated()) {
    if (expected.has(out)) continue;
    if (check) problems.push(`${rel(out)} has no CSS file`);
    else remove(out);
  }
  return problems;
}

function write(out, content) {
  fs.mkdirSync(path.dirname(out), { recursive: true });
  fs.writeFileSync(out, content);
}

/** Removes a generated module, and its `_generated` folder once empty. */
function remove(out) {
  fs.rmSync(out, { force: true });
  const dir = path.dirname(out);
  if (fs.existsSync(dir) && fs.readdirSync(dir).length === 0) fs.rmdirSync(dir);
}

function rel(file) {
  return path.relative(packageRoot, file);
}

/** Regenerates one file's module when its CSS changes or goes away. */
export function updateFile(file) {
  if (!file.startsWith(srcRoot) || !file.endsWith('.css')) return;
  if (fs.existsSync(file)) {
    const { out, content } = compile(file);
    if (!fs.existsSync(out) || fs.readFileSync(out, 'utf8') !== content) write(out, content);
  } else {
    remove(outputOf(file));
  }
}

/** A Vite plugin (Storybook) that keeps the generated modules current. */
export function cssModulesPlugin() {
  return {
    name: 'msp-css-modules',
    buildStart() {
      const problems = generateAll();
      if (problems.length) this.error(problems.join('\n'));
    },
    configureServer(server) {
      const onChange = (file) => {
        try {
          updateFile(path.resolve(file));
        } catch (e) {
          server.config.logger.error(String(e));
        }
      };
      server.watcher.on('add', onChange);
      server.watcher.on('change', onChange);
      server.watcher.on('unlink', onChange);
    },
  };
}

if (process.argv[1] === fileURLToPath(import.meta.url)) {
  const check = process.argv.includes('--check');
  const problems = generateAll({ check });
  if (problems.length) {
    console.error(problems.join('\n'));
    process.exit(1);
  }
  if (process.argv.includes('--watch')) {
    console.log('Watching src/**/*.css');
    fs.watch(srcRoot, { recursive: true }, (_event, name) => {
      if (!name) return;
      try {
        updateFile(path.join(srcRoot, name));
      } catch (e) {
        console.error(String(e));
      }
    });
  }
}
