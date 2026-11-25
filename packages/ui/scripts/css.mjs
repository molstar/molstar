/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Compiles every `.css` file under `src/` into a module in the `_generated/` folder next to it (`x.css` →
// `_generated/x.css.ts`) that exports its class names (for `*.module.css`) and registers its CSS with
// `styles/registry.ts` (styles §8). The generated modules are committed, so tsc, Vitest, bundlers and browsers need no
// CSS support.
//
//   node scripts/css.mjs            regenerate everything
//   node scripts/css.mjs --check    fail if a generated module is missing, stale or orphaned, or spacing is off the grid
//   node scripts/css.mjs --watch    regenerate on change
//
// A package built on @molstar/ui runs the same script on itself: `--root <package> --prefix <name>`. Its modules import
// the registry from `@molstar/ui/styles/registry` and name their classes `msui-<name>-<file>-<class>`.

import fs from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';
import { transform } from 'lightningcss';

const uiRoot = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '..');

const HEADER = fs.readFileSync(fileURLToPath(import.meta.url), 'utf8').split('*/')[0] + '*/';

const GENERATED = '_generated';

/** A package whose CSS is built: @molstar/ui itself (the default), or one built on it, with its class-name prefix. */
export function cssPackage(root = uiRoot, prefix = undefined) {
  root = path.resolve(root);
  const own = root === uiRoot;
  if (!own && !prefix) throw new Error(`${root}: a package built on @molstar/ui needs a class-name prefix (--prefix)`);
  return {
    root,
    srcRoot: path.join(root, 'src'),
    prefix: own ? undefined : prefix,
    /** Where the generated modules import `defineStyles` from: a path in @molstar/ui, its package export elsewhere. */
    registry: own ? path.join(root, 'src', 'styles', 'registry.js') : '@molstar/ui/styles/registry',
  };
}

const ui = cssPackage();

function outputOf(file) {
  return path.join(path.dirname(file), GENERATED, `${path.basename(file)}.ts`);
}

function listCss(dir, out = []) {
  for (const entry of fs.readdirSync(dir, { withFileTypes: true })) {
    const full = path.join(dir, entry.name);
    if (entry.isDirectory()) {
      if (entry.name !== GENERATED) listCss(full, out);
    } else if (entry.name.endsWith('.css')) out.push(full);
  }
  return out;
}

function listGenerated(dir, out = []) {
  for (const entry of fs.readdirSync(dir, { withFileTypes: true })) {
    const full = path.join(dir, entry.name);
    if (entry.isDirectory()) listGenerated(full, out);
    else if (path.basename(dir) === GENERATED && entry.name.endsWith('.css.ts')) out.push(full);
  }
  return out;
}

/** The class-name prefix of a module: its file name, with its folder when they differ; in a package built on the UI,
 * its file name after the package's prefix. */
function modulePrefix(pkg, file) {
  const base = path.basename(file, '.module.css');
  if (pkg.prefix) return `${pkg.prefix}-${base}`;
  const dir = path.basename(path.dirname(file));
  return base === dir ? base : `${dir}-${base}`;
}

function escapeTemplate(text) {
  return text.replace(/\\/g, '\\\\').replace(/`/g, '\\`').replace(/\$\{/g, '\\${');
}

/** The generated module for a CSS file. */
export function compile(file, pkg = ui) {
  const isModule = file.endsWith('.module.css');
  const local = path.relative(pkg.srcRoot, file).split(path.sep).join('/');
  const id = pkg.prefix ? `${pkg.prefix}/${local}` : local;
  const prefix = isModule ? modulePrefix(pkg, file) : undefined;
  const result = transform({
    filename: id,
    code: fs.readFileSync(file),
    // Grid area and container names stay global: they are layout contracts between modules.
    cssModules: isModule ? { pattern: `msui-${prefix}-[local]`, grid: false, container: false } : false,
    minify: false,
  });
  let css = result.code.toString().trim();
  // Names that already carry the `msui-` prefix are global (the primitives, shared keyframes): a module refers to them
  // as they are; its own names never start with it.
  if (isModule) css = css.replaceAll(`msui-${prefix}-msui-`, 'msui-');
  if (isModule) css = `@layer msui.components {\n${css}\n}`;
  const out = outputOf(file);
  let importPath = pkg.registry;
  if (path.isAbsolute(importPath)) {
    importPath = path.relative(path.dirname(out), importPath).split(path.sep).join('/');
    if (!importPath.startsWith('.')) importPath = `./${importPath}`;
  }
  const lines = [HEADER, '', `// Generated from ${path.basename(file)} by scripts/css.mjs; do not edit.`, ''];
  lines.push(`import { defineStyles } from '${importPath}';`, '');
  if (isModule) {
    const classes = Object.entries(result.exports ?? {})
      .filter(([name]) => !name.startsWith('msui-'))
      .sort(([a], [b]) => (a < b ? -1 : a > b ? 1 : 0))
      .map(([name, e]) => `  '${name}': '${e.name}',`);
    lines.push('export const classes = {', ...classes, '} as const;', '');
  }
  lines.push(`defineStyles('${id}', \`${escapeTemplate(css)}\n\`);`, '');
  return { out, content: lines.join('\n'), prefix };
}

const Spacing = /(?:^|[\s;{])(padding|margin|gap|row-gap|column-gap|padding-[a-z-]+|margin-[a-z-]+)\s*:([^;{}]*)/g;

/**
 * Spacing off the grid (styles §3): a padding, margin or gap in px that is not a multiple of 4 (1px hairlines aside),
 * unless the line says why (`off-grid: …`). Prefer the tokens (`--msui-space-*`, `--msui-inset`, `--msui-inner`).
 */
function spacingProblems(file, rel) {
  const problems = [];
  const lines = fs.readFileSync(file, 'utf8').split('\n');
  lines.forEach((line, i) => {
    if (line.includes('off-grid') || /^\s*(\/\*|\*)/.test(line)) return;
    for (const m of line.matchAll(Spacing)) {
      const value = m[2].replace(/var\([^()]*\)/g, '');
      for (const px of value.matchAll(/(?<![\w.-])(-?\d+(?:\.\d+)?)px\b/g)) {
        const n = Math.abs(Number(px[1]));
        if (n === 0 || n === 1 || n % 4 === 0) continue;
        problems.push(`${rel(file)}:${i + 1}: spacing off the 4px grid: ${m[0].trim()}`);
        break;
      }
    }
  });
  return problems;
}

/** Regenerates (or, with `check`, verifies) all generated modules; returns the problems found. */
export function generateAll({ check = false, pkg = ui } = {}) {
  const problems = [];
  const prefixes = new Map();
  const expected = new Set();
  const rel = (file) => path.relative(pkg.root, file);
  for (const file of listCss(pkg.srcRoot)) {
    if (check) problems.push(...spacingProblems(file, rel));
    const { out, content, prefix } = compile(file, pkg);
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
  for (const out of listGenerated(pkg.srcRoot)) {
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

function remove(out) {
  fs.rmSync(out, { force: true });
  const dir = path.dirname(out);
  if (fs.existsSync(dir) && fs.readdirSync(dir).length === 0) fs.rmdirSync(dir);
}

/** Regenerates one file's module when its CSS changes or goes away. */
export function updateFile(file, pkg = ui) {
  if (!file.startsWith(pkg.srcRoot + path.sep) || !file.endsWith('.css')) return;
  if (fs.existsSync(file)) {
    const { out, content } = compile(file, pkg);
    if (!fs.existsSync(out) || fs.readFileSync(out, 'utf8') !== content) write(out, content);
  } else {
    remove(outputOf(file));
  }
}

/** A Vite plugin (Storybook) that keeps the generated modules of these packages (@molstar/ui by default) current. */
export function cssModulesPlugin(packages = [ui]) {
  return {
    name: 'msui-css-modules',
    buildStart() {
      const problems = packages.flatMap((pkg) => generateAll({ pkg }));
      if (problems.length) this.error(problems.join('\n'));
    },
    configureServer(server) {
      const onChange = (file) => {
        try {
          for (const pkg of packages) updateFile(path.resolve(file), pkg);
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
  const arg = (name) => {
    const i = process.argv.indexOf(name);
    return i > 0 ? process.argv[i + 1] : undefined;
  };
  const pkg = cssPackage(arg('--root') ?? uiRoot, arg('--prefix'));
  const check = process.argv.includes('--check');
  const problems = generateAll({ check, pkg });
  if (problems.length) {
    console.error(problems.join('\n'));
    process.exit(1);
  }
  if (process.argv.includes('--watch')) {
    console.log('Watching src/**/*.css');
    fs.watch(pkg.srcRoot, { recursive: true }, (_event, name) => {
      if (!name) return;
      try {
        updateFile(path.join(pkg.srcRoot, name), pkg);
      } catch (e) {
        console.error(String(e));
      }
    });
  }
}
