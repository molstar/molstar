/**
 * Value-import graph check for the plugin composition work (.v6/designs/plugin-composition.md §8, §9).
 *
 * Classification comes from `catalogs.json` (catalog modules, default-composition modules, base entry points, the
 * package layer order, and an explicit allowlist of planned violations). Apps are everything under `apps/`,
 * `examples/`, `cli/`, `servers/` and `smoke/`; any other module is a base module.
 *
 * Rules:
 *   a. A value import of a catalog module is allowed only from a catalog module, a default-composition module, or an app.
 *   b. Base entry points must not reach a default-composition module by value (transitively).
 *   c. Base entry points must not reach any catalog module by value (transitively). The search does not continue past a
 *      catalog or default-composition module: the first one on each path is reported (default-composition modules by b).
 *   d. A type-only import must not point to a higher package (manifest `layers` order).
 *   e. The slim example entry (`slim-exclusions.json`, plan step 4) must not reach an excluded module by value
 *      (transitively), except for the known leaks listed there. The first excluded module on each path is reported
 *      with its import chain.
 *   f. The boundary modules of `slim-exclusions.json` (external color themes and structure queries, spec §5.4) must not
 *      reach a transform module, `PluginContext` or any catalog / default-composition module by value.
 *
 * Known violations live in the manifest allowlist with a reason and the plan step that removes them. The check fails
 * on any violation that is not allowlisted and on allowlist entries that no longer match a violation.
 */
import fs from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';
import ts from '@typescript/typescript6';
import { resolveExport } from './exports.mjs';
import { importsFrom } from './imports.mjs';
import { exclusionErrors, findExcluded, globToRegExp, loadExclusions } from './slim-exclusions.mjs';

const repoRoot = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '../..');
const APP_ROOTS = ['apps/', 'examples/', 'cli/', 'servers/', 'smoke/'];
const SCAN_ROOTS = ['packages/', 'extensions/', ...APP_ROOTS];
const SKIP_DIRS = new Set(['node_modules', 'lib', '_test', '.claude', 'dist', 'build', '.git']);
const SOURCE_EXTENSIONS = ['.ts', '.tsx'];

export function loadManifest(root = repoRoot) {
  return JSON.parse(fs.readFileSync(path.join(root, 'scripts/workspace/catalogs.json'), 'utf8'));
}

const toPosix = (file) => file.split(path.sep).join('/');
const isSource = (name) => /\.tsx?$/.test(name) && !name.endsWith('.d.ts') && !/\.test\.tsx?$/.test(name);

function walk(dir, files = []) {
  if (!fs.existsSync(dir)) return files;
  for (const ent of fs.readdirSync(dir, { withFileTypes: true })) {
    if (ent.isDirectory()) {
      if (!SKIP_DIRS.has(ent.name)) walk(path.join(dir, ent.name), files);
    } else if (ent.isFile() && isSource(ent.name)) files.push(path.join(dir, ent.name));
  }
  return files;
}

function isFile(file) {
  return fs.existsSync(file) && fs.statSync(file).isFile();
}

/** Resolves a file path (possibly with a `.js` extension or without one) to an existing TS source file. */
function resolveSourceFile(file) {
  const candidates = [];
  if (/\.[cm]?jsx?$/.test(file)) {
    const base = file.replace(/\.[cm]?jsx?$/, '');
    for (const ext of SOURCE_EXTENSIONS) candidates.push(base + ext);
  } else if (SOURCE_EXTENSIONS.some((ext) => file.endsWith(ext))) candidates.push(file);
  else {
    for (const ext of SOURCE_EXTENSIONS) candidates.push(file + ext);
    for (const ext of SOURCE_EXTENSIONS) candidates.push(path.join(file, 'index' + ext));
  }
  return candidates.find((candidate) => isFile(candidate) && isSource(path.basename(candidate)));
}

function packageExportTarget(exports, key) {
  const value = resolveExport(exports ?? {}, key);
  if (typeof value === 'string') return value;
  if (value && typeof value === 'object') {
    if (typeof value['molstar-src'] === 'string') return value['molstar-src'];
    for (const field of ['import', 'default', 'types']) {
      if (typeof value[field] === 'string') return value[field].replace(/^\.\/lib\//u, './src/');
    }
  }
  return undefined;
}

function pathPrefixMatch(file, dir) {
  return file === dir || file.startsWith(dir + path.sep);
}

/** Builds the value-import and type-import edges of every scanned source file. */
export function buildImportGraph({ root, packages }) {
  const scanned = packages.filter((pkg) => SCAN_ROOTS.some((prefix) => pkg.path.startsWith(prefix)));
  const byName = new Map(packages.map((pkg) => [pkg.name, pkg]));
  const manifests = new Map();
  const manifestOf = (pkg) => {
    if (!manifests.has(pkg.name)) {
      const file = path.join(root, pkg.path, 'package.json');
      manifests.set(pkg.name, isFile(file) ? JSON.parse(fs.readFileSync(file, 'utf8')) : {});
    }
    return manifests.get(pkg.name);
  };
  const packageDirs = packages.map((pkg) => ({ pkg, dir: path.resolve(root, pkg.path) }));
  const ownerOf = (file) => {
    let best;
    for (const entry of packageDirs)
      if (pathPrefixMatch(file, entry.dir) && (!best || entry.dir.length > best.dir.length)) best = entry;
    return best?.pkg;
  };
  const rel = (file) => toPosix(path.relative(root, file));

  function resolveSpecifier(from, specifier) {
    if (specifier.startsWith('.')) return resolveSourceFile(path.resolve(path.dirname(from), specifier));
    if (!specifier.startsWith('@molstar/')) return undefined;
    const parts = specifier.split('/');
    const pkg = byName.get(parts.slice(0, 2).join('/'));
    if (!pkg) return undefined;
    const key = parts.length > 2 ? `./${parts.slice(2).join('/')}` : '.';
    const target = packageExportTarget(manifestOf(pkg).exports, key);
    if (!target) return undefined;
    return resolveSourceFile(path.resolve(root, pkg.path, target));
  }

  const files = new Set();
  const value = new Map();
  const types = [];
  for (const pkg of scanned)
    for (const file of walk(path.join(root, pkg.path))) {
      // Nested packages (e.g. packages/mvs/runtime inside packages/mvs) own their own files.
      if (ownerOf(file)?.name === pkg.name) files.add(file);
    }
  for (const file of files) {
    const from = rel(file);
    const text = fs.readFileSync(file, 'utf8');
    const source = ts.createSourceFile(
      file,
      text,
      ts.ScriptTarget.Latest,
      true,
      file.endsWith('x') ? ts.ScriptKind.TSX : ts.ScriptKind.TS,
    );
    const edges = new Set();
    for (const imported of importsFrom(file, source)) {
      const target = resolveSpecifier(file, imported.specifier);
      if (!target) continue;
      const to = rel(target);
      if (imported.typeOnly) {
        const fromPkg = ownerOf(file);
        const toPkg = ownerOf(target);
        if (fromPkg && toPkg && fromPkg.name !== toPkg.name)
          types.push({ from, to, fromPkg: fromPkg.name, toPkg: toPkg.name });
      }
      // Only `import type` is erased (verbatimModuleSyntax); `import { type A }` still evaluates the module.
      if (!imported.erased && to !== from) edges.add(to);
    }
    value.set(from, edges);
  }
  return { files: new Set([...files].map(rel)), value, types };
}

export function classify(file, manifest) {
  if (manifest.catalogs.includes(file)) return 'catalog';
  if (manifest.defaultComposition.includes(file)) return 'default-composition';
  if (APP_ROOTS.some((prefix) => file.startsWith(prefix))) return 'app';
  return 'base';
}

function shortestPath(value, start, isGoal, isBarrier = () => false) {
  const previous = new Map([[start, undefined]]);
  const queue = [start];
  const found = new Map();
  for (let i = 0; i < queue.length; i++) {
    const current = queue[i];
    if (current !== start && isBarrier(current)) continue;
    for (const next of value.get(current) ?? []) {
      if (previous.has(next)) continue;
      previous.set(next, current);
      if (isGoal(next)) {
        const chain = [];
        for (let node = next; node !== undefined; node = previous.get(node)) chain.unshift(node);
        found.set(next, chain);
      }
      queue.push(next);
    }
  }
  return found;
}

/** Rules e and f (see the header): the slim example's excluded modules and the query/theme boundary. */
export function findSlimViolations(graph, manifest, exclusions) {
  const violations = [];
  const files = [...graph.files];
  const { violations: excluded } = findExcluded(files, exclusions);
  if (!graph.files.has(exclusions.entry))
    return [
      {
        rule: 'e',
        from: exclusions.entry,
        to: exclusions.entry,
        message: `slim example entry does not exist or is not scanned: ${exclusions.entry}`,
      },
    ];
  // Known leaks are not barriers: what they import is reported too, so it is a conscious decision to list it.
  const reached = shortestPath(
    graph.value,
    exclusions.entry,
    (file) => excluded.has(file),
    (file) => excluded.has(file),
  );
  for (const [to, chain] of reached)
    violations.push({
      rule: 'e',
      from: exclusions.entry,
      to,
      message: `slim example ${exclusions.entry} reaches excluded module ${to} (${excluded.get(to)}) by value: ${chain.join(' -> ')}`,
    });
  const everything = shortestPath(graph.value, exclusions.entry, () => true);
  for (const leak of exclusions.knownLeaks ?? []) {
    const re = globToRegExp(leak.path);
    if (![...everything.keys()].some((file) => re.test(file)))
      violations.push({
        rule: 'e',
        from: exclusions.entry,
        to: leak.path,
        message: `slim-exclusions.json: stale known leak (the slim example no longer reaches it, remove it): ${leak.path}`,
      });
  }

  const boundary = exclusions.boundary;
  if (boundary) {
    const forbidden = Object.entries(boundary.forbidden ?? {}).map(([name, globs]) => ({
      name,
      regExps: globs.map(globToRegExp),
    }));
    const forbiddenName = (file) => {
      if (manifest.catalogs.includes(file)) return 'catalog';
      if (manifest.defaultComposition.includes(file)) return 'default composition';
      return forbidden.find((f) => f.regExps.some((re) => re.test(file)))?.name;
    };
    const sources = files.filter(
      (file) =>
        boundary.sources.some((glob) => globToRegExp(glob).test(file)) &&
        // a catalog in the boundary directory is the catalog itself
        !manifest.catalogs.includes(file),
    );
    for (const glob of boundary.sources) {
      const re = globToRegExp(glob);
      if (!files.some((file) => re.test(file)))
        violations.push({
          rule: 'f',
          from: glob,
          to: glob,
          message: `slim-exclusions.json: boundary source matches no module: ${glob}`,
        });
    }
    for (const source of sources) {
      const found = shortestPath(
        graph.value,
        source,
        (file) => !!forbiddenName(file),
        (file) => !!forbiddenName(file),
      );
      for (const [to, chain] of found)
        violations.push({
          rule: 'f',
          from: source,
          to,
          message: `${source} reaches ${forbiddenName(to)} module ${to} by value: ${chain.join(' -> ')}`,
        });
    }
  }
  return violations;
}

/** Finds all rule violations without applying the allowlist. Each has a stable `key` used by the allowlist. */
export function findViolations(graph, manifest, exclusions) {
  const violations = [];
  for (const [from, edges] of graph.value) {
    if (['catalog', 'default-composition', 'app'].includes(classify(from, manifest))) continue;
    for (const to of edges)
      if (manifest.catalogs.includes(to))
        violations.push({
          rule: 'a',
          from,
          to,
          message: `${from} value-imports catalog module ${to} (only catalogs, default-composition modules and apps may)`,
        });
  }
  for (const entry of manifest.baseEntryPoints) {
    const reached = shortestPath(graph.value, entry, (file) => manifest.defaultComposition.includes(file));
    for (const [to, chain] of reached)
      violations.push({
        rule: 'b',
        from: entry,
        to,
        message: `base entry point ${entry} reaches default-composition module ${to} by value: ${chain.join(' -> ')}`,
      });
    const catalogs = shortestPath(
      graph.value,
      entry,
      (file) => manifest.catalogs.includes(file),
      (file) => manifest.catalogs.includes(file) || manifest.defaultComposition.includes(file),
    );
    for (const [to, chain] of catalogs)
      violations.push({
        rule: 'c',
        from: entry,
        to,
        message: `base entry point ${entry} reaches catalog module ${to} by value: ${chain.join(' -> ')}`,
      });
  }
  const layer = new Map(manifest.layers.map((name, index) => [name, index]));
  for (const edge of graph.types) {
    if (layer.has(edge.fromPkg) && layer.has(edge.toPkg) && layer.get(edge.toPkg) > layer.get(edge.fromPkg))
      violations.push({
        rule: 'd',
        from: edge.from,
        to: edge.to,
        message: `${edge.from} (${edge.fromPkg}) type-imports ${edge.to} from the higher package ${edge.toPkg}`,
      });
  }
  if (exclusions) violations.push(...findSlimViolations(graph, manifest, exclusions));
  const keyOf = (violation) => `${violation.rule} ${violation.from} -> ${violation.to}`;
  for (const violation of violations) violation.key = keyOf(violation);
  return violations;
}

function manifestErrors(root, graph, manifest) {
  const errors = [];
  const pending = new Set(manifest.pendingDefaultComposition ?? []);
  const lists = {
    catalogs: manifest.catalogs,
    defaultComposition: manifest.defaultComposition,
    baseEntryPoints: manifest.baseEntryPoints,
  };
  for (const [name, list] of Object.entries(lists)) {
    if (!Array.isArray(list)) {
      errors.push(`catalogs.json: ${name} must be an array`);
      continue;
    }
    if (new Set(list).size !== list.length) errors.push(`catalogs.json: ${name} contains duplicates`);
    for (const file of list)
      if (!graph.files.has(file)) errors.push(`catalogs.json: ${name} entry does not exist or is not scanned: ${file}`);
  }
  for (const file of pending)
    if (isFile(path.join(root, file)))
      errors.push(`catalogs.json: pendingDefaultComposition ${file} now exists; move it to defaultComposition`);
  for (const file of manifest.catalogs ?? [])
    if ((manifest.defaultComposition ?? []).includes(file))
      errors.push(`catalogs.json: ${file} is listed as both catalog and default-composition module`);
  for (const file of graph.files)
    if (path.posix.basename(file) === 'catalog.ts' && !(manifest.catalogs ?? []).includes(file))
      errors.push(`catalogs.json: catalog module is missing from the manifest: ${file}`);
  return errors;
}

/** The slim-plugin exclusions, or `undefined` when the workspace has none (rules e and f are then skipped). */
function loadOptionalExclusions(root) {
  return fs.existsSync(path.join(root, 'scripts/workspace/slim-exclusions.json')) ? loadExclusions(root) : undefined;
}

/** Runs the whole check and returns the list of error messages (empty when the graph is valid). */
export function checkImportGraph({
  root = repoRoot,
  packages,
  manifest = loadManifest(root),
  exclusions = loadOptionalExclusions(root),
}) {
  const graph = buildImportGraph({ root, packages });
  const errors = manifestErrors(root, graph, manifest);
  if (exclusions) errors.push(...exclusionErrors(exclusions));
  const violations = findViolations(graph, manifest, exclusions);
  const allowlist = manifest.allowlist ?? [];
  const entries = new Map();
  for (const entry of allowlist) {
    const key = `${entry.rule} ${entry.from} -> ${entry.to}`;
    if (!entry.reason || !entry.step) errors.push(`catalogs.json: allowlist entry needs a reason and a step: ${key}`);
    if (entries.has(key)) errors.push(`catalogs.json: duplicate allowlist entry: ${key}`);
    entries.set(key, entry);
  }
  const used = new Set();
  for (const violation of violations) {
    if (entries.has(violation.key)) used.add(violation.key);
    else errors.push(`import-graph rule ${violation.rule}: ${violation.message}`);
  }
  for (const key of entries.keys())
    if (!used.has(key)) errors.push(`import-graph: stale allowlist entry (the violation is gone, remove it): ${key}`);
  return { errors, violations, graph };
}

if (process.argv[1] && path.resolve(process.argv[1]) === fileURLToPath(import.meta.url)) {
  const started = Date.now();
  const inventory = JSON.parse(fs.readFileSync(path.join(repoRoot, 'scripts/workspace/inventory.json'), 'utf8'));
  const { errors, violations, graph } = checkImportGraph({ packages: inventory.packages ?? [] });
  const verbose = process.argv.includes('--list');
  if (verbose) for (const violation of violations) console.log(violation.key);
  if (errors.length) {
    console.error(errors.join('\n'));
    process.exitCode = 1;
  } else
    console.log(
      `Import graph is valid (${graph.files.size} modules, ${violations.length} allowlisted violations, ${Date.now() - started} ms).`,
    );
}
