/**
 * Keeps the v6 migration records from going stale.
 *
 * - `.v6/plans/migration-map.json`: `{ "<v5 path>": "<v6 path>" | ["<v6 path>", ...] | null }`. Every non-null target
 *   must exist.
 * - `.v6/plans/migration-symbols.json`: `{ "<v5 path>": { "<Symbol>": "<v6 path>" | null } }`. Every non-null target
 *   must exist and export the symbol. A dotted symbol (`Namespace.Member`) is checked by its last segment, except
 *   `Namespace.BuiltIn`: that namespace value was removed (the type stays), so only the target file must exist.
 *
 * Export detection is a syntactic scan of top-level declarations, named re-exports and relative `export *`.
 */
import fs from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';
import ts from '@typescript/typescript6';

const repoRoot = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '../..');

/** Names exported by a TypeScript module (declarations, `export { }`, `export * as`, relative `export *`). */
export function exportNames(file, seen = new Set()) {
  const names = new Set();
  if (seen.has(file) || !fs.existsSync(file)) return names;
  seen.add(file);
  const source = ts.createSourceFile(
    file,
    fs.readFileSync(file, 'utf8'),
    ts.ScriptTarget.Latest,
    true,
    file.endsWith('x') ? ts.ScriptKind.TSX : ts.ScriptKind.TS,
  );
  for (const node of source.statements) {
    const modifiers = ts.canHaveModifiers(node) ? (ts.getModifiers(node) ?? []) : [];
    const exported = modifiers.some((m) => m.kind === ts.SyntaxKind.ExportKeyword);
    if (ts.isVariableStatement(node)) {
      if (exported) for (const d of node.declarationList.declarations) addBindingNames(d.name, names);
    } else if (
      exported &&
      node.name &&
      (ts.isFunctionDeclaration(node) ||
        ts.isClassDeclaration(node) ||
        ts.isInterfaceDeclaration(node) ||
        ts.isTypeAliasDeclaration(node) ||
        ts.isEnumDeclaration(node) ||
        ts.isModuleDeclaration(node))
    ) {
      names.add(node.name.text);
    } else if (ts.isExportDeclaration(node)) {
      const clause = node.exportClause;
      if (!clause) {
        const specifier = node.moduleSpecifier?.text;
        if (specifier?.startsWith('.')) {
          const target = resolveRelative(path.dirname(file), specifier);
          if (target) for (const name of exportNames(target, seen)) names.add(name);
        }
      } else if (ts.isNamedExports(clause)) for (const e of clause.elements) names.add(e.name.text);
      else names.add(clause.name.text);
    }
  }
  return names;
}

function addBindingNames(name, names) {
  if (ts.isIdentifier(name)) names.add(name.text);
  else for (const element of name.elements) if (!ts.isOmittedExpression(element)) addBindingNames(element.name, names);
}

function resolveRelative(dir, specifier) {
  const base = path.resolve(dir, specifier);
  const stem = base.replace(/\.[cm]?js$/u, '');
  return [`${stem}.ts`, `${stem}.tsx`, path.join(base, 'index.ts')].find((f) => fs.existsSync(f));
}

function readJson(file, errors) {
  try {
    return JSON.parse(fs.readFileSync(file, 'utf8'));
  } catch (error) {
    errors.push(`${file}: ${error.message}`);
  }
}

export function checkMigrationRecords({
  root = repoRoot,
  mapFile = path.join(root, '.v6/plans/migration-map.json'),
  symbolsFile = path.join(root, '.v6/plans/migration-symbols.json'),
} = {}) {
  const errors = [];
  const rel = (file) => path.relative(root, file).split(path.sep).join('/');
  const exists = (target) => fs.existsSync(path.join(root, target)) && fs.statSync(path.join(root, target)).isFile();

  const map = readJson(mapFile, errors);
  for (const [source, value] of Object.entries(map ?? {})) {
    const targets = Array.isArray(value) ? value : [value];
    if (Array.isArray(value) && (value.length < 2 || new Set(value).size !== value.length))
      errors.push(`${rel(mapFile)}: ${source}: an array needs at least two distinct targets`);
    for (const target of targets) {
      if (target === null && !Array.isArray(value)) continue;
      if (typeof target !== 'string')
        errors.push(`${rel(mapFile)}: ${source}: invalid target ${JSON.stringify(target)}`);
      else if (!exists(target)) errors.push(`${rel(mapFile)}: ${source}: target does not exist: ${target}`);
    }
  }

  const symbols = readJson(symbolsFile, errors);
  const exportsCache = new Map();
  for (const [source, entries] of Object.entries(symbols ?? {})) {
    for (const [symbol, target] of Object.entries(entries)) {
      const where = `${rel(symbolsFile)}: ${source}: ${symbol}`;
      if (target === null) continue;
      if (typeof target !== 'string') {
        errors.push(`${where}: invalid target ${JSON.stringify(target)}`);
        continue;
      }
      if (!exists(target)) {
        errors.push(`${where}: target does not exist: ${target}`);
        continue;
      }
      const segments = symbol.split('.');
      if (segments.length > 1 && segments.at(-1) === 'BuiltIn') continue;
      if (!exportsCache.has(target)) exportsCache.set(target, exportNames(path.join(root, target)));
      if (!exportsCache.get(target).has(segments.at(-1))) errors.push(`${where}: not exported by ${target}`);
    }
  }
  return { errors };
}
