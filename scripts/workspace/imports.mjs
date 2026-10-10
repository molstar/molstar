// AST inspection only; builds and declaration checks use the native TypeScript 7 CLI.
import ts from '@typescript/typescript6';

/**
 * Lists the import/export/require/dynamic-import specifiers of a source file.
 *
 * `typeOnly` flags imports whose bindings are all used in type positions. `erased` flags what the compiler removes
 * under `verbatimModuleSyntax`, which is only `import type` / `export type` declarations and import types: every other
 * import declaration stays in the output, including `import { type A }` (emitted as `import {} from '...'`) and a
 * binding used only as a type. A bundler therefore evaluates the module unless `erased` is true.
 */
export function importsFrom(file, source) {
  const result = [];
  function isTypeOnlyImport(node) {
    const clause = node.importClause;
    if (clause?.isTypeOnly) return true;
    const names = [];
    if (clause?.name) names.push(clause.name.text);
    if (clause?.namedBindings && ts.isNamespaceImport(clause.namedBindings)) names.push(clause.namedBindings.name.text);
    if (clause?.namedBindings && ts.isNamedImports(clause.namedBindings))
      names.push(...clause.namedBindings.elements.map((element) => element.name.text));
    if (!names.length) return false;
    const uses = new Map(names.map((name) => [name, []]));
    const findUses = (current) => {
      let declarationParent = current;
      while (declarationParent && declarationParent !== source && declarationParent !== node.importClause)
        declarationParent = declarationParent.parent;
      if (ts.isIdentifier(current) && uses.has(current.text) && declarationParent !== node.importClause) {
        let parent = current.parent,
          typePosition = false;
        while (parent && parent !== source) {
          if (ts.isTypeNode(parent)) {
            typePosition = true;
            break;
          }
          if (ts.isExpression(parent)) break;
          parent = parent.parent;
        }
        uses.get(current.text).push(typePosition);
      }
      ts.forEachChild(current, findUses);
    };
    findUses(source);
    const found = [...uses.values()].flat();
    return found.length > 0 && found.every(Boolean);
  }
  const visit = (node) => {
    if (ts.isImportDeclaration(node) && node.moduleSpecifier && ts.isStringLiteral(node.moduleSpecifier))
      result.push({
        specifier: node.moduleSpecifier.text,
        typeOnly: isTypeOnlyImport(node),
        erased: !!node.importClause?.isTypeOnly,
      });
    if (ts.isExportDeclaration(node) && node.moduleSpecifier && ts.isStringLiteral(node.moduleSpecifier))
      result.push({ specifier: node.moduleSpecifier.text, typeOnly: node.isTypeOnly, erased: node.isTypeOnly });
    if (
      ts.isImportEqualsDeclaration(node) &&
      ts.isExternalModuleReference(node.moduleReference) &&
      node.moduleReference.expression &&
      ts.isStringLiteral(node.moduleReference.expression)
    )
      result.push({ specifier: node.moduleReference.expression.text, typeOnly: false, erased: false });
    if (
      ts.isCallExpression(node) &&
      node.expression.kind === ts.SyntaxKind.ImportKeyword &&
      node.arguments.length === 1 &&
      ts.isStringLiteral(node.arguments[0])
    )
      result.push({ specifier: node.arguments[0].text, typeOnly: false, erased: false });
    if (
      ts.isCallExpression(node) &&
      ts.isIdentifier(node.expression) &&
      node.expression.text === 'require' &&
      node.arguments.length === 1 &&
      ts.isStringLiteral(node.arguments[0])
    )
      result.push({ specifier: node.arguments[0].text, typeOnly: false, erased: false });
    if (ts.isImportTypeNode(node) && ts.isLiteralTypeNode(node.argument) && ts.isStringLiteral(node.argument.literal))
      result.push({ specifier: node.argument.literal.text, typeOnly: true, erased: true });
    ts.forEachChild(node, visit);
  };
  visit(source);
  return result;
}
