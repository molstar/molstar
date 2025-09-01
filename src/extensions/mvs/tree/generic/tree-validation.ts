import { PluginContext } from '../../../../mol-plugin/context';
import { onelinerJsonString } from '../../../../mol-util/json';
import { isPlainObject } from '../../../../mol-util/object';
import { Field } from './field-schema';
import { SimpleParamsSchema, paramsValidationIssues } from './params-schema';
import { getChildren, getParams, Tree, TreeSchema } from './tree-schema';
import { treeToString } from './tree-utils';

/** Return `undefined` if a tree conforms to the given schema,
 * return validation issues (as a list of lines) if it does not conform.
 * If `options.requireAll`, all parameters (including optional) must have a value provided.
 * If `options.noExtra` is true, presence of any extra parameters is treated as an issue.
 * If `options.anyRoot` is true, the kind of the root node is not enforced.
 */
export function treeValidationIssues(schema: TreeSchema, tree: Tree, options: { requireAll?: boolean, noExtra?: boolean, anyRoot?: boolean, parent?: string } = {}): string[] | undefined {
    if (!isPlainObject(tree)) return [`Node must be an object, not ${tree}`];
    if (!options.anyRoot && tree.kind !== schema.rootKind) return [`Invalid root node kind "${tree.kind}", root must be of kind "${schema.rootKind}"`];
    const nodeSchema = schema.nodes[tree.kind];
    if (!nodeSchema) return [`Unknown node kind "${tree.kind}"`];
    if (nodeSchema.parent && (options.parent !== undefined) && !nodeSchema.parent.includes(options.parent)) {
        return [`Node of kind "${tree.kind}" cannot appear as a child of "${options.parent}". Allowed parents for "${tree.kind}" are: ${nodeSchema.parent.map(s => `"${s}"`).join(', ')}`];
    }
    const issues = paramsValidationIssues(nodeSchema.params, getParams(tree), options);
    if (issues) return [`Invalid parameters for node of kind "${tree.kind}":`, ...issues.map(s => '  ' + s)];
    if (tree.custom !== undefined && (typeof tree.custom !== 'object' || tree.custom === null)) {
        return [`Invalid "custom" for node of kind "${tree.kind}": must be an object, not ${tree.custom}.`];
    }
    for (const child of getChildren(tree)) {
        const issues = treeValidationIssues(schema, child, { ...options, anyRoot: true, parent: tree.kind });
        if (issues) return issues;
    }
    return undefined;
}

/** Validate a tree against the given schema.
 * Do nothing if OK; print validation issues on console and throw an error is the tree does not conform.
 * Include `label` in the printed output. */
export function validateTree(schema: TreeSchema, tree: Tree, label: string, plugin: PluginContext): void {
    const issues = treeValidationIssues(schema, tree, { noExtra: true });
    if (issues) {
        console.warn(`Invalid ${label} tree:\n${treeToString(tree)}`);
        console.error(`${label} tree validation issues:`);
        plugin.log.error(`${label} tree validation issues:`);
        for (const line of issues) {
            console.error(' ', line);
            plugin.log.error(line);
        }
        throw new Error('FormatError');
    }
}

/** Return documentation for a tree schema as plain text */
export function treeSchemaToString<S extends TreeSchema>(schema: S): string {
    return treeSchemaToString_(schema, false);
}
/** Return documentation for a tree schema as markdown text */
export function treeSchemaToMarkdown<S extends TreeSchema>(schema: S): string {
    return treeSchemaToString_(schema, true);
}
function treeSchemaToString_<S extends TreeSchema>(schema: S, markdown: boolean = false): string {
    const out: string[] = [];
    const bold = (str: string) => markdown ? `**${str}**` : str;
    const code = (str: string) => markdown ? `\`${str}\`` : str;
    const h1 = markdown ? '## ' : '  - ';
    const p1 = markdown ? '' : '    ';
    const h2 = markdown ? '- ' : '      - ';
    const p2 = markdown ? '  ' : '        ';
    const h3 = markdown ? '  - ' : '          - ';
    const p3 = markdown ? '    ' : '            ';
    const newline = markdown ? '\n\n' : '\n';
    out.push(`Tree schema:`);
    for (const kind in schema.nodes) {
        const { description, params, parent } = schema.nodes[kind];
        out.push(`${h1}${code(kind)}`);
        if (kind === schema.rootKind) {
            out.push(`${p1}[Root of the tree must be of this kind]`);
        }
        if (description) {
            out.push(`${p1}${description}`);
        }
        out.push(`${p1}Parent: ${!parent ? 'any' : parent.length === 0 ? 'none' : parent.map(code).join(' or ')}`);
        out.push(`${p1}Params:${Object.keys(params).length > 0 ? '' : ' none'}`);
        if (params.type === 'simple') {
            formatSimpleParams(out, params, { h: h2, p: p2, code, bold });
        } else {
            const key = params.discriminator;
            const casesStr = Object.keys(params.cases).join(' | ');
            out.push(`${h2}${bold(code(key + ': '))}${code(casesStr)}`);
            if (params.discriminatorDescription) {
                out.push(`${p2}${params.discriminatorDescription}`);
            }
            out.push(`${p2}[This parameter determines the rest of parameters]`);
            for (const case_ in params.cases) {
                const caseStr = `${params.discriminator}: "${case_}"`;
                out.push(`${p2}${bold(`Case ${code(caseStr)}:`)}`);
                formatSimpleParams(out, params.cases[case_], { h: h3, p: p3, code, bold });
            }
        }
    }
    return out.join(newline);
}

function formatSimpleParams(out: string[], params: SimpleParamsSchema, formatting: { h: string, p: string, code: (str: string) => string, bold: (str: string) => string }): void {
    const { h, p, code, bold } = formatting;
    for (const key in params.fields) {
        const field = params.fields[key];
        out.push(`${h}${bold(code(key + (field.required ? ': ' : '?: ')))}${code(formatFieldType(field))}`);
        const defaultValue = field.required ? undefined : field.default;
        if (field.description) {
            out.push(`${p}${field.description}`);
        }
        if (defaultValue !== undefined) {
            out.push(`${p}Default: ${code(onelinerJsonString(defaultValue))}`);
        }
    }
}

function formatFieldType(field: Field): string {
    const typeString = field.type.name;
    if (typeString.startsWith('(') && typeString.endsWith(')')) {
        return typeString.slice(1, -1);
    } else {
        return typeString;
    }
}
