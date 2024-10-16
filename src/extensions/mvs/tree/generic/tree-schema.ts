/**
 * Copyright (c) 2023-2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */

import { onelinerJsonString } from '../../../../mol-util/json';
import { isPlainObject, mapObjectMap } from '../../../../mol-util/object';
import { AllRequired, DefaultsFor, ParamsSchema, ValuesFor, paramsValidationIssues } from './params-schema';
import { treeToString } from './tree-utils';


/** Type of "custom" of a tree node (key-value storage with arbitrary JSONable values) */
export type CustomProps = Partial<Record<string, any>>

/** Tree node without children */
export type Node<TKind extends string = string, TParams extends {} = {}> =
    {} extends TParams ? {
        kind: TKind,
        params?: TParams, // params can be dropped if {} is valid value for params
        custom?: CustomProps,
        ref?: string,
    } : {
        kind: TKind,
        params: TParams, // params must be here if {} is not valid value for params
        custom?: CustomProps,
        ref?: string,
    }

/** Kind type for a tree node */
export type Kind<TNode extends Node> = TNode['kind']

/** Params type for a tree node */
export type Params<TNode extends Node> = NonNullable<TNode['params']>


/** Tree (i.e. a node with optional children) where the root node is of type `TRoot` and other nodes are of type `TNode` */
export type Tree<TNode extends Node<string, {}> = Node<string, {}>, TRoot extends TNode = TNode> =
    TRoot & {
        children?: Tree<TNode, TNode>[],
    }

/** Type of any subtree that can occur within given `TTree` tree type */
export type Subtree<TTree extends Tree> = NonNullable<TTree['children']>[number]

/** Type of any subtree that can occur within given `TTree` tree type and has kind type `TKind` */
export type SubtreeOfKind<TTree extends Tree, TKind extends Kind<Subtree<TTree>> = Kind<Subtree<TTree>>> = RootOfKind<Subtree<TTree>, TKind>

type RootOfKind<TTree extends Tree, TKind extends Kind<TTree>> = Extract<TTree, Tree<any, Node<TKind>>>

/** Params type for a given kind type within a tree */
export type ParamsOfKind<TTree extends Tree, TKind extends Kind<Subtree<TTree>> = Kind<Subtree<TTree>>> = NonNullable<SubtreeOfKind<TTree, TKind>['params']>


/** Get params from a tree node */
export function getParams<TNode extends Node>(node: TNode): Params<TNode> {
    return node.params ?? {};
}
/** Get custom properties from a tree node */
export function getCustomProps<TCustomProps extends CustomProps = CustomProps>(node: Node): TCustomProps {
    return (node.custom ?? {}) as TCustomProps;
}
/** Get children from a tree node */
export function getChildren<TTree extends Tree>(tree: TTree): Subtree<TTree>[] {
    return tree.children ?? [];
}


type ParamsSchemas = { [kind: string]: ParamsSchema }

/** Definition of tree type, specifying allowed node kinds, types of their params, required kind for the root, and allowed parent-child kind combinations */
export interface TreeSchema<TParamsSchemas extends ParamsSchemas = ParamsSchemas, TRootKind extends keyof TParamsSchemas = string> {
    /** Required kind of the root node */
    rootKind: TRootKind,
    /** Definition of allowed node kinds */
    nodes: {
        [kind in keyof TParamsSchemas]: {
            /** Params schema for this node kind */
            params: TParamsSchemas[kind],
            /** Documentation for this node kind */
            description?: string,
            /** Node kinds that can serve as parent for this node kind (`undefined` means the parent can be of any kind) */
            parent?: (string & keyof TParamsSchemas)[],
        }
    },
}
export function TreeSchema<P extends ParamsSchemas = ParamsSchemas, R extends keyof P = string>(schema: TreeSchema<P, R>): TreeSchema<P, R> {
    return schema as any;
}

/** ParamsSchemas per node kind */
type ParamsSchemasOf<TTreeSchema extends TreeSchema> = TTreeSchema extends TreeSchema<infer TParamsSchema, any> ? TParamsSchema : never;

/** Variation of params schemas where all param fields are required */
type ParamsSchemasWithAllRequired<TParamsSchemas extends ParamsSchemas> = { [kind in keyof TParamsSchemas]: AllRequired<TParamsSchemas[kind]> }

/** Variation of a tree schema where all param fields are required */
export type TreeSchemaWithAllRequired<TTreeSchema extends TreeSchema> = TreeSchema<ParamsSchemasWithAllRequired<ParamsSchemasOf<TTreeSchema>>, TTreeSchema['rootKind']>
export function TreeSchemaWithAllRequired<TTreeSchema extends TreeSchema>(schema: TTreeSchema): TreeSchemaWithAllRequired<TTreeSchema> {
    return {
        ...schema,
        nodes: mapObjectMap(schema.nodes, node => ({ ...node, params: AllRequired(node.params) })) as any,
    };
}

/** Type of tree node which can occur as the root of a tree conforming to tree schema `TTreeSchema` */
export type RootFor<TTreeSchema extends TreeSchema> = NodeFor<TTreeSchema, TTreeSchema['rootKind']>

/** Type of tree node which can occur anywhere in a tree conforming to tree schema `TTreeSchema`,
 * optionally narrowing down to a given node kind */
export type NodeFor<TTreeSchema extends TreeSchema, TKind extends keyof ParamsSchemasOf<TTreeSchema> = keyof ParamsSchemasOf<TTreeSchema>>
    = { [key in keyof ParamsSchemasOf<TTreeSchema>]: Node<key & string, ValuesFor<ParamsSchemasOf<TTreeSchema>[key]>> }[TKind]

/** Type of tree which conforms to tree schema `TTreeSchema` */
export type TreeFor<TTreeSchema extends TreeSchema> = Tree<NodeFor<TTreeSchema>, RootFor<TTreeSchema> & NodeFor<TTreeSchema>>

/** Type of default parameter values for each node kind in a tree schema `TTreeSchema` */
export type DefaultsForTree<TTreeSchema extends TreeSchema> = { [kind in keyof TTreeSchema['nodes']]: DefaultsFor<TTreeSchema['nodes'][kind]['params']> }


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
export function validateTree(schema: TreeSchema, tree: Tree, label: string): void {
    const issues = treeValidationIssues(schema, tree, { noExtra: true });
    if (issues) {
        console.warn(`Invalid ${label} tree:\n${treeToString(tree)}`);
        console.error(`${label} tree validation issues:`);
        for (const line of issues) {
            console.error(' ', line);
        }
        throw new Error('FormatError');
    }
}

/** Return documentation for a tree schema as plain text */
export function treeSchemaToString<S extends TreeSchema>(schema: S, defaults?: DefaultsForTree<S>): string {
    return treeSchemaToString_(schema, defaults, false);
}
/** Return documentation for a tree schema as markdown text */
export function treeSchemaToMarkdown<S extends TreeSchema>(schema: S, defaults?: DefaultsForTree<S>): string {
    return treeSchemaToString_(schema, defaults, true);
}
function treeSchemaToString_<S extends TreeSchema>(schema: S, defaults?: DefaultsForTree<S>, markdown: boolean = false): string {
    const out: string[] = [];
    const bold = (str: string) => markdown ? `**${str}**` : str;
    const code = (str: string) => markdown ? `\`${str}\`` : str;
    const h1 = markdown ? '## ' : '  - ';
    const p1 = markdown ? '' : '    ';
    const h2 = markdown ? '- ' : '      - ';
    const p2 = markdown ? '  ' : '        ';
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
        for (const key in params) {
            const field = params[key];
            let typeString = field.type.name;
            if (typeString.startsWith('(') && typeString.endsWith(')')) {
                typeString = typeString.slice(1, -1);
            }
            out.push(`${h2}${bold(code(key + (field.required ? ': ' : '?: ')))}${code(typeString)}`);
            const defaultValue = (defaults?.[kind] as any)?.[key];
            if (field.description) {
                out.push(`${p2}${field.description}`);
            }
            if (defaultValue !== undefined) {
                out.push(`${p2}Default: ${code(onelinerJsonString(defaultValue))}`);
            }
        }
    }
    return out.join(newline);
}
