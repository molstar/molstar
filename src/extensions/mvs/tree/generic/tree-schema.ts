/**
 * Copyright (c) 2023-2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */

import { mapObjectMap } from '../../../../mol-util/object';
import { AllRequired, ParamsSchema, ValuesFor } from './params-schema';

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
    return schema;
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
