/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { RuntimeContext } from 'mol-task';
import { Structure, StructureElement } from '../structure';

export interface QueryContextView {
    readonly element: StructureElement;
}

export class QueryContext implements QueryContextView {
    private currentStack: StructureElement[] = [];

    readonly structure: Structure;
    readonly taskCtx: RuntimeContext;

    /** Current element */
    readonly element: StructureElement = StructureElement.create();

    pushCurrentElement(): StructureElement {
        this.currentStack[this.currentStack.length] = this.element;
        (this.element as StructureElement) = StructureElement.create();
        return this.element;
    }

    popCurrentElement() {
        (this.element as StructureElement) = this.currentStack.pop()!;
    }

    constructor(structure: Structure, taskCtx: RuntimeContext) {
        this.structure = structure;
        this.taskCtx = taskCtx;
    }
}

export interface QueryPredicate { (ctx: QueryContextView): boolean }
export interface QueryFn<T = any> { (ctx: QueryContextView): T }