/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Structure, StructureElement } from '../structure';

export interface QueryContextView {
    readonly element: StructureElement;
    readonly currentStructure: Structure;
}

export class QueryContext implements QueryContextView {
    private currentElementStack: StructureElement[] = [];
    private currentStructureStack: Structure[] = [];

    readonly inputStructure: Structure;

    /** Current element */
    readonly element: StructureElement = StructureElement.create();
    currentStructure: Structure = void 0 as any;

    pushCurrentElement(): StructureElement {
        this.currentElementStack[this.currentElementStack.length] = this.element;
        (this.element as StructureElement) = StructureElement.create();
        return this.element;
    }

    popCurrentElement() {
        (this.element as StructureElement) = this.currentElementStack.pop()!;
    }

    pushCurrentStructure() {
        if (this.currentStructure) this.currentStructureStack.push(this.currentStructure);
    }

    popCurrentStructure() {
        if (this.currentStructureStack.length) (this.currentStructure as Structure) = this.currentStructureStack.pop()!;
        else (this.currentStructure as Structure) = void 0 as any;
    }

    constructor(structure: Structure) {
        this.inputStructure = structure;
    }
}

export interface QueryPredicate { (ctx: QueryContextView): boolean }
export interface QueryFn<T = any> { (ctx: QueryContextView): T }