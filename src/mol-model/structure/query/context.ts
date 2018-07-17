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
    private currentStack: StructureElement[] = [];
    private currentLocked = false;

    readonly inputStructure: Structure;

    /** Current element */
    readonly element: StructureElement = StructureElement.create();
    readonly currentStructure: Structure = void 0 as any;

    pushCurrentElement(): StructureElement {
        this.currentStack[this.currentStack.length] = this.element;
        (this.element as StructureElement) = StructureElement.create();
        return this.element;
    }

    popCurrentElement() {
        (this.element as StructureElement) = this.currentStack.pop()!;
    }

    lockCurrentStructure(structure: Structure) {
        if (this.currentLocked) throw new Error('Current structure already locked.');
        this.currentLocked = true;
        (this.currentStructure as Structure) = structure;
    }

    unlockCurrentStructure() {
        this.currentLocked = false;
        (this.currentStructure as any) = void 0;
    }

    constructor(structure: Structure) {
        this.inputStructure = structure;
    }
}

export interface QueryPredicate { (ctx: QueryContextView): boolean }
export interface QueryFn<T = any> { (ctx: QueryContextView): T }