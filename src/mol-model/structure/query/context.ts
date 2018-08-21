/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Structure, StructureElement } from '../structure';
import { now } from 'mol-task';

export interface QueryContextView {
    readonly element: StructureElement;
    readonly currentStructure: Structure;
}

export class QueryContext implements QueryContextView {
    private currentElementStack: StructureElement[] = [];
    private currentStructureStack: Structure[] = [];
    private inputStructureStack: Structure[] = [];

    private timeCreated = now();

    readonly timeoutMs: number;
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

    pushInputStructure(structure: Structure) {
        this.inputStructureStack.push(this.inputStructure);
        (this.inputStructure as any) = structure;
    }

    popInputStructure() {
        if (this.inputStructureStack.length === 0) throw new Error('Must push before pop.');
        (this.inputStructure as any) = this.inputStructureStack.pop();
    }

    throwIfTimedOut() {
        if (this.timeoutMs === 0) return;
        if (now() - this.timeCreated > this.timeoutMs) {
            throw new Error(`The query took too long to execute (> ${this.timeoutMs / 1000}s).`);
        }
    }

    constructor(structure: Structure, timeoutMs = 0) {
        this.inputStructure = structure;
        this.timeoutMs = timeoutMs;
    }
}

export interface QueryPredicate { (ctx: QueryContext): boolean }
export interface QueryFn<T = any> { (ctx: QueryContext): T }