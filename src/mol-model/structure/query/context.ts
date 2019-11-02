/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Structure, StructureElement, Unit } from '../structure';
import { now } from '../../../mol-util/now';
import { ElementIndex } from '../model';
import { LinkType } from '../model/types';
import { StructureSelection } from './selection';
import { defaultLinkTest } from './queries/internal';

export interface QueryContextView {
    readonly element: StructureElement.Location;
    readonly currentStructure: Structure;
}

export class QueryContext implements QueryContextView {
    private currentElementStack: StructureElement.Location[] = [];
    private currentAtomicLinkStack: QueryContextLinkInfo<Unit.Atomic>[] = [];
    private currentStructureStack: Structure[] = [];
    private inputStructureStack: Structure[] = [];

    private timeCreated = now();

    readonly timeoutMs: number;
    readonly inputStructure: Structure;

    /** Current element */
    readonly element = StructureElement.Location.create();
    currentStructure: Structure = void 0 as any;

    /** Current link between atoms */
    readonly atomicLink = new QueryContextLinkInfo<Unit.Atomic>();

    /** Supply this from the outside. Used by the internal.generator.current symbol */
    currentSelection: StructureSelection | undefined = void 0;

    setElement(unit: Unit, e: ElementIndex) {
        this.element.unit = unit;
        this.element.element = e;
    }

    pushCurrentElement(): StructureElement.Location {
        this.currentElementStack[this.currentElementStack.length] = this.element;
        (this.element as StructureElement.Location) = StructureElement.Location.create();
        return this.element;
    }

    popCurrentElement() {
        (this.element as StructureElement.Location) = this.currentElementStack.pop()!;
    }

    pushCurrentLink() {
        if (this.atomicLink) this.currentAtomicLinkStack.push(this.atomicLink);
        (this.atomicLink as QueryContextLinkInfo<Unit.Atomic>) = new QueryContextLinkInfo();
        return this.atomicLink;
    }

    popCurrentLink() {
        if (this.currentAtomicLinkStack.length > 0) {
            (this.atomicLink as QueryContextLinkInfo<Unit.Atomic>) = this.currentAtomicLinkStack.pop()!;
        } else {
            (this.atomicLink as any) = void 0;
        }
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

    tryGetCurrentSelection() {
        if (!this.currentSelection) throw new Error('The current selection is not assigned.');
        return this.currentSelection;
    }

    constructor(structure: Structure, options?: QueryContextOptions) {
        this.inputStructure = structure;
        this.timeoutMs = (options && options.timeoutMs) || 0;
        this.currentSelection = options && options.currentSelection;
    }
}

export interface QueryContextOptions {
    timeoutMs?: number,
    currentSelection?: StructureSelection
}

export interface QueryPredicate { (ctx: QueryContext): boolean }
export interface QueryFn<T = any> { (ctx: QueryContext): T }

export class QueryContextLinkInfo<U extends Unit = Unit> {
    a: StructureElement.Location<U> = StructureElement.Location.create();
    aIndex: StructureElement.UnitIndex = 0 as StructureElement.UnitIndex;
    b: StructureElement.Location<U> = StructureElement.Location.create();
    bIndex: StructureElement.UnitIndex = 0 as StructureElement.UnitIndex;
    type: LinkType = LinkType.Flag.None;
    order: number = 0;

    private testFn: QueryPredicate = defaultLinkTest;

    setTestFn(fn?: QueryPredicate) {
        this.testFn = fn || defaultLinkTest;
    }

    test(ctx: QueryContext, trySwap: boolean) {
        if (this.testFn(ctx)) return true;
        if (trySwap) {
            this.swap();
            return this.testFn(ctx);
        }
        return false;
    }

    private swap() {
        const idxA = this.aIndex;
        this.aIndex = this.bIndex;
        this.bIndex = idxA;

        const unitA = this.a.unit;
        this.a.unit = this.b.unit;
        this.b.unit = unitA;

        const eA = this.a.element;
        this.a.element = this.b.element;
        this.b.element = eA;
    }

    get length() {
        return StructureElement.Location.distance(this.a, this. b);
    }
}