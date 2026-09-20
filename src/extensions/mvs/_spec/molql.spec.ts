/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { componentPropsFromSelector, prettyNameFromSelector } from '../load-helpers';
import { MVSData } from '../mvs-data';
import { createMVSBuilder } from '../tree/mvs/mvs-builder';
import { MolQLExpressionT, PrimitiveMolQLExpressionT, isMolQLExpression, isPrimitiveMolQLExpression } from '../tree/mvs/param-types';
import { buildStory } from '../../../examples/mvs-stories/stories/molql';
import { compile } from '../../../mol-script/runtime/query/base';
import { MolScriptBuilder as MS } from '../../../mol-script/language/builder';
import { parse } from '../../../mol-script/transpile';

describe('MVS MolQL selectors', () => {
    const molql = MS.struct.generator.atomGroups({});

    it('requires a compilable molql application wrapper', () => {
        const malformed = { head: { name: 'not-implemented' } };
        expect(MolQLExpressionT.decode({ molql })._tag).toEqual('Right');
        expect(MolQLExpressionT.decode({ molql })._tag).toEqual('Right');
        expect(MolQLExpressionT.decode({ molql: 'not-a-query' })._tag).toEqual('Left');
        expect(MolQLExpressionT.decode({ molql: malformed })._tag).toEqual('Left');
        expect(MolQLExpressionT.decode({ molql, label_seq_id: 1 })._tag).toEqual('Left');
        expect(MolQLExpressionT.decode({})._tag).toEqual('Left');
        expect(isMolQLExpression({ molql })).toEqual(true);
        expect(isMolQLExpression({})).toEqual(false);
        expect(isMolQLExpression([])).toEqual(false);
    });

    it('allows primitive MolQL positions to reference another structure', () => {
        const position = { molql, structure_ref: 'other-structure' };
        expect(PrimitiveMolQLExpressionT.decode(position)._tag).toEqual('Right');
        expect(PrimitiveMolQLExpressionT.decode({ molql, structure_ref: 1 })._tag).toEqual('Left');
        expect(isPrimitiveMolQLExpression(position)).toEqual(true);
        expect(isMolQLExpression({ molql, structure_ref: 1 })).toEqual(false);
        expect(isPrimitiveMolQLExpression({ molql, structure_ref: 1 })).toEqual(false);
    });

    it('passes wrapped expressions through unchanged and preserves legacy selectors', () => {
        const wrapped = { molql };
        expect(componentPropsFromSelector(wrapped)).toEqual({ name: 'expression', params: molql });
        expect(componentPropsFromSelector('protein')).toEqual({ name: 'static', params: 'protein' });
        expect(componentPropsFromSelector({ label_asym_id: 'A' })).toMatchObject({ name: 'expression' });
        expect(componentPropsFromSelector([{ label_asym_id: 'A' }])).toMatchObject({ name: 'expression' });
        expect(prettyNameFromSelector(wrapped)).toEqual('MolQL Selection');
    });

    it('rejects invalid selector and primitive addresses during MVS validation', () => {
        const builder = createMVSBuilder();
        builder.download({ url: 'example.bcif' }).parse({ format: 'bcif' }).modelStructure()
            .component({ selector: { molql: 'not-a-query' } as any });
        expect(MVSData.validationIssues(builder.getState())?.join('\n')).toContain('MolQL expression must be a MolScript application');

        const mixedBuilder = createMVSBuilder();
        mixedBuilder.download({ url: 'example.bcif' }).parse({ format: 'bcif' }).modelStructure()
            .component({ selector: { molql, label_seq_id: 1 } as any });
        expect(MVSData.validationIssues(mixedBuilder.getState())?.join('\n')).toContain('Expected an object with only a molql property');

        const primitiveBuilder = createMVSBuilder();
        const structure = primitiveBuilder.download({ url: 'example.bcif' }).parse({ format: 'bcif' }).modelStructure();
        structure.primitives().label({ position: { molql, structure_ref: 1 } as any, text: 'Invalid reference' });
        expect(MVSData.validationIssues(primitiveBuilder.getState())?.join('\n')).toContain('Expected an object with only molql and an optional string structure_ref');

        const referencedBuilder = createMVSBuilder();
        const referenced = referencedBuilder.download({ url: 'other.bcif' }).parse({ format: 'bcif' }).modelStructure({ ref: 'struct1' });
        referenced.primitives().label({ position: { structure_ref: 'struct1', label_asym_id: 'A' }, text: 'Legacy reference' });
        expect(MVSData.validationIssues(referencedBuilder.getState())).toEqual(undefined);
    });

    it('builds a valid story from MolScriptBuilder and a PyMOL-transpiled expression', () => {
        const story = buildStory();
        expect(MVSData.validationIssues(story)).toEqual(undefined);
        expect(story.snapshots).toHaveLength(4);
        expect(findNodes(story, 'camera')).toHaveLength(0);

        const componentSelectors = findNodes(story, 'component').map(node => node.params.selector);
        expect(componentSelectors.some(selector => selector?.molql)).toEqual(true);

        const colorSelectors = findNodes(story, 'color').map(node => node.params.selector);
        expect(colorSelectors.filter(selector => selector?.molql)).toHaveLength(2);

        const primitivePositions = findNodes(story, 'primitive').flatMap(node => [node.params.position, node.params.start, node.params.end]);
        expect(primitivePositions.filter(position => position?.molql)).toHaveLength(2);

        const serialized = MVSData.toMVSJ(story);
        expect(serialized).toContain('"kind":"component"');
        expect(serialized).toContain('"kind":"color"');
        expect(serialized).toContain('"kind":"primitive"');
        expect(serialized).toContain('"selector":{"molql"');
        expect(() => compile(parse('pymol', 'byres polymer within 5 of resn STI'))).not.toThrow();
        expect(() => compile(MS.struct.generator.atomGroups({
            'chain-test': MS.core.rel.eq([MS.struct.atomProperty.macromolecular.label_asym_id(), 'G']),
            'atom-test': MS.core.rel.eq([MS.struct.atomProperty.macromolecular.label_atom_id(), 'N13']),
        }))).not.toThrow();
        expect(() => compile(parse('pymol', 'chain A and resi 315 and name OG1'))).not.toThrow();
    });
});

function findNodes(tree: any, kind: string): any[] {
    const nodes: any[] = [];
    const visit = (node: any) => {
        if (node.kind === kind) nodes.push(node);
        for (const child of node.children ?? []) visit(child);
    };

    for (const snapshot of tree.snapshots ?? []) visit(snapshot.root);
    return nodes;
}
