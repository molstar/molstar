/**
 * Copyright (c) 2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */

import { Structure, StructureSelection } from '../../../mol-model/structure';
import { StructureQueryHelper } from '../../../mol-plugin-state/helpers/structure-query';
import { PluginStateObject as SO } from '../../../mol-plugin-state/objects';
import { StateObject, StateTransformer } from '../../../mol-state';
import { deepEqual } from '../../../mol-util';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { rowsToExpression } from '../helpers/selections';
import { omitObjectKeys } from '../helpers/utils';
import { getAnnotationForStructure } from './annotation-prop';


/** Parameter definition for `AnnotationStructureComponent` transformer */
export const AnnotationStructureComponentParams = {
    annotationId: PD.Text('', { description: 'Reference to "Annotation" custom model property' }),
    fieldName: PD.Text('component', { description: 'Annotation field (column) from which to take component identifier' }),
    fieldValues: PD.MappedStatic('all', {
        all: PD.EmptyGroup(),
        selected: PD.ObjectList({
            value: PD.Text(),
        }, obj => obj.value),
    }),
    nullIfEmpty: PD.Optional(PD.Boolean(true, { isHidden: false })),
    label: PD.Text('', { isHidden: false }),
};

/** Parameter values for `AnnotationStructureComponent` transformer */
export type AnnotationStructureComponentProps = PD.ValuesFor<typeof AnnotationStructureComponentParams>


/** Transformer builder for MVS extension */
export const MVSTransform = StateTransformer.builderFactory('mvs');

/** Transformer for creating a structure component based on custom model property "Annotations" */
export type AnnotationStructureComponent = typeof AnnotationStructureComponent
export const AnnotationStructureComponent = MVSTransform({
    name: 'mvs-structure-component-from-annotation',
    display: { name: 'Annotation Component', description: 'A molecular structure component defined by annotation data.' },
    from: SO.Molecule.Structure,
    to: SO.Molecule.Structure,
    params: AnnotationStructureComponentParams,
})({
    apply({ a, params }) {
        return createAnnotationStructureComponent(a.data, params);
    },
    update: ({ a, b, oldParams, newParams }) => {
        return updateAnnotationStructureComponent(a.data, b, oldParams, newParams);
    },
    dispose({ b }) {
        b?.data.customPropertyDescriptors.dispose();
    }
});


/** Create a substructure based on `AnnotationStructureComponentProps` */
export function createAnnotationSubstructure(structure: Structure, params: AnnotationStructureComponentProps): Structure {
    const { annotation } = getAnnotationForStructure(structure, params.annotationId);
    if (annotation) {
        let rows = annotation.getRows();
        if (params.fieldValues.name === 'selected') {
            const selectedValues = new Set<string | undefined>(params.fieldValues.params.map(obj => obj.value));
            rows = rows.filter((row, i) => selectedValues.has(annotation.getValueForRow(i, params.fieldName)));
        }
        const expression = rowsToExpression(rows);

        const { selection } = StructureQueryHelper.createAndRun(structure, expression);
        return StructureSelection.unionStructure(selection);
    } else {
        return Structure.Empty;
    }
}

/** Create a substructure PSO based on `AnnotationStructureComponentProps` */
export function createAnnotationStructureComponent(structure: Structure, params: AnnotationStructureComponentProps) {
    const component = createAnnotationSubstructure(structure, params);

    if (params.nullIfEmpty && component.elementCount === 0) return StateObject.Null;

    let label = params.label;
    if (label === undefined || label === '') {
        if (params.fieldValues.name === 'selected' && params.fieldValues.params.length > 0) {
            const ellipsis = params.fieldValues.params.length > 1 ? '+...' : '';
            label = `${params.fieldName}: "${params.fieldValues.params[0].value}"${ellipsis}`;
        } else {
            label = 'Component from Annotation';
        }
    }

    const props = { label, description: Structure.elementDescription(component) };
    return new SO.Molecule.Structure(component, props);
}

/** Update a substructure PSO based on `AnnotationStructureComponentProps` */
export function updateAnnotationStructureComponent(a: Structure, b: SO.Molecule.Structure, oldParams: AnnotationStructureComponentProps, newParams: AnnotationStructureComponentProps) {
    const change = !deepEqual(newParams, oldParams);
    const needsRecreate = !deepEqual(omitObjectKeys(newParams, ['label']), omitObjectKeys(oldParams, ['label']));
    if (!change) {
        return StateTransformer.UpdateResult.Unchanged;
    }
    if (!needsRecreate) {
        b.label = newParams.label || b.label;
        return StateTransformer.UpdateResult.Updated;
    }
    return StateTransformer.UpdateResult.Recreate;
}
