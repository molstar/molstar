/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { PluginContext } from 'mol-plugin/context';
import { PluginBehavior } from '../behavior';
import { ParamDefinition as PD } from 'mol-util/param-definition'
import { Mat4, Vec3 } from 'mol-math/linear-algebra';
import { PluginStateObject as SO, PluginStateObject } from '../../state/objects';
import { StateObjectCell, State, StateSelection } from 'mol-state';
import { RuntimeContext } from 'mol-task';
import { Shape } from 'mol-model/shape';
import { Text } from 'mol-geo/geometry/text/text';
import { ShapeRepresentation } from 'mol-repr/shape/representation';
import { ColorNames } from 'mol-util/color/tables';
import { TextBuilder } from 'mol-geo/geometry/text/text-builder';
import { Unit, StructureElement, StructureProperties } from 'mol-model/structure';
import { SetUtils } from 'mol-util/set';
import { arrayEqual } from 'mol-util';
import { MoleculeType } from 'mol-model/structure/model/types';
import { getElementMoleculeType } from 'mol-model/structure/util';

// TODO
// - support more object types than structures
// - tether label to the element nearest to the bounding sphere center
// - [Started] multiple levels of labels: structure, polymer, ligand
// - show structure/unit label only when there is a representation with sufficient overlap
// - support highlighting
// - better support saccharides (use data available after re-mediation)
// - size based on min bbox dimension (to avoid huge labels for very long but narrow polymers)
// - fixed size labels (invariant to zoom) [needs feature in text geo]
// - ??? max label length
// - ??? multi line labels [needs feature in text geo]
// - ??? use prevalent (how to define) color of representations of a structure to color the label
// - completely different approach (render not as 3d objects): overlay free layout in screenspace with occlusion info from bboxes

export type SceneLabelsLevels = 'structure' | 'polymer' | 'ligand'

export const SceneLabelsParams = {
    ...Text.Params,

    background: PD.Boolean(true),
    backgroundMargin: PD.Numeric(0.2, { min: 0, max: 1, step: 0.01 }),
    backgroundColor: PD.Color(ColorNames.snow),
    backgroundOpacity: PD.Numeric(0.9, { min: 0, max: 1, step: 0.01 }),

    levels: PD.MultiSelect([] as SceneLabelsLevels[], [
        ['structure', 'structure'], ['polymer', 'polymer'], ['ligand', 'ligand']
    ] as [SceneLabelsLevels, string][]),
}
export type SceneLabelsParams = typeof SceneLabelsParams
export type SceneLabelsProps = PD.Values<typeof SceneLabelsParams>

interface LabelsData {
    transforms: Mat4[]
    texts: string[]
    positions: Vec3[]
    sizes: number[]
    depths: number[]
}

function getLabelsText(data: LabelsData, props: PD.Values<Text.Params>, text?: Text) {
    const { texts, positions, depths } = data
    const textBuilder = TextBuilder.create(props, texts.length * 10, texts.length * 10 / 2, text)
    for (let i = 0, il = texts.length; i < il; ++i) {
        const p = positions[i]
        textBuilder.add(texts[i], p[0], p[1], p[2], depths[i], i)
    }
    return textBuilder.getText()
}

export const SceneLabels = PluginBehavior.create<SceneLabelsProps>({
    name: 'scene-labels',
    category: 'representation',
    display: { name: 'Scene Labels' },
    canAutoUpdate: () => true,
    ctor: class extends PluginBehavior.Handler<SceneLabelsProps> {
        private data: LabelsData = {
            transforms: [Mat4.identity()],
            texts: [],
            positions: [],
            sizes: [],
            depths: []
        }
        private repr: ShapeRepresentation<LabelsData, Text, SceneLabelsParams>
        private geo = Text.createEmpty()
        private structures = new Set<SO.Molecule.Structure>()

        constructor(protected ctx: PluginContext, protected params: SceneLabelsProps) {
            super(ctx, params)
            this.repr = ShapeRepresentation(this.getLabelsShape, Text.Utils)
            ctx.events.state.object.created.subscribe(this.triggerUpdate)
            ctx.events.state.object.removed.subscribe(this.triggerUpdate)
            ctx.events.state.object.updated.subscribe(this.triggerUpdate)
            ctx.events.state.cell.stateUpdated.subscribe(this.triggerUpdate)
        }

        private triggerUpdate = async () => {
            await this.update(this.params)
        }

        private getColor = () => ColorNames.dimgrey
        private getSize = (groupId: number) => this.data.sizes[groupId]
        private getLabel = () => ''

        private getLabelsShape = (ctx: RuntimeContext, data: LabelsData, props: SceneLabelsProps, shape?: Shape<Text>) => {
            this.geo = getLabelsText(data, props, this.geo)
            return Shape.create('Scene Labels', data, this.geo, this.getColor, this.getSize, this.getLabel, data.transforms)
        }

        /** Update structures to be labeled, returns true if changed */
        private updateStructures(p: SceneLabelsProps) {
            const state = this.ctx.state.dataState
            const structures = state.selectQ(q => q.rootsOfType(PluginStateObject.Molecule.Structure));
            const rootStructures = new Set<SO.Molecule.Structure>()
            for (const s of structures) {
                const rootStructure = getRootStructure(s, state)
                if (!rootStructure || !SO.Molecule.Structure.is(rootStructure.obj)) continue
                if (!state.cellStates.get(s.transform.ref).isHidden) {
                    rootStructures.add(rootStructure.obj)
                }
            }
            if (!SetUtils.areEqual(rootStructures, this.structures)) {
                this.structures = rootStructures
                return true
            } else {
                return false
            }
        }

        private updateLabels(p: SceneLabelsProps) {
            const l = StructureElement.create()

            const { texts, positions, sizes, depths } = this.data
            texts.length = 0
            positions.length = 0
            sizes.length = 0
            depths.length = 0

            this.structures.forEach(structure => {
                if (p.levels.includes('structure')) {
                    texts.push(`${structure.data.model.label}`)
                    positions.push(structure.data.boundary.sphere.center)
                    sizes.push(structure.data.boundary.sphere.radius / 10)
                    depths.push(structure.data.boundary.sphere.radius)
                }

                for (let i = 0, il = structure.data.units.length; i < il; ++i) {
                    let label = ''
                    const u = structure.data.units[i]
                    l.unit = u
                    l.element = u.elements[0]

                    if (p.levels.includes('polymer') && u.polymerElements.length) {
                        label = `${StructureProperties.entity.pdbx_description(l).join(', ')} (${getAsymId(u)(l)})`
                    }

                    if (p.levels.includes('ligand') && !u.polymerElements.length) {
                        const moleculeType = getElementMoleculeType(u, u.elements[0])
                        if (moleculeType === MoleculeType.other || moleculeType === MoleculeType.saccharide) {
                            label = `${StructureProperties.entity.pdbx_description(l).join(', ')} (${getAsymId(u)(l)})`
                        }
                    }

                    if (label) {
                        texts.push(label)
                        const { center, radius } = u.lookup3d.boundary.sphere
                        const transformedCenter = Vec3.transformMat4(Vec3.zero(), center, u.conformation.operator.matrix)
                        positions.push(transformedCenter)
                        sizes.push(Math.max(2, radius / 10))
                        depths.push(radius)
                    }
                }
            })
        }

        register(): void { }

        async update(p: SceneLabelsProps) {
            // console.log('update')
            let updated = false
            if (this.updateStructures(p) || !arrayEqual(this.params.levels, p.levels)) {
                // console.log('update with data')
                this.updateLabels(p)
                await this.repr.createOrUpdate(p, this.data).run()
                updated = true
            } else if (!PD.areEqual(SceneLabelsParams, this.params, p)) {
                // console.log('update props only')
                await this.repr.createOrUpdate(p).run()
                updated = true
            }
            if (updated) {
                Object.assign(this.params, p)
                this.ctx.canvas3d.add(this.repr)
            }
            return updated;
        }

        unregister() {

        }
    },
    params: () => SceneLabelsParams
});

//

function getRootStructure(root: StateObjectCell, state: State) {
    let parent: StateObjectCell | undefined
    while (true) {
        const _parent = StateSelection.findAncestorOfType(state.tree, state.cells, root.transform.ref, [PluginStateObject.Molecule.Structure])
        if (_parent) {
            parent = _parent
            root = _parent
        } else {
            break
        }
    }
    return parent ? parent :
        SO.Molecule.Structure.is(root.obj) ? root : undefined
}

function getAsymId(unit: Unit): StructureElement.Property<string> {
    switch (unit.kind) {
        case Unit.Kind.Atomic:
            return StructureProperties.chain.auth_asym_id
        case Unit.Kind.Spheres:
        case Unit.Kind.Gaussians:
            return StructureProperties.coarse.asym_id
    }
}