/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { PluginContext } from 'mol-plugin/context';
import { PluginBehavior } from '../behavior';
import { ParamDefinition as PD } from 'mol-util/param-definition'
import { degToRad } from 'mol-math/misc';
import { Mat4, Vec3 } from 'mol-math/linear-algebra';
import { PluginStateObject as SO, PluginStateObject } from '../../state/objects';
import { StateSelection } from 'mol-state/state/selection';
import { StateObjectCell, State } from 'mol-state';
import { StructureUnitTransforms } from 'mol-model/structure/structure/util/unit-transforms';
import { UUID } from 'mol-util';

const StructureAnimationParams = {
    rotate: PD.Boolean(false),
    rotateValue: PD.Numeric(0, { min: 0, max: 360, step: 0.1 }),
    explode: PD.Boolean(false),
    explodeValue: PD.Numeric(0, { min: 0, max: 100, step: 0.1 }),
}
type StructureAnimationProps = PD.Values<typeof StructureAnimationParams>

/**
 * TODO
 * - animation class is just for testing purposes, needs better API
 * - allow per-unit transform `unitTransform: { [unitId: number]: Mat4 }`
 */
export const StructureAnimation = PluginBehavior.create<StructureAnimationProps>({
    name: 'structure-animation',
    display: { name: 'Structure Animation', group: 'Animation' },
    ctor: class extends PluginBehavior.Handler<StructureAnimationProps> {
        private tmpMat = Mat4.identity()
        private rotMat = Mat4.identity()
        private transMat = Mat4.identity()
        private animMat = Mat4.identity()
        private transVec = Vec3.zero()
        private rotVec = Vec3.create(0, 1, 0)
        private centerVec = Vec3.zero()

        private rotateAnimHandle = -1
        private explodeAnimHandle = -1

        private updatedUnitTransforms = new Set<SO.Molecule.Structure>()
        private structureUnitTransforms = new Map<UUID, StructureUnitTransforms>()

        constructor(protected ctx: PluginContext, protected params: StructureAnimationProps) {
            super(ctx, params)
            this.update(params)
        }

        private getUnitTransforms(structure: SO.Molecule.Structure) {
            let unitTransforms = this.structureUnitTransforms.get(structure.id)
            if (!unitTransforms) {
                unitTransforms = new StructureUnitTransforms(structure.data)
                this.structureUnitTransforms.set(structure.id, unitTransforms)
            }
            return unitTransforms
        }

        rotate(rad: number) {
            this.updatedUnitTransforms.clear()
            const state = this.ctx.state.dataState
            const reprs = state.select(q => q.rootsOfType(PluginStateObject.Molecule.Representation3D))
            Mat4.rotate(this.rotMat, this.tmpMat, rad, this.rotVec)
            for (const r of reprs) {
                if (!SO.isRepresentation3D(r.obj)) return
                const structure = getRootStructure(r, state)
                if (!structure || !SO.Molecule.Structure.is(structure.obj)) continue

                const unitTransforms = this.getUnitTransforms(structure.obj)

                if (!this.updatedUnitTransforms.has(structure.obj)) {
                    for (let i = 0, il = structure.obj.data.units.length; i < il; ++i) {
                        const u = structure.obj.data.units[i]
                        Vec3.transformMat4(this.centerVec, u.lookup3d.boundary.sphere.center, u.conformation.operator.matrix)

                        Vec3.negate(this.transVec, Vec3.copy(this.transVec, this.centerVec))
                        Mat4.fromTranslation(this.transMat, this.transVec)
                        Mat4.mul(this.animMat, this.rotMat, this.transMat)

                        Vec3.copy(this.transVec, this.centerVec)
                        Mat4.fromTranslation(this.transMat, this.transVec)
                        Mat4.mul(this.animMat, this.transMat, this.animMat)

                        unitTransforms.setTransform(this.animMat, u)
                    }
                    this.updatedUnitTransforms.add(structure.obj)
                }

                r.obj.data.setState({ unitTransforms })
                this.ctx.canvas3d.add(r.obj.data)
            }
            this.ctx.canvas3d.requestDraw(true)
        }

        animateRotate(play: boolean) {
            if (play) {
                const animateRotate = (t: number) => {
                    this.rotate(degToRad((t / 10) % 360))
                    this.rotateAnimHandle = requestAnimationFrame(animateRotate)
                }
                this.rotateAnimHandle = requestAnimationFrame(animateRotate)
            } else {
                cancelAnimationFrame(this.rotateAnimHandle)
            }
        }

        explode(p: number) {
            this.updatedUnitTransforms.clear()
            const state = this.ctx.state.dataState
            const reprs = state.select(q => q.rootsOfType(PluginStateObject.Molecule.Representation3D));
            for (const r of reprs) {
                if (!SO.isRepresentation3D(r.obj)) return
                const structure = getRootStructure(r, state)
                if (!structure || !SO.Molecule.Structure.is(structure.obj)) continue

                const unitTransforms = this.getUnitTransforms(structure.obj)
                const d = structure.obj.data.boundary.sphere.radius * (p / 100)

                if (!this.updatedUnitTransforms.has(structure.obj)) {
                    for (let i = 0, il = structure.obj.data.units.length; i < il; ++i) {
                        const u = structure.obj.data.units[i]
                        Vec3.transformMat4(this.centerVec, u.lookup3d.boundary.sphere.center, u.conformation.operator.matrix)

                        Vec3.sub(this.transVec, this.centerVec, structure.obj.data.boundary.sphere.center)
                        Vec3.setMagnitude(this.transVec, this.transVec, d)
                        Mat4.fromTranslation(this.animMat, this.transVec)

                        unitTransforms.setTransform(this.animMat, u)
                    }
                    this.updatedUnitTransforms.add(structure.obj)
                }

                r.obj.data.setState({ unitTransforms })
                this.ctx.canvas3d.add(r.obj.data)
            }
            this.ctx.canvas3d.requestDraw(true)
        }

        animateExplode(play: boolean) {
            if (play) {
                const animateExplode = (t: number) => {
                    this.explode((Math.sin(t * 0.001) + 1) * 50)
                    this.explodeAnimHandle = requestAnimationFrame(animateExplode)
                }
                this.explodeAnimHandle = requestAnimationFrame(animateExplode)
            } else {
                cancelAnimationFrame(this.explodeAnimHandle)
            }
        }

        register(): void { }

        update(p: StructureAnimationProps) {
            let updated = !PD.areEqual(StructureAnimationParams, this.params, p)
            if (this.params.rotate !== p.rotate) {
                this.params.rotate = p.rotate
                this.animateRotate(this.params.rotate)
            }
            if (this.params.explode !== p.explode) {
                this.params.explode = p.explode
                this.animateExplode(this.params.explode)
            }
            if (this.params.rotateValue !== p.rotateValue) {
                this.params.rotateValue = p.rotateValue
                this.rotate(degToRad(this.params.rotateValue))
            }
            if (this.params.explodeValue !== p.explodeValue) {
                this.params.explodeValue = p.explodeValue
                this.explode(this.params.explodeValue)
            }
            return updated;
        }

        unregister() {
            cancelAnimationFrame(this.rotateAnimHandle)
            cancelAnimationFrame(this.explodeAnimHandle)
        }
    },
    params: () => StructureAnimationParams
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