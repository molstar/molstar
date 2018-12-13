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

const StructureAnimationParams = {
    rotate: PD.Boolean(false),
    rotateValue: PD.Numeric(0, { min: 0, max: 360, step: 0.1 }),
    explode: PD.Boolean(false),
    explodeValue: PD.Numeric(0, { min: 0, max: 100, step: 0.1 }),
}
type StructureAnimationProps = PD.Values<typeof StructureAnimationParams>

function getStructure(root: StateObjectCell, state: State) {
    const parent = StateSelection.findAncestorOfType(state.tree, state.cells, root.transform.ref, [PluginStateObject.Molecule.Structure])
    return parent && parent.obj ? parent.obj as PluginStateObject.Molecule.Structure : undefined
}

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
    return parent && parent.obj ? parent.obj as PluginStateObject.Molecule.Structure : undefined
}

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

        private rotateAnimHandle = -1
        private explodeAnimHandle = -1

        constructor(protected ctx: PluginContext, protected params: StructureAnimationProps) {
            super(ctx, params)
            this.update(params)
        }

        rotate(rad: number) {
            const state = this.ctx.state.dataState
            const reprs = state.select(q => q.rootsOfType(PluginStateObject.Molecule.Representation3D));
            Mat4.rotate(this.rotMat, this.tmpMat, rad, this.rotVec)
            for (const r of reprs) {
                if (!SO.isRepresentation3D(r.obj)) return
                const structure = getRootStructure(r, state)
                if (!structure) continue

                Vec3.negate(this.transVec, Vec3.copy(this.transVec, structure.data.boundary.sphere.center))
                Mat4.fromTranslation(this.transMat, this.transVec)
                Mat4.mul(this.animMat, this.rotMat, this.transMat)

                Vec3.copy(this.transVec, structure.data.boundary.sphere.center)
                Mat4.fromTranslation(this.transMat, this.transVec)
                Mat4.mul(this.animMat, this.transMat, this.animMat)

                r.obj.data.setState({ transform: this.animMat })
                this.ctx.canvas3d.add(r.obj.data)
                this.ctx.canvas3d.requestDraw(true)
            }
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

        explode(d: number) {
            const state = this.ctx.state.dataState
            const reprs = state.select(q => q.rootsOfType(PluginStateObject.Molecule.Representation3D));
            for (const r of reprs) {
                if (!SO.isRepresentation3D(r.obj)) return
                const structure = getStructure(r, state)
                if (!structure) continue
                const rootStructure = getRootStructure(r, state)
                if (!rootStructure) continue

                Vec3.sub(this.transVec, structure.data.boundary.sphere.center, rootStructure.data.boundary.sphere.center)
                Vec3.setMagnitude(this.transVec, this.transVec, d)
                Mat4.fromTranslation(this.animMat, this.transVec)

                r.obj.data.setState({ transform: this.animMat })
                this.ctx.canvas3d.add(r.obj.data)
                this.ctx.canvas3d.requestDraw(true)
            }
        }

        animateExplode(play: boolean) {
            if (play) {
                const animateExplode = (t: number) => {
                    this.explode((Math.sin(t * 0.001) + 1) * 5)
                    this.explodeAnimHandle = requestAnimationFrame(animateExplode)
                }
                this.explodeAnimHandle = requestAnimationFrame(animateExplode)
            } else {
                cancelAnimationFrame(this.explodeAnimHandle)
            }
        }

        register(): void { }

        update(p: StructureAnimationProps) {
            let updated = PD.areEqual(StructureAnimationParams, this.params, p)
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