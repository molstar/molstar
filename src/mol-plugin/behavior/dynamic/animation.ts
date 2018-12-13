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
    explode: PD.Boolean(false)
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

// TODO this is just for testing purposes
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

        rotate(play: boolean) {
            if (play) {
                const state = this.ctx.state.dataState
                const reprs = state.select(q => q.rootsOfType(PluginStateObject.Molecule.Representation3D));
                const rotate = (t: number) => {
                    const rad = degToRad((t / 10) % 360)
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
                    this.rotateAnimHandle = requestAnimationFrame(rotate)
                }
                this.rotateAnimHandle = requestAnimationFrame(rotate)
            } else {
                cancelAnimationFrame(this.rotateAnimHandle)
            }
        }

        explode(play: boolean) {
            if (play) {
                const state = this.ctx.state.dataState
                const reprs = state.select(q => q.rootsOfType(PluginStateObject.Molecule.Representation3D));
                const explode = (t: number) => {
                    const d = (Math.sin(t * 0.001) + 1) * 5
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
                    this.explodeAnimHandle = requestAnimationFrame(explode)
                }
                this.explodeAnimHandle = requestAnimationFrame(explode)
            } else {
                cancelAnimationFrame(this.explodeAnimHandle)
            }
        }

        register(): void { }

        update(p: StructureAnimationProps) {
            let updated = PD.areEqual(StructureAnimationParams, this.params, p)
            if (this.params.rotate !== p.rotate) {
                this.params.rotate = p.rotate
                this.rotate(this.params.rotate)
            }
            if (this.params.explode !== p.explode) {
                this.params.explode = p.explode
                this.explode(this.params.explode)
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