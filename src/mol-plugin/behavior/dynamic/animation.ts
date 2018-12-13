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
    rotate: PD.Boolean(false)
}
type StructureAnimationProps = PD.Values<typeof StructureAnimationParams>

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
    if (!parent || !parent.obj) return
    return parent.obj as PluginStateObject.Molecule.Structure
}

// TODO this is just for testing purposes
export const StructureAnimation = PluginBehavior.create<StructureAnimationProps>({
    name: 'structure-animation',
    display: { name: 'Structure Animation', group: 'Animation' },
    ctor: class extends PluginBehavior.Handler<StructureAnimationProps> {
        private tmpMat = Mat4.identity()
        private rotMat = Mat4.identity()
        private transMat = Mat4.identity()
        private rotAnimMat = Mat4.identity()
        private transVec = Vec3.zero()
        private rotVec = Vec3.create(0, 1, 0)

        private rotateAnimHandle = -1

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
                        Mat4.mul(this.rotAnimMat, this.rotMat, this.transMat)

                        Vec3.copy(this.transVec, structure.data.boundary.sphere.center)
                        Mat4.fromTranslation(this.transMat, this.transVec)
                        Mat4.mul(this.rotAnimMat, this.transMat, this.rotAnimMat)

                        r.obj.data.setState({ transform: this.rotAnimMat })
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

        register(): void { }

        update(p: StructureAnimationProps) {
            let updated = this.params.rotate !== p.rotate
            this.params.rotate = p.rotate
            if (updated) {
                this.rotate(this.params.rotate)
            }
            return updated;
        }

        unregister() {
            cancelAnimationFrame(this.rotateAnimHandle)
        }
    },
    params: () => StructureAnimationParams
});