import { Camera } from '../../../mol-canvas3d/camera';
import { Mat3, Vec3 } from '../../../mol-math/linear-algebra';
import { PluginCommands } from '../../../mol-plugin/commands';
import { PluginContext } from '../../../mol-plugin/context';
import { PluginStateObject } from '../../objects';


export function resetAxes(plugin: PluginContext) {
    if (!plugin.canvas3d) return;

    PluginCommands.Camera.SetSnapshot(plugin, { snapshot: cameraSetRotation(plugin.canvas3d.camera.getSnapshot(), Mat3.identity()) }); 
    // will probably not work for Headless plugin because of asynchronicity!

    // adjustCamera(plugin, old => cameraSetRotation(old, Mat3.identity()), 0);
}

export function orientAxes(plugin: PluginContext) {
    const structs = plugin.state.data.selectQ(q => q.ofType(PluginStateObject.Molecule.Structure));
    const rootStructs = structs.filter(s => s.obj && !s.obj?.data.parent);
    console.log('structs:', structs);
    console.log('rootStructs:', rootStructs);
    adjustCamera;
}



/** Return a new camera snapshot with the same target and camera distance from the target as `old`
 * but with diferent orientation. 
 * The actual rotation applied to the camera is the invert of `rotation`, 
 * which creates the same effect as if `rotation` were applied to the whole scene without moving the camera.
 * The rotation is relative to the default camera orientation (not to the current orientation). */
function cameraSetRotation(old: Camera.Snapshot, rotation: Mat3): Camera.Snapshot {
    const cameraRotation = Mat3.invert(Mat3(), rotation);
    const dist = Vec3.distance(old.position, old.target);
    const relPosition = Vec3.transformMat3(Vec3(), Vec3.create(0, 0, dist), cameraRotation);
    const newUp = Vec3.transformMat3(Vec3(), Vec3.create(0, 1, 0), cameraRotation);
    const newPosition = Vec3.add(Vec3(), old.target, relPosition);
    return { ...old, position: newPosition, up: newUp };
}

/** Apply `change` to the camera snapshot (i.e. target, position, orientation) in a plugin.
 * The `change` function will get the current camera snapshot and the result of the function will be used as the new snapshot. */
function adjustCamera(plugin: PluginContext, change: (old: Camera.Snapshot) => Camera.Snapshot, durationMs?: number) {
    if (!plugin.canvas3d) throw new Error('plugin.canvas3d is undefined');
    plugin.canvas3d.commit(true);
    const oldSnapshot = plugin.canvas3d.camera.getSnapshot();
    const newSnapshot = change(oldSnapshot);
    // plugin.canvas3d.camera.setState(newSnapshot, durationMs);
    plugin.managers.camera.setSnapshot(newSnapshot, durationMs);
    const checkSnapshot = plugin.canvas3d.camera.getSnapshot();
    if (durationMs && oldSnapshot.radius > 0 && !Camera.areSnapshotsEqual(newSnapshot, checkSnapshot)) {
        console.error('Error: The camera has not been adjusted correctly.');
        console.error('Required:');
        console.error(newSnapshot);
        console.error('Real:');
        console.error(checkSnapshot);
        throw new Error(`AssertionError: The camera has not been adjusted correctly.`);
    }
}
