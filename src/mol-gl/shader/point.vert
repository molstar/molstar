precision mediump float;
uniform mat4 projection, view;
attribute vec3 position;

varying vec3 vPosition;

void main(){
    gl_PointSize = 20.0;
    gl_Position = projection * view * vec4(position, 1.0);
}