#version 330 compatibility

uniform float	uTime;		// "Time", from Animate( )
uniform float uA;

out vec2  	vST;		// texture coords

const float PI = 	3.14159265;
const float AMP = 	0.2;		// amplitude
const float W = 	2.;		// frequency


out  vec3  vN;		// normal vector
out  vec3  vL;		// vector from point to light
out  vec3  vE;		// vector from point to eye

vec3 LightPosition = vec3(  0., 5., 5. );
float a = .5;
float b = .4;
float c = .6;
void
main( )
{ 

	vec3 vert = gl_Vertex.xyz;
	vST = gl_MultiTexCoord0.st;

	//<< change vert to perform vertex distortion >>  ??
	if( vert.x < a)
		vert.x = vert.x + (vert.x-5.)*uA;
	else
		vert.x = vert.x + (vert.x+5.)*uA;

	if( vert.y > b)
		vert.y = vert.y + (vert.y-5.)*uA;
	else
		vert.y = vert.y + (vert.y+5.)*uA;

	vert.z = vert.z + (vert.z-5.)*uA;




	vec4 ECposition = gl_ModelViewMatrix * vec4( vert, 1. );
	vN = normalize( gl_NormalMatrix * gl_Normal );	// normal vector
	vL = LightPosition - ECposition.xyz;		// vector from the point
							// to the light position
	vE = vec3( 0., 0., 0. ) - ECposition.xyz;	// vector from the point
							// to the eye position 
	//gl_Position = gl_ModelViewProjectionMatrix * vec4( vert, 1. );

	gl_Position = gl_ModelViewProjectionMatrix * vec4( vert, 1. );
}
