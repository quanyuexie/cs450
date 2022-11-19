#version 330 compatibility

uniform float	uTime;		// "Time", from Animate( )


uniform float uC;
uniform float   uKa =.2;
uniform float   uKd=.6;
uniform float   uKs = .2;		// coefficients of each type of lighting -- make sum to 1.0
uniform vec3    uColor = {1,1,0};			// object color
uniform vec3    uSpecularColor={1,1,1};		// light color
uniform float   uShininess=40;		// specular exponent

in  vec2  vST;			// texture coords
in  vec3  vN;			// normal vector
in  vec3  vL;			// vector from point to light
in  vec3  vE;			// vector from point to eye


void
main( )
{
	vec3 Normal = normalize(vN);
	vec3 Light     = normalize(vL);
	vec3 Eye        = normalize(vE);

	vec3 myColor = uColor;

	if((.4 -.1*uC)>vST.t && (.2+.1*uC)< vST.t ||
	(.8 -.1*uC)> vST.t && (.6+.1*uC) <vST.t )
	{
		myColor = vec3( 1., 0., 1. );
	}

	vec3 ambient = uKa * myColor;

	float d = max( dot(Normal,Light), 0. );       // only do diffuse if the light can see the point
	vec3 diffuse = uKd * d * myColor;

	float s = 0.;
	if( dot(Normal,Light) > 0. )	          // only do specular if the light can see the point
	{
		vec3 ref = normalize(  reflect( -Light, Normal )  );
		s = pow( max( dot(Eye,ref),0. ), uShininess );
	}
	vec3 specular = uKs * s * uSpecularColor;
	gl_FragColor = vec4( ambient + diffuse + specular,  1. );



}