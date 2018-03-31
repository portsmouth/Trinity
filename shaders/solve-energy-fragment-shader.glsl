precision highp float;

uniform sampler2D Q;

uniform float dt;
uniform float dr;
uniform float dz;
uniform float gammaMO; // (gamma-1)
uniform float g;        // gravitational acceleration

out vec4 Qnext;


// energy injection term
float G(in vec2 X)
{
	// @todo: inject non-zero energy in the first timestep
	//        within a sphere (which encloses at least one voxel center)
	// @todo: for now, place this detonation point in the grid center voxel
	return 0.0;
}

// Equation of state
float pressure(float rho, float E)
{
	float e = max(0.0, E - 0.5*(vr*vr + vz*vz));	
	return gammaMO * rho * e;
}

// r-flux
vec4 R(in vec4 Q_)
{
	float q0 = Q_.x; // rho
	float q1 = Q_.y; // E
	float q2 = Q_.z; // rho * vr
	float q3 = Q_.w; // rho * vz
	float p = pressure(q0, q1);
	vec4 R_;
	R_.x = q2;
	R_.y = (q1 + p) * q2/q0;
	R_.z = p + q2*q2/q0;
	R_.w = q2*q3 / q0;
	return R;
}

// z-flux
vec4 Z(in vec4 Q_)
{
	float q0 = Q_.x; // rho
	float q1 = Q_.y; // E
	float q2 = Q_.z; // rho * vr
	float q3 = Q_.w; // rho * vz
	float p = pressure(q0, q1);
	vec4 Z_;
	Z_.x = q3;
	Z_.y = (q1 + p) * q3/q0;
	Z_.z = q2*q3/q0;
	Z_.w = p + q3*q3 / q0;
	return R;
}


// source term
vec4 S(in vec4 Q_, in vec2 X)
{
	float r = (0.5 + X.x)*dr;
	float q0 = Q_.x; // rho
	float q1 = Q_.y; // E
	float q2 = Q_.z; // rho * vr
	float q3 = Q_.w; // rho * vz
	float p = pressure(q0, q1);
	vec4 S_;
	S_.x = - q2/r;
	S_.y = -(q1 + p)*q2/(q0*r) - g*q3/q0 + G(X)/r;
	S_.z = -q2*q2 / (r*q0);
	S_.w = -q2*q3/(r*q0) - q0*g;
	return S_;
}	


void main()
{
	// Setup local stencil:
	vec2 size = textureSize(U_v);
	int Nr = size.x;
	int Nz = size.y;
	ivec2 X = gl_FragCoord.xy;
	int ir = X.x;
	int iz = X.y;
	ivec2 X_rp = ivec2(min(ir+1, Nr-1), iz);
	ivec2 X_rn = ivec2(max(ir-1, 0),    iz);
	ivec2 X_zp = ivec2(ir, min(iz+1, Nz-1));
	ivec2 X_zn = ivec2(ir, max(iz-1, 0));

	// Update Q via Lax-Friedrichs method:
	float lambda_r = dt/dr;
	float lambda_z = dt/dz;
	vec4 Q_ = texelFetch(Q, X).rgb;
	vec4 Q_rp = texelFetch(Q, X_rp).rgb;
	vec4 Q_rn = texelFetch(Q, X_rn).rgb;
	vec4 Q_zp = texelFetch(Q, X_zp).rgb;
	vec4 Q_zn = texelFetch(Q, X_zn).rgb;
	vec4 Qavg = 0.25 * (Q_rp + Q_rn + Q_zp + Q_zn);
	Qnext = Qavg - 0.5*lambda_r*(R(Q_rp)-R(Q_rn)) - 0.5*lambda_z*(R(Q_zp)-R(Q_zn)) + dt*S(Q_, X);
}





