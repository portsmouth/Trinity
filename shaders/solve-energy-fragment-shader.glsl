precision highp float;

uniform sampler2D Q;

uniform int Nr;
uniform int Ny;
uniform float dt;
uniform float dr;
uniform float dy;
uniform float gammaMO; // (gamma-1)
uniform float g;        // gravitational acceleration

out vec4 Qnext;

#define DENOM_EPS 1.0e-7

// energy injection term
float G(in ivec2 X)
{
	// @todo: inject non-zero energy in the first timestep
	//        within a sphere (which encloses at least one voxel center)
	// @todo: for now, place this detonation point in the grid center voxel
	return 0.0;
}

// Equation of state
float pressure(in vec4 Q_)
{
	float rho = Q_.x; // rho
	float E   = Q_.y; // total energy density
	float vr  = Q_.z / max(DENOM_EPS, rho);
	float vz  = Q_.w / max(DENOM_EPS, rho);
	// p = (gamma-1) * rho * e, where e = E/rho - v^2/2 (specific internal energy)
	float p = gammaMO * max(0.0, E - 0.5*rho*(vr*vr + vz*vz));	
	return p;
}

// r-flux
vec4 R(in vec4 Q_)
{
	float q0 = Q_.x; // rho
	float q1 = Q_.y; // E
	float q2 = Q_.z; // rho * vr
	float q3 = Q_.w; // rho * vz
	float p = pressure(Q_);
	vec4 R_;
	R_.x = q2;
	R_.y = (q1 + p) * q2/q0;
	R_.z = p + q2*q2/q0;
	R_.w = q2*q3 / q0;
	return R_;
}

// y-flux
vec4 Y(in vec4 Q_)
{
	float q0 = Q_.x; // rho
	float q1 = Q_.y; // E
	float q2 = Q_.z; // rho * vr
	float q3 = Q_.w; // rho * vz
	float p = pressure(Q_);
	vec4 Y_;
	Y_.x = q3;
	Y_.y = (q1 + p) * q3/q0;
	Y_.z = q2*q3/q0;
	Y_.w = p + q3*q3 / q0;
	return Y_;
}


// source term
vec4 S(in vec4 Q_, in ivec2 X)
{
	float r = (0.5 + float(X.x))*dr;
	float q0 = Q_.x; // rho
	float q1 = Q_.y; // E
	float q2 = Q_.z; // rho * vr
	float q3 = Q_.w; // rho * vz
	float p = pressure(Q_);
	vec4 S_;
	S_.x = -q2/r;
	S_.y = -(q1 + p)*q2/(q0*r) - g*q3/q0 + G(X)/r;
	S_.z = -q2*q2 / (r*q0);
	S_.w = -q2*q3/(r*q0) - q0*g;
	return S_;
}	


void main()
{
	// Setup local stencil:
	ivec2 X = ivec2(gl_FragCoord.xy);
	int ir = X.x;
	int iy = X.y;
	ivec2 X_rp = ivec2(min(ir+1, Nr-1), iy);
	ivec2 X_rn = ivec2(max(ir-1, 0),    iy);
	ivec2 X_yp = ivec2(ir, min(iy+1, Ny-1));
	ivec2 X_yn = ivec2(ir, max(iy-1, 0));

	// Update Q via Lax-Friedrichs method:
	float lambda_r = dt/dr;
	float lambda_y = dt/dy;
	vec4 Q_   = texelFetch(Q, X,    0);
	vec4 Q_rp = texelFetch(Q, X_rp, 0);
	vec4 Q_rn = texelFetch(Q, X_rn, 0);
	vec4 Q_yp = texelFetch(Q, X_yp, 0);
	vec4 Q_yn = texelFetch(Q, X_yn, 0);
	vec4 Qavg = 0.25 * (Q_rp + Q_rn + Q_yp + Q_yn);

	Qnext = Qavg - 0.5 * lambda_r * ( R(Q_rp) - R(Q_rn) ) 
	             - 0.5 * lambda_y * ( Y(Q_yp) - Y(Q_yn) );
	             + dt * S(Q_, X);
}





