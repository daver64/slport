/***
*
* License: BSD.
*
*/
#include "template_3dmath.h"
#include <random>


const double_t PI = 3.14159265358979323846264338327950288;
const double_t TAU = 2.0 * PI;
const double_t PIDIV2 = PI / 2.0;
double_t ROUNDING_ERROR = std::numeric_limits<double_t>::epsilon();

const float_t PIf = 3.14159265358979323846;
const float_t TAUf = 2.0f * PIf;
const float_t PIDIV2f = PIf / 2.0f;
float_t ROUNDING_ERRORf = std::numeric_limits<float_t>::epsilon();
float_t vec3_distance(Vec3<int32_t> v1, Vec3<int32_t> v2)
{
	const Vec3<float_t> fv1{ (float_t)v1.x,(float_t)v1.y,(float_t)v1.z };
	const Vec3<float_t> fv2{ (float_t)v2.x,(float_t)v2.y,(float_t)v2.z };
	return (float_t)vec3_distance(fv1, fv2);
}
Vec3d vec3i_to_vec3d(const Vec3i& v)
{
	Vec3d result{ (double_t)v.x,(double_t)v.y,(double_t)v.z };
	return result;
}
Vec3d vec3f_to_vec3d(const Vec3f& v)
{
	Vec3d result{ (double_t)v.x,(double_t)v.y,(double_t)v.z };
	return result;
}

Vec3i vec3d_to_vec3i(const Vec3d& v)
{
	Vec3i result{ (int32_t)v.x,(int32_t)v.y,(int32_t)v.z };
	return result;
}
Vec3i vec3f_to_vec3i(const Vec3f& v)
{
	Vec3i result{ (int32_t)v.x,(int32_t)v.y,(int32_t)v.z };
	return result;
}

Vec3f vec3i_to_vec3f(const Vec3i& v)
{
	Vec3f result{ (float_t)v.x,(float_t)v.y,(float_t)v.z };
	return result;
}
Vec3f vec3d_to_vec3f(const Vec3d& v)
{
	Vec3f result{ (float_t)v.x,(float_t)v.y,(float_t)v.z };
	return result;
}

void quadratic_plot(Vec2f startpos, Vec2f controlpos, Vec2f endpos, std::vector< Vec2f>& plotresult, int32_t numsegments)
{
	for (int32_t i = 0; i <= numsegments; ++i)
	{
		const double_t t = (double_t)i / (double_t)numsegments;
		const double_t a = std::pow((1.0 - t), 2.0);
		const double_t b = 2.0 * t * (1.0 - t);
		const double_t c = std::pow((float_t)t, 2.0);
		const double_t x = a * startpos.x + b * controlpos.x + c * endpos.x;
		const double_t y = a * startpos.y + b * controlpos.y + c * endpos.y;
		plotresult.emplace_back(Vec2f((float_t)x, (float_t)y));
	}
}

void cubic_plot(Vec2f startpos, Vec2f controlpos1, Vec2f endpos, Vec2f controlpos2, std::vector< Vec2f>& plotresult, int32_t numsegments)
{
	for (int32_t i = 0; i <= numsegments; ++i)
	{
		const float_t t = (float_t)i / (float_t)numsegments;
		const float_t a = powf((1.0f - t), 3.0f);
		const float_t b = 3.0f * t * powf((1.0f - t), 2.0f);
		const float_t c = 3.0f * powf(t, 2.0f) * (1.0f - t);
		const float_t d = powf(t, 3.0f);
		const float_t x = a * startpos.x + b * controlpos1.x + d * endpos.x + c * controlpos2.x;
		const float_t y = a * startpos.y + b * controlpos1.y + d * endpos.y + c * controlpos2.y;
		plotresult.emplace_back(Vec2f((float_t)x, (float_t)y));
	}
}

float_t rand_gauss(float_t mu, float_t sigma)
{
	static const float_t eps = 0.001f;
	static float_t z1;
	static int32_t generate = false;
	generate = !generate;
	if (!generate)
		return z1 * sigma * mu;
	float_t u1, u2;
	do
	{
		u1 = rand_range(0, 1);
		u2 = rand_range(0, 1);
	} while (u1 < eps);
	float_t z0;
	z0 = sqrtf(-2.0f * logf(u1)) * cosf((float_t)TAU * u2);
	z1 = sqrtf(-2.0f * logf(u1)) * sinf((float_t)TAU * u2);
	return z0 * sigma + mu;

}

struct xorshift64s_state
{
	uint64_t a;
};

uint64_t xorshift64s(xorshift64s_state& xorstate)
{
	uint64_t x = xorstate.a;	/* The state must be seeded with a nonzero value. */
	x ^= x >> 12; // a
	x ^= x << 25; // b
	x ^= x >> 27; // c
	xorstate.a = x;
	return x * 0x2545F4914F6CDD1DULL;
}
struct xorshift128p_state
{
	uint64_t a, b;
};
uint64_t xorshiftp(xorshift128p_state& xorstate)
{
	uint64_t t = xorstate.a;
	uint64_t const s = xorstate.b;
	xorstate.a = s;
	t ^= t << 23;		// a
	t ^= t >> 17;		// b
	t ^= s ^ (s >> 26);	// c
	xorstate.b = t;
	return t + s;
}
uint32_t xorshift(uint32_t& xorstate)
{
	uint32_t x = xorstate;
	x ^= x << 13;
	x ^= x >> 17;
	x ^= x << 5;
	return xorstate = x;
}


float_t rand_range(float_t a, float_t b)
{
	static std::random_device rd; 
	static std::mt19937 eng(rd()); 
	std::uniform_real_distribution<> distr(a, b);
	return (float_t)distr(eng);
}
int32_t int_rand_range(int32_t a, int32_t b)
{
	static std::random_device rd;
	static std::mt19937 eng(rd());
	std::uniform_int_distribution<> distr(a, b);
	return distr(eng);
}

float_t rand_range_gauss(float_t a, float_t stddev)
{
	static std::random_device rd{};
	static std::mt19937 eng{ rd() };
	std::normal_distribution<float_t> distr{ a, stddev };
	return distr(eng);
}
int32_t int_rand_range_gauss(int32_t a, float_t stddev)
{
	static std::random_device rd{};
	static std::mt19937 eng{ rd() };
	std::normal_distribution<float_t> distr{ (float_t)a, stddev };
	return (int32_t)distr(eng);
}
Vec2f get_line_intersection(Vec2f v0, Vec2f v1, Vec2f v2, Vec2f v3,bool& success)
{
	Vec2f vr;
	success = get_line_intersection(v0.x, v0.y, v1.x, v1.y, v2.x, v2.y, v3.x, v3.y, vr.x, vr.y);
	return vr;
}

int32_t get_line_intersection(float_t p0_x, float_t p0_y,
	float_t p1_x, float_t p1_y,
	float_t p2_x, float_t p2_y,
	float_t p3_x, float_t p3_y,
	float_t& i_x, float_t& i_y)
{
	float_t s1_x, s1_y, s2_x, s2_y;
	s1_x = p1_x - p0_x;     s1_y = p1_y - p0_y;
	s2_x = p3_x - p2_x;     s2_y = p3_y - p2_y;
	float_t s, t;
	s = (-s1_y * (p0_x - p2_x) + s1_x * (p0_y - p2_y)) / (-s2_x * s1_y + s1_x * s2_y);
	t = (s2_x * (p0_y - p2_y) - s2_y * (p0_x - p2_x)) / (-s2_x * s1_y + s1_x * s2_y);
	if (s >= 0 && s <= 1 && t >= 0 && t <= 1)
	{
		i_x = p0_x + (t * s1_x);
		i_y = p0_y + (t * s1_y);
		return true;
	}
	i_x = NAN;
	i_y = NAN;
	return false;
}

bool box_box_intersect(AABBoxd& abox, AABBoxd& bbox)
{
	return (
		abox.min_edge.x <= bbox.max_edge.x &&
		abox.max_edge.x >= bbox.min_edge.x &&
		abox.min_edge.y <= bbox.max_edge.y &&
		abox.max_edge.y >= bbox.min_edge.y &&
		abox.min_edge.z <= bbox.max_edge.z &&
		abox.max_edge.z >= bbox.min_edge.z
		);
}
bool box_box_intersect(AABBoxf& abox, AABBoxf& bbox)
{
	return (
		abox.min_edge.x <= bbox.max_edge.x &&
		abox.max_edge.x >= bbox.min_edge.x &&
		abox.min_edge.y <= bbox.max_edge.y &&
		abox.max_edge.y >= bbox.min_edge.y &&
		abox.min_edge.z <= bbox.max_edge.z &&
		abox.max_edge.z >= bbox.min_edge.z
		);
}
bool box_sphere_intersect(AABBoxd& abox, Vec3d position, double_t radius)
{
	double_t x = std::max(abox.min_edge.x, std::min(position.x, abox.max_edge.x));
	double_t y = std::max(abox.min_edge.y, std::min(position.y, abox.max_edge.y));
	double_t z = std::max(abox.min_edge.z, std::min(position.z, abox.max_edge.z));
	double_t distance = sqrt(
		((x - position.x) * (x - position.x)) +
		((y - position.y) * (y - position.y)) +
		((z - position.z) * (z - position.z))
	);
	return distance < radius;
}
bool box_sphere_intersect(AABBoxf& abox, Vec3f position, float_t radius)
{
	double_t x = std::max(abox.min_edge.x, std::min(position.x, abox.max_edge.x));
	double_t y = std::max(abox.min_edge.y, std::min(position.y, abox.max_edge.y));
	double_t z = std::max(abox.min_edge.z, std::min(position.z, abox.max_edge.z));
	double_t distance = sqrt(
		((x - position.x) * (x - position.x)) +
		((y - position.y) * (y - position.y)) +
		((z - position.z) * (z - position.z))
	);
	return distance < radius;
}


double_t gcircle_distance(double_t lon1, double_t lat1, double_t lon2, double_t lat2, double_t R)
{
	const double_t theta1 = degtorad(lat1);
	const double_t theta2 = degtorad(lat2);
	const double_t delta_theta = degtorad(lat2 - lat1);
	const double_t delta_tau = degtorad(lon2 - lon1);
	const double_t a = sin(delta_theta / 2.0) * sin(delta_theta / 2.0) +
		cos(theta1) * cos(theta2) * sin(delta_tau / 2) * sin(delta_tau / 2.0);
	const double_t c = 2 * atan2(sqrt(a), sqrt(1 - a));
	double_t d = R * c;
	return d;
}

double_t gcircle_bearing(double_t from_lon, double_t from_lat, double_t to_lon, double_t to_lat)
{
	const double_t y = sin(degtorad(to_lon - from_lon)) * cos(degtorad(to_lat));
	const double_t x = cos(degtorad(from_lat)) * sin(degtorad(to_lat)) - sin(degtorad(from_lat)) * cos(degtorad(to_lat)) * cos(degtorad(to_lon - from_lon));
	double_t b = radtodeg(atan2(y, x));
	if (b < 0)
		b = 360 + b;
	return b;
}

void gcircle_destination(double_t from_lon, double_t from_lat, double_t bearing, double_t distance, double_t& to_lon, double_t& to_lat, double_t R)
{
	const double_t to_lat1 = asin(sin(degtorad(from_lat)) * cos(distance / R) +
		cos(degtorad(from_lat)) * sin(distance / R) * cos(degtorad(bearing)));
	const double_t lon = degtorad(from_lon) + atan2(sin(degtorad(bearing)) *
		sin(distance / R) * cos(degtorad(from_lat)),
		cos(distance / R) - sin(degtorad(from_lat)) * sin(degtorad(to_lat1)));
	to_lat = radtodeg(to_lat1);
	to_lon = radtodeg(lon);
}

