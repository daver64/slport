/***
 * Copyright 2004-2021, Dave Rowbotham and Toni Ylisirnio
 * All rights reserved.
 *
 * License: BSD.
 *
 */

#pragma once
#include <limits>
#include <cmath>
#include <cstdint>
#include <algorithm>
#include <vector>
#include <map>
#include <xmmintrin.h>
extern const double_t PI;
extern const double_t TAU;
extern const double_t PIDIV2;
extern double ROUNDING_ERROR;

extern const float_t PIf;
extern const float_t TAUf;
extern const float_t PIDIV2f;
extern float_t ROUNDING_ERRORf;

/*
util.
*/
inline float_t celsius_to_kelvin(float_t val) { return val + 273.15f; }
inline double_t celsius_to_kelvin(double_t val) { return val + 273.15; }
inline float_t lerp(float_t v1, float_t v2, float_t a) { return (v1 * (1.0f - a) + v2 * a); }
inline double_t lerp(double_t v1, double_t v2, double_t a) { return (v1 * (1.0 - a) + v2 * a); }

inline double_t ms_to_knots(double_t val) { return val * 1.9438444924406; }
inline float_t ms_to_knots(float_t val) { return val * 1.9438444924406f; }
inline double_t knots_to_ms(double_t val) { return val * 0.51444; }
inline float_t knots_to_ms(float_t val) { return val * 0.51444f; }
inline double_t metres_to_feet(double_t val) { return val * 3.28084; }
inline float_t metres_to_feet(float_t val) { return val * 3.2808f; }
inline double_t feet_to_metres(double_t val) { return val * 0.3048; }
inline float_t feet_to_metres(float_t val) { return val * 0.3048f; }
inline double_t meters_to_nm(double_t val) { return val * 0.000539957; }
inline float_t meters_to_nm(float_t val) { return val * 0.000539957f; }
inline float_t nm_to_meters(float_t val) { return val * 1852.0f; }
inline double_t nm_to_meters(double_t val) { return val * 1852.0; }
inline double_t degtorad(double_t val) { return val * 0.0174532925; }
inline float_t degtorad(float_t val) { return val * 0.0174532925f; }
inline double_t radtodeg(double_t val) { return val * 57.2957795; }
inline float_t radtodeg(float_t val) { return val * 57.2957795f; }

template <typename T, typename U, typename V>
inline T clamp(T a, U low, V high)
{
	T val = a < low ? low : a;
	return val > high ? high : val;
}

template <typename T = double_t>
inline T normalise_deg(T angle)
{
	angle = std::fmod(angle, 360.0);
	if (angle < 0)
		angle += 360.0;
	return angle;
}
template <typename T = int32_t>
inline double normalise_deg(T angle)
{
	return normalise_deg((double)angle);
}

template <typename T = double_t>
inline T normalize_rad(T angle)
{
	angle = std::fmod(angle, TAU);
	if (angle < 0.0)
		angle += TAU;
	return angle;
}
template <typename T = double_t>
inline bool dbl_equals(const T &a,
					   const T &b, const T &tolerance = ROUNDING_ERROR)
{
	return (a + tolerance >= b) && (a - tolerance <= b);
}
/*
 * Round down if direction (sign) is negative, round up if positive.
 * Purpose is to find the next integer transition based on
 * the delta direction.
 */
template <typename T = double_t>
inline T sign_round(const T &v, const T &sign)
{
	T rv{0.f}; // Note corner cases:
	// at eg. v==2.0,
	// next + transition is at floor(3.0) => 3.0 (correct)
	if (sign > 0)
		rv = std::floor(v + (T)1.0);
	// at eg. v==2.0, next - transition is
	// at ceil(1.0) => 1.0 (correct)
	else
		rv = std::ceil(v - (T)1.0);
	return rv;
}
template <typename T = double_t>
inline bool is_zero(const T &a, const T &tolerance = ROUNDING_ERROR)
{
	return std::fabs(a) <= tolerance;
}
template <typename T = double_t>
inline T alpha_function(const T &a, const T &x)
{
	return x * std::exp(-a * x);
}
template <typename T = double_t>
inline T gompertz_function(const T &k,
						   const T &a, const T &x, const T &y)
{
	T kf = (k - k * std::exp(-a * x)) / a;
	return (y * std::exp(kf));
}

template <typename T>
class Vec3;

template <typename T = double_t>
class Vec2
{
public:
	union
	{
		T d[2] = {0, 0};
		struct
		{
			T x, y;
		};
		struct
		{
			T u, v;
		};
		struct
		{
			T w, h;
		};
	};
	Vec2(const Vec2 &f) : x(f.x), y(f.y) {}
	Vec2(const T dat[2]) : x(dat[0]), y(dat[1]) {}
	Vec2(const T x = 0, const T y = 0) : x(x), y(y) {}
	Vec2(const Vec3<T> v) : x(v.x), y(v.y) {}
	T &operator[](size_t index)
	{
		return d[index];
	}
	const T length() const
	{
		return std::sqrt(x * x + y * y);
	}
	Vec2 &normalise()
	{
		const T len = std::sqrt(x * x + y * y);
		if (len != 0)
		{
			x /= len;
			y /= len;
		}
		return *this;
	}
	T dot(const Vec2 &v)
	{
		return x * v.x + y * v.y;
	}
	Vec2 operator+(const Vec2 &v) const
	{
		Vec2 res(*this);
		res.x += v.x;
		res.y += v.y;
		return res;
	}
	Vec2 operator-(const Vec2 &v) const
	{
		Vec2 res(*this);
		res.x -= v.x;
		res.y -= v.y;
		return res;
	}
	Vec2 operator-() const
	{
		Vec2 res(*this);
		res.x = -res.x;
		res.y = -res.y;
		return res;
	}
	Vec2 operator*(T v) const
	{
		Vec2 res(*this);
		res.x *= v;
		res.y *= v;
		return res;
	}
	Vec2 operator/(const Vec2 &v)
	{
		Vec2 res(*this);
		if (!is_zero(v.x))
			res.x /= v.x;
		if (!is_zero(v.y))
			res.y /= v.y;
		return res;
	}
	Vec2 operator=(const Vec2 &v)
	{
		x = v.x;
		y = v.y;
		return *this;
	}
	Vec2 operator+=(const Vec2 &v)
	{
		x += v.x;
		y += v.y;
		return *this;
	}
	Vec2 operator-=(const Vec2 &v)
	{
		x -= v.x;
		y -= v.y;
		return *this;
	}
	Vec2 operator*=(const T s)
	{
		x *= s;
		y *= s;
		return *this;
	}
	Vec2 operator/=(const T s)
	{
		if (!is_zero(s))
		{
			x /= s;
			y /= s;
		}
		return *this;
	}
	const bool operator<=(const Vec2 &v) const
	{
		return x <= v.x && y <= v.y;
	}
	const bool operator>=(const Vec2 &v) const
	{
		return x >= v.x && y >= v.y;
	}
	const bool operator<(const Vec2 &v) const
	{
		return x < v.x && y < v.y;
	}
	const bool operator>(const Vec2 &v) const
	{
		return x > v.x && y > v.y;
	}
	const bool operator==(const Vec2 &v) const
	{
		return this->equals(v);
	}
	const bool operator!=(const Vec2 &v) const
	{
		return !this->equals(v);
	}
	const bool equals(const Vec2 &v, const T tolerance = 0.0001) const
	{
		return dbl_equals(x, v.x, tolerance) &&
			   dbl_equals(y, v.y, tolerance);
	}
	void set(T a, T b)
	{
		x = a;
		y = b;
	}
	const T length_squared() const
	{
		return x * x + y * y;
	}
	const T manhattan(const Vec2 &v)
	{
		return x - v.x + y - v.y;
	}
	const bool is_between_points(const Vec2 &begin, const Vec2 &end) const
	{
		T f = (end - begin).length_squared();
		return get_distance_from_squared(begin) < f &&
			   get_distance_from_squared(end) < f;
	}
	const T get_distance_from_squared(const Vec2 &v) const
	{
		return Vec2(x - v.x, y - v.y).length_squared();
	}
	Vec2 &invert()
	{
		x *= 1.;
		y *= -1.;
		return *this;
	}
};
template <typename T = double_t>
class Vec3
{
public:
	// data is sized specifically to play nice with cache lines. this library itself
	// doesn't use the fourth float, it's for padding only.
	union
	{
		T d[4] = {0, 0, 0, 0};
		struct
		{
			T x, y, z, w;
		};
		struct
		{
			T h, p, b, w;
		};
		struct
		{
			T width, height, depth, w;
		};
	};
	Vec3() : x(0), y(0), z(0), w(0) {}
	Vec3(const Vec2<T> v) : x(v.x), y(v.y), z(0), w(0) {}
	Vec3(const Vec3<T> &v = Vec3<T>(0, 0, 0)) : x(v.x), y(v.y), z(v.z), w(0) {}
	Vec3(const T dat[3]) : x(dat[0]), y(dat[1]), z(dat[2]), w(0) {}
	Vec3(const T x, const T y, const T z = 0) : x(x), y(y), z(z), w(0) {}

	// Vec3(int32_t x, int32_t y,int32_t z ) :x(x), y(y), z(z), w(0) {}
	T &operator[](size_t index)
	{
		return d[index];
	}
	const double_t length() const
	{
		return (double_t)std::sqrt((double_t)(x * x + y * y + z * z));
	}
	const float_t lengthf() const
	{
		return (float_t)std::sqrt((float_t)(x * x + y * y + z * z));
	}
	Vec3 &normalise()
	{
		const T len = std::sqrt(x * x + y * y + z * z);
		if (len != 0)
		{
			x /= len;
			y /= len;
			z /= len;
		}
		return *this;
	}
	T dot(const Vec3 &v) const
	{
		return x * v.x + y * v.y + z * v.z;
	}
	Vec2<> to_Vec2()
	{
		Vec2<> res;
		res.x = x;
		res.y = y;
		return res;
	}
	Vec3 cross(const Vec3 &v) const
	{
		Vec3 res;
		res.x = y * v.z - z * v.y;
		res.y = z * v.x - x * v.z;
		res.z = x * v.y - y * v.x;
		return res;
	}
	Vec3 operator+(const Vec3 &v) const
	{
		Vec3 res(*this);
		res.x += v.x;
		res.y += v.y;
		res.z += v.z;
		return res;
	}
	Vec3 operator-(const Vec3 &v) const
	{
		Vec3 res(*this);
		res.x -= v.x;
		res.y -= v.y;
		res.z -= v.z;
		return res;
	}
	Vec3 operator-() const
	{
		Vec3 res(*this);
		res.x = -res.x;
		res.y = -res.y;
		res.z = -res.z;
		return res;
	}
	Vec3 operator*(T v) const
	{
		Vec3 res(*this);
		res.x *= v;
		res.y *= v;
		res.z *= v;
		return res;
	}
	Vec3 operator/(const Vec3 &v)
	{
		Vec3 res(*this);
		if (!is_zero(v.x))
			res.x /= v.x;
		if (!is_zero(v.y))
			res.y /= v.y;
		if (!is_zero(v.z))
			res.z /= v.z;
		return res;
	}
	Vec3 operator/(T v) const
	{
		Vec3 res(*this);
		if (!is_zero(v))
		{
			res.x /= v;
			res.y /= v;
			res.z /= v;
		}
		return res;
	}
	Vec3 operator+=(const Vec3 &v)
	{
		x += v.x;
		y += v.y;
		z += v.z;
		return *this;
	}
	Vec3 operator-=(const Vec3 &v)
	{
		x -= v.x;
		y -= v.y;
		z -= v.z;
		return *this;
	}
	Vec3 operator*=(const T s)
	{
		x *= s;
		y *= s;
		z *= s;
		return *this;
	}
	Vec3 operator/=(const T s)
	{
		if (!is_zero(s))
		{
			x /= s;
			y /= s;
			z /= s;
		}
		return *this;
	}
	const bool operator<=(const Vec3 &v) const
	{
		return x <= v.x && y <= v.y && z <= v.z;
	}
	const bool operator>=(const Vec3 &v) const
	{
		return x >= v.x && y >= v.y && z >= v.z;
	}
	const bool operator<(const Vec3 &v) const
	{
		return x < v.x && y < v.y && z < v.z;
	}
	const bool operator>(const Vec3 &v) const
	{
		return x > v.x && y > v.y && z > v.z;
	}
	const bool operator==(const Vec3 &v) const
	{
		return this->equals(v);
	}
	const bool operator!=(const Vec3 &v) const
	{
		return !this->equals(v);
	}
	const bool equals(const Vec3 &v, const T tolerance = 0.0001) const
	{
		return dbl_equals(x, v.x, tolerance) &&
			   dbl_equals(y, v.y, tolerance) &&
			   dbl_equals(z, v.z, tolerance);
	}
	void set(T a, T b, T c)
	{
		x = a;
		y = b;
		z = c;
	}
	const T length_squared() const
	{
		return x * x + y * y + z * z;
	}
	const T manhattan(const Vec3 &v)
	{
		return x - v.x + y - v.y + z - v.z;
	}
	const bool is_between_points(const Vec3 &begin, const Vec3 &end) const
	{
		T f = (end - begin).length_squared();
		return get_distance_from_squared(begin) < f &&
			   get_distance_from_squared(end) < f;
	}
	const T get_distance_from_squared(const Vec3 &v) const
	{
		return Vec3(x - v.x, y - v.y, z - v.z).length_squared();
	}
	Vec3 &invert()
	{
		x *= 1.;
		y *= -1.;
		z *= -1.;
		return *this;
	}
};

template <typename T = double_t>
inline Vec3<T> make_vector(const Vec3<T> &v1, const Vec3<T> &v2)
{
	return v2 - v1;
}
template <typename T = double_t>
inline Vec2<T> make_vector2(const Vec2<T> &v1, const Vec2<T> &v2)
{
	return v2 - v1;
}
template <typename T = double_t>
inline Vec3<T> mid_point(const Vec3<T> &v1, const Vec3<T> &v2)
{
	return (v1 + v2) * 0.5;
}
template <typename T = double_t>
inline Vec2<T> mid_point2(const Vec2<T> &v1, const Vec2<T> &v2)
{
	return (v1 + v2) * 0.5;
}
template <typename T = double_t>
inline Vec3<T> face_normal(const Vec3<T> &v1,
						   const Vec3<T> &v2, const Vec3<T> &v3)
{
	return (v1 - v2).cross((v1 - v3)).normalise();
}
template <typename T = double_t>
inline T vec2_distance(const Vec2<T> &v1, const Vec2<T> &v2)
{
	Vec2<T> r = v1 - v2;
	return r.length();
}
template <typename T = double_t, typename S = double_t>
inline T vec3_distance(const Vec3<T> &v1, const Vec3<S> &v2)
{
	Vec3<T> r = v1 - v2;
	return (T)r.length();
}
float_t vec3_distance(Vec3<int32_t> v1, Vec3<int32_t> v2);
extern float_t rand_range(float_t, float_t);
template <typename T = double_t>
inline Vec3<T> rand_vec3()
{
	const T z = rand_range(-1.0f, 1.0f);
	const T a = rand_range(0.0f, (float)TAU);
	const T r = (T)std::sqrt(1.0 - z * z);
	const T x = r * std::cos(a);
	const T y = r * std::sin(a);
	return Vec3<T>(x, y, z);
}

template <typename T = double_t>
class Polar
{
public:
	union
	{
		T d[4] = {0, 0, 0, 0};
		struct
		{
			T ra, dec, rad, w;
		};
		struct
		{
			T lon, lat, alt, p;
		};
	};
	Polar(const T &lon = 0,
		  const T &lat = 0,
		  const T &alt = 0) : lon(lon), lat(lat), alt(alt) {}
	Polar(const T p[3]) : lon(p[0]), lat(p[1]), alt(p[0]) {}
	Polar(const Vec3<T> &pos)
	{
		rad = std::sqrt(pos.x * pos.x + pos.y * pos.y + pos.z * pos.z);
		lon = std::atan(pos.y / pos.x);
		lat = std::atan(std::sqrt(pos.x * pos.x + pos.y * pos.y) / pos.z);
	}
	T &operator[](size_t index)
	{
		return d[index];
	}
	Vec3<T> polar_to_cartesian(const Polar &pos)
	{
		return Vec3<T>(std::cos(pos.lon) * std::cos(pos.lat) * pos.rad,
					   std::sin(pos.lon) * std::cos(pos.lat) * pos.rad,
					   std::sin(pos.lat) * pos.rad);
	}
	Polar<T> cartesian_to_polar(const Vec3<T> &pos)
	{
		return Polar(pos);
	}
	Vec3<T> to_cartesian()
	{
		return polar_to_cartesian(*this);
	}

	void from_cartesian(const Vec3<> &pos)
	{
		*this = cartesian_to_polar(pos);
	}
};

template <typename T = double_t>
class Line
{
public:
	Vec2<T> start;
	Vec2<T> end;
	Line() : start(0, 0), end(0, 0) {}
	Line operator+(const Vec2<T> &point) const
	{
		return Line(start + point, end + point);
	}
};

template <typename T = double_t>
class Line3
{
public:
	Vec3<T> start;
	Vec3<T> end;
	Line3() : start(0, 0, 0), end(1, 1, 1) {}
	Line3(double_t xa, double_t ya, double_t za, double_t xb, double_t yb, double_t zb)
		: start(xa, ya, za), end(xb, yb, zb) {}
	Line3(const Vec3<T> &start, const Vec3<T> &end) : start(start), end(end) {}
	Line3 operator+(const Vec3<T> &point) const
	{
		return Line3(start + point, end + point);
	}
	Line3 &operator+=(const Vec3<T> &point)
	{
		start += point;
		end += point;
		return *this;
	}
	Line3 operator-(const Vec3<T> &point) const
	{
		return Line3(start - point, end - point);
	}
	Line3 &operator-=(const Vec3<T> &point)
	{
		start -= point;
		end -= point;
		return *this;
	}
	bool operator==(const Line3 &other) const
	{
		return !(start.equals(other.start) && end.equals(other.end)) || (end.equals(other.start) && start.equals(other.end));
	}
	bool operator!=(const Line3 &other) const
	{
		return !(start.equals(other.start) && end.equals(other.end)) || (end.equals(other.start) && start.equals(other.end));
	}

	void set_line(const double_t &xa, const double_t &ya, const double_t &za,
				  const double_t &xb, const double_t &yb, const double_t &zb)
	{
		start = Vec3<T>(xa, ya, za);
		end = Vec3<T>(xb, yb, zb);
	}
	void set_line(const Vec3<T> &nstart, const Vec3<T> &nend)
	{
		start = nstart;
		end = nend;
	}
	void set_line(const Line3 &line)
	{
		start = line.start;
		end = line.end;
	}
	double length() const
	{
		return make_vector(start, end).length();
	}
	double length_squared() const
	{
		return make_vector(start, end).length_squared();
	}
	Vec3<T> get_middle() const
	{
		return (start + end) * 0.5;
	}
	Vec3<T> get_vector() const
	{
		return end - start;
	}
	bool is_point_between_start_and_end(const Vec3<T> &point) const
	{
		return point.is_between_points(start, end);
	}
	Vec3<T> get_closest_point(const Vec3<T> &point) const
	{
		Vec3<T> c = point - start;
		Vec3<T> v = end - start;
		const T d = v.length();
		if (is_zero(d))
		{
			return start;
		}
		v /= d;
		const T t = v.dot(c);
		if (t < 0.0)
			return start;
		if (t > d)
			return end;
		v *= t;
		return start + v;
	}
	bool get_intersection_with_sphere(Vec3<T> sorigin,
									  T sradius, T &outdistance) const
	{
		Vec3<T> q = sorigin - start;
		const T c = q.length();
		const T v = q.dot(get_vector().normalise());
		const T d = sradius * sradius - (c * c - v * v);
		if (d < 0.0)
			return false;
		outdistance = v - std::sqrt(d);
		return true;
	}
};

enum EIntersectionRelation3D
{
	ISREL3D_FRONT = 0,
	ISREL3D_BACK,
	ISREL3D_PLANAR,
	ISREL3D_SPANNING,
	ISREL3D_CLIPPED
};

template <typename T = double_t>
class Plane
{
public:
	Plane()
	{
		normal = Vec3<>(0, 1, 0);
		recalculate(Vec3<>(0, 0, 0));
	}
	Plane(const Vec3<> &mpoint, const Vec3<> &nml)
	{
		normal = nml;
		recalculate(mpoint);
	}
	Plane(T px, T py, T pz, T nx, T ny, T nz)
	{
		normal.x = nx;
		normal.y = ny;
		normal.z = nz;
		recalculate(Vec3<>(px, py, pz));
	}
	Plane(const Vec3<> &point1, const Vec3<> &point2,
		  const Vec3<> &point3)
	{
		set_plane(point1, point2, point3);
	}
	bool operator==(const Plane<> &other) const
	{
		return (D == other.D && normal.equals(other.normal));
	}
	bool operator!=(const Plane<> &other) const
	{
		return !(D == other.D && normal.equals(other.normal));
	}
	void set_plane(
		const Vec3<> &point, const Vec3<> &nvector)
	{
		normal = nvector;
		recalculate(point);
	}
	void set_plane(const Vec3<> &nvect, T d)
	{
		normal = nvect;
		D = d;
	}
	void set_plane(
		const Vec3<> &point1, const Vec3<> &point2, const Vec3<> &point3)
	{
		normal = (point2 - point1).cross(point3 - point1);
		normal.normalise();
		recalculate(point1);
	}
	bool get_intersection_with_line(
		const Vec3<> &linePoint,
		const Vec3<> &lineVect,
		Vec3<> &outIntersection) const
	{
		T t2 = normal.dot(lineVect);
		if (t2 == 0)
			return false;
		T t = -(normal.dot(linePoint) + D) / t2;
		outIntersection = linePoint + (lineVect * t);
		return true;
	}
	T get_known_intersection_with_line(
		const Vec3<> &linePoint1, const Vec3<> &linePoint2) const
	{
		Vec3<> vect = linePoint2 - linePoint1;
		const T t2 = normal.dot(vect);
		return -((normal.dot(linePoint1) + D) / t2);
	}
	bool get_intersection_with_limited_line(
		const Vec3<> &linePoint1,
		const Vec3<> &linePoint2,
		Vec3<> &outIntersection) const
	{
		bool tmp1 = get_intersection_with_line(linePoint1,
											   linePoint2 - linePoint1,
											   outIntersection);
		bool tmp2 = outIntersection.is_between_points(linePoint1,
													  linePoint2);
		return (tmp1 && tmp2);
	}
	EIntersectionRelation3D classify_point_relation(
		const Vec3<> &point) const
	{
		const T d = normal.dot(point) + D;
		if (d < -ROUNDING_ERROR)
			return ISREL3D_BACK;
		if (d > ROUNDING_ERROR)
			return ISREL3D_FRONT;
		return ISREL3D_PLANAR;
	}
	void recalculate(const Vec3<> &MPoint)
	{
		D = -MPoint.dot(normal);
	}
	Vec3<> get_member_point() const
	{
		return normal * -D;
	}
	bool exists_intersection(const Plane<> &other) const
	{
		Vec3<> cross = other.normal.cross(normal);
		return static_cast<bool>(cross.length());
	}
	bool get_intersection_with_plane(
		Plane<> &other, Vec3<> &outLinePoint, Vec3<> &outLineVect)
	{
		const T fn00 = normal.length();
		const T fn01 = normal.dot(other.normal);
		const T fn11 = other.normal.length();
		const T det = fn00 * fn11 - fn01 * fn01;
		if (std::fabs(det) < ROUNDING_ERROR)
			return false;
		const T invdet = 1.0 / det;
		const T fc0 = (fn11 * -D + fn01 * other.D) * invdet;
		const T fc1 = (fn00 * -other.D + fn01 * D) * invdet;
		outLineVect = normal.cross(other.normal);
		outLinePoint = normal * fc0 + other.normal * fc1;
		return true;
	}
	bool get_intersection_with_planes(Plane<> &o1, Plane &o2, Vec3<> &outPoint)
	{
		Vec3<> linePoint, lineVect;
		if (get_intersection_with_plane(o1, linePoint, lineVect))
			return o2.get_intersection_with_line(linePoint, lineVect, outPoint);
		return false;
	}
	bool is_front_facing(const Vec3<> &lookDirection) const
	{
		const T d = normal.dot(lookDirection);
		return (d <= 0.0);
	}
	T get_distance_to(const Vec3<> &point) const
	{
		return point.dot(normal) + D;
	}
	Vec3<> normal;
	T D;
};

template <typename T = double_t>
class Matrix
{

#ifndef GLenum
	typedef unsigned int GLenum;
#endif

public:
	union
	{
		T m[16];
		T d[4][4];
		struct
		{
			T x[3];
			T dpx;
			T y[3];
			T dpy;
			T z[3];
			T dpz;
			union
			{
				struct
				{
					T tx, ty, tz, tw;
				};
				T t[3];
			};
		};
	};
	Matrix()
	{
		set_identity();
	}
	Matrix(GLenum pname)
	{
		glGetdoublev(pname, d[0]);
	}
	Matrix(const Matrix &ma)
	{
		set(ma);
	}
	Matrix(const Vec3<> &tl)
	{
		set_identity();
		tx = tl.x;
		ty = tl.y;
		tz = tl.z;
	}
	Matrix(T yrot, T xrot, T zrot)
	{
		// T cx, cy, cz;
		// T sx, sy, sz;
		const T cx = std::cos(xrot);
		const T cy = std::cos(yrot);
		const T cz = std::cos(zrot);
		const T sx = std::sin(xrot);
		const T sy = std::sin(yrot);
		const T sz = std::sin(zrot);

		d[0][0] = cy * cz;
		d[0][1] = cy * sz;
		d[0][2] = -sy;
		d[1][0] = sx * sy * cz - cx * sz;
		d[1][1] = sx * sy * sz + cx * cz;
		d[1][2] = sx * cy;
		d[2][0] = cx * sy * cz + sx * sz;
		d[2][1] = cx * sy * sz - sx * cz;
		d[2][2] = cx * cy;

		dpx = dpy = dpz = 0;
		tx = ty = tz = 0;
		tw = 1.0;
	}
	T &operator[](size_t index)
	{
		return m[index];
	}
	void get_glu_lookat_params(Vec3<> *eye, Vec3<> *center, Vec3<> *up)
	{
		center->x = eye->x + z[0];
		center->y = eye->y + z[1];
		center->z = eye->z + z[2];
		up->x = y[0];
		up->y = y[1];
		up->z = y[2];
	}
	void set_identity()
	{
		static T ident[16] =
			{
				1.0, 0.0, 0.0, 0.0,
				0.0, 1.0, 0.0, 0.0,
				0.0, 0.0, 1.0, 0.0,
				0.0, 0.0, 0.0, 1.0};
		for (auto i = 0; i < 16; i++)
			m[i] = ident[i];
	}
	Matrix &operator+(const Matrix &o)
	{
		for (auto i = 0; i < 16; i++)
			m[i] += o.m[i];
		return *this;
	}
	Matrix &operator-(const Matrix &o)
	{
		for (auto i = 0; i < 16; i++)
			m[i] -= o.m[i];
		return *this;
	}
	Matrix &operator*(const Matrix &om)
	{
		Matrix result = (*this).mult(om);
		for (auto i = 0; i < 16; i++)
			m[i] = result.m[i];
		return *this;
	}
	Matrix &operator*(const T &d)
	{
		for (auto i = 0; i < 16; i++)
			m[i] *= d;
		return *this;
	}

	T &operator()(const size_t row, const size_t col)
	{
		return m[row * 4 + col];
	}
	Matrix &operator=(const Matrix &ma)
	{
		return set(ma);
	}
	const bool operator==(const Matrix &ma) const
	{
		for (auto i = 0; i < 16; i++)
		{
			if (!dbl_equals(m[i], ma.m[i], 0.00001f))
				return false;
		}
		return true;
	}
	void mult(Vec3<> &out, const Vec3<> &v)
	{
		for (auto j = 0; j < 4; j++)
		{
			T sum = 0;
			for (int k = 0; k < 4; k++)
				sum += v.d[k] * d[k][j];
			out.d[j] = sum;
		}
	}
	Matrix &set(const Matrix &ma)
	{
		for (auto i = 0; i < 16; i++)
			m[i] = ma.m[i];
		return *this;
	}
	Matrix &mult(const Matrix &m)
	{
		return set(mult(m, *this));
	}
	Matrix mult(const Matrix &m1, const Matrix &om)
	{
		Matrix result;
		result.m[0] = m1.m[0] * om.m[0] + m1.m[4] * om.m[1] + m1.m[8] * om.m[2] + m1.m[12] * om.m[3];
		result.m[4] = m1.m[0] * om.m[4] + m1.m[4] * om.m[5] + m1.m[8] * om.m[6] + m1.m[12] * om.m[7];
		result.m[8] = m1.m[0] * om.m[8] + m1.m[4] * om.m[9] + m1.m[8] * om.m[10] + m1.m[12] * om.m[11];
		result.m[12] = m1.m[0] * om.m[12] + m1.m[4] * om.m[13] + m1.m[8] * om.m[14] + m1.m[12] * om.m[15];

		result.m[1] = m1.m[1] * om.m[0] + m1.m[5] * om.m[1] + m1.m[9] * om.m[2] + m1.m[13] * om.m[3];
		result.m[5] = m1.m[1] * om.m[4] + m1.m[5] * om.m[5] + m1.m[9] * om.m[6] + m1.m[13] * om.m[7];
		result.m[9] = m1.m[1] * om.m[8] + m1.m[5] * om.m[9] + m1.m[9] * om.m[10] + m1.m[13] * om.m[11];
		result.m[13] = m1.m[1] * om.m[12] + m1.m[5] * om.m[13] + m1.m[9] * om.m[14] + m1.m[13] * om.m[15];

		result.m[2] = m1.m[2] * om.m[0] + m1.m[6] * om.m[1] + m1.m[10] * om.m[2] + m1.m[14] * om.m[3];
		result.m[6] = m1.m[2] * om.m[4] + m1.m[6] * om.m[5] + m1.m[10] * om.m[6] + m1.m[14] * om.m[7];
		result.m[10] = m1.m[2] * om.m[8] + m1.m[6] * om.m[9] + m1.m[10] * om.m[10] + m1.m[14] * om.m[11];
		result.m[14] = m1.m[2] * om.m[12] + m1.m[6] * om.m[13] + m1.m[10] * om.m[14] + m1.m[14] * om.m[15];

		result.m[3] = m1.m[3] * om.m[0] + m1.m[7] * om.m[1] + m1.m[11] * om.m[2] + m1.m[15] * om.m[3];
		result.m[7] = m1.m[3] * om.m[4] + m1.m[7] * om.m[5] + m1.m[11] * om.m[6] + m1.m[15] * om.m[7];
		result.m[11] = m1.m[3] * om.m[8] + m1.m[7] * om.m[9] + m1.m[11] * om.m[10] + m1.m[15] * om.m[11];
		result.m[15] = m1.m[3] * om.m[12] + m1.m[7] * om.m[13] + m1.m[11] * om.m[14] + m1.m[15] * om.m[15];
		return result;
	}
	Matrix &translate(const Vec3<T> &tl)
	{
		tx += tl.x;
		ty += tl.y;
		tz += tl.z;
		return *this;
	}
	Matrix &translate(const Vec2<T> &tl)
	{
		tx += tl.x;
		ty += tl.y;
		return *this;
	}
	Vec3<T> transform(const Vec3<T> &v) const
	{
		Vec3<T> res;
		res.x = v.x * d[0][0] + v.y * d[1][0] + v.z * d[2][0] + t[0];
		res.y = v.x * d[0][1] + v.y * d[1][1] + v.z * d[2][1] + t[1];
		res.z = v.x * d[0][2] + v.y * d[1][2] + v.z * d[2][2] + t[2];
		return res;
	}
	Vec2<T> transform(const Vec2<T> &v) const
	{
		Vec2<T> res;
		res.x = v.x * d[0][0] + v.y * d[1][0] + t[0];
		res.y = v.x * d[0][1] + v.y * d[1][1] + t[1];
		return res;
	}
	Vec3<T> rotate(const Vec3<T> &v) const
	{
		Vec3<T> res;
		res.x = v.x * d[0][0] + v.y * d[1][0] + v.z * d[2][0];
		res.y = v.x * d[0][1] + v.y * d[1][1] + v.z * d[2][1];
		res.z = v.x * d[0][2] + v.y * d[1][2] + v.z * d[2][2];
		return res;
	}
	Vec2<T> rotate(const Vec2<T> &v) const
	{
		Vec2<T> res;
		res.x = v.x * d[0][0] + v.y * d[1][0];
		res.y = v.x * d[0][1] + v.y * d[1][1];
		return res;
	}
	Matrix &invert()
	{
		Vec3<T> pos(-tx, -ty, -tz);
		std::swap(d[0][1], d[1][0]);
		std::swap(d[0][2], d[2][0]);
		std::swap(d[1][2], d[2][1]);
		tx = ty = tz = 0.0;
		pos = transform(pos);
		tx = pos.x;
		ty = pos.y;
		tz = pos.z;
		return *this;
	}
	void sub_matrix(Matrix &mr, Matrix &mb, int i, int j)
	{
		for (auto di = 0; di < 3; di++)
		{
			for (auto dj = 0; dj < 3; dj++)
			{
				const int si = di + ((di >= i) ? 1 : 0);
				const int sj = dj + ((dj >= j) ? 1 : 0);
				mb.m[di * 3 + dj] = mr.m[si * 4 + sj];
			}
		}
	}
	T m3_determinant(Matrix &ma)
	{
		return ma.m[0] * (ma.m[4] * ma.m[8] - ma.m[7] * ma.m[5]) - ma.m[1] * (ma.m[3] * ma.m[8] - ma.m[6] * ma.m[5]) + ma.m[2] * (ma.m[3] * ma.m[7] - ma.m[6] * ma.m[4]);
	}
	const bool m3_inverse(Matrix &mr, Matrix &ma) const
	{
		const T result = m3_determinant(ma);
		if (fabs(result) < ROUNDING_ERROR)
		{
			mr.set_identity();
			return false;
		}
		mr.m[0] = ma.m[4] * ma.m[8] - ma.m[5] * ma.m[7] / result;
		mr.m[1] = -(ma.m[1] * ma.m[8] - ma.m[7] * ma.m[2]) / result;
		mr.m[2] = ma.m[1] * ma.m[5] - ma.m[4] * ma.m[2] / result;
		mr.m[3] = -(ma.m[3] * ma.m[8] - ma.m[5] * ma.m[6]) / result;
		mr.m[4] = ma.m[0] * ma.m[8] - ma.m[6] * ma.m[2] / result;
		mr.m[5] = -(ma.m[0] * ma.m[5] - ma.m[3] * ma.m[2]) / result;
		mr.m[6] = ma.m[3] * ma.m[7] - ma.m[6] * ma.m[4] / result;
		mr.m[7] = -(ma.m[0] * ma.m[7] - ma.m[6] * ma.m[1]) / result;
		mr.m[8] = ma.m[0] * ma.m[4] - ma.m[1] * ma.m[3] / result;
		return true;
	}
	const bool invert_matrix(Matrix &mresult, Matrix &ma) const
	{
		const T mdet = ma.determinant();
		if (std::fabs(mdet) < ROUNDING_ERROR)
		{
			mresult.set_identity();
			return false;
		}
		for (auto i = 0; i < 4; i++)
		{
			for (auto j = 0; j < 4; j++)
			{
				Matrix mtemp;
				int32_t sign = 1 - ((i + j) % 2) * 2;
				sub_matrix(ma, mtemp, i, j);
				mresult.m[i + j * 4] = (m3_determinant(mtemp) * sign) / mdet;
			}
		}
		return true;
	}
	T determinant()
	{
		T det = 0.;
		T result = 0.;
		T i = 1.;
		Matrix ms3;
		for (auto n = 0; n < 4; n++, i *= -1.)
		{
			sub_matrix(*this, ms3, 0, n);
			det = m3_determinant(ms3);
			result += m[n] * det * i;
		}
		return result;
	}
	const Matrix &transpose_matrix()
	{
		Matrix<> result;
		for (auto i = 0; i < 4; i++)
		{
			for (auto j = 0; j < 4; j++)
				result.m[i * 4 + j] = m[j * 4 + i];
		}
		for (auto i = 0; i < 16; i++)
			m[i] = result.m[i];
		return *this;
	}
	const Matrix &axis_rotation_matrix(
		const T angle, const Vec3<> &axis)
	{
		set_identity();
		Matrix<> result;
		const T rcos = std::cos(angle);
		const T rsin = std::sin(angle);
		result.m[0] = rcos + axis.x * axis.x * (1 - rcos);
		result.m[1] = axis.z * rsin + axis.y * axis.x * (1 - rcos);
		result.m[2] = -axis.y * rsin + axis.z * axis.x * (1 - rcos);
		result.m[4] = -axis.z * rsin + axis.x * axis.y * (1 - rcos);
		result.m[5] = rcos + axis.y * axis.y * (1 - rcos);
		result.m[6] = axis.x * rsin + axis.z * axis.y * (1 - rcos);
		result.m[8] = axis.y * rsin + axis.x * axis.z * (1 - rcos);
		result.m[9] = -axis.x * rsin + axis.y * axis.z * (1 - rcos);
		result.m[10] = rcos + axis.z * axis.z * (1 - rcos);
		for (auto i = 0; i < 16; i++)
			m[i] = result.m[i];
		return *this;
	}
	Matrix &scaling_matrix(const Vec3<T> &scale)
	{
		set_identity();
		Matrix<> result;
		m[0] = scale.x;
		m[5] = scale.y;
		m[10] = scale.z;
		for (auto i = 0; i < 16; i++)
			m[i] = result.m[i];
		return *this;
	}
	Matrix &x_axis_rotation_matrix(T angle)
	{
		set_identity();
		m[5] = std::cos(angle);
		m[6] = -std::sin(angle);
		m[9] = std::sin(angle);
		m[10] = std::cos(angle);
		return *this;
	}
	Matrix &y_axis_rotation_matrix(T angle)
	{
		set_identity();
		m[0] = std::cos(angle);
		m[2] = std::sin(angle);
		m[8] = -std::sin(angle);
		m[10] = std::cos(angle);
		return *this;
	}
	Matrix &z_axis_rotation_matrix(T angle)
	{
		set_identity();
		m[0] = std::cos(angle);
		m[1] = -std::sin(angle);
		m[4] = std::sin(angle);
		m[5] = std::cos(angle);
		return *this;
	}
	void apply()
	{
		glMultMatrixd(&m[0]);
	}
	const Matrix &from_euler(T heading, T pitch, T roll)
	{
		set_identity();
		Matrix result;
		const T A = std::cos(pitch);
		const T B = std::sin(pitch);
		const T C = std::cos(heading);
		const T D = std::sin(heading);
		const T E = std::cos(roll);
		const T F = std::sin(roll);
		const T AD = A * D;
		const T BD = B * D;
		result.m[0] = C * E;
		result.m[1] = -C * F;
		result.m[2] = D;
		result.m[4] = BD * E + A * F;
		result.m[5] = -BD * E + A * E;
		result.m[6] = -B * C;
		result.m[8] = -AD * E + B * F;
		result.m[9] = AD * F + B * E;
		result.m[10] = A * C;
		result.m[3] = result.m[7] = result.m[11] =
			result.m[12] = result.m[13] = result.m[14] = 0.;
		result.m[15] = 1.;
		for (auto i = 0; i < 16; i++)
			m[i] = result.m[i];
		return *this;
	}
	void to_euler(T *heading, T *pitch, T *bank)
	{
		T ax = 0.;
		T az = 0.;
		T ay = m[3] = std::asin(m[2]);
		T c = m[2] = std::cos(ay);

		if (fabs(c) > ROUNDING_ERROR)
		{
			T tr_x = m[10] / c;
			T tr_y = -m[6] / c;
			ax = std::atan2(tr_y, tr_x);
			tr_x = m[0] / c;
			tr_y = -m[1] / c;
			az = std::atan2(tr_y, tr_x);
		}
		else
		{
			ax = 0.;
			T tr_x = m[5];
			T tr_y = m[4];
			az = std::atan2(tr_y, tr_x);
		}
		if (ax < 0)
			ax += 2 * PI;
		if (ay < 0)
			ay += 2 * PI;
		if (az < 0)
			az += 2 * PI;

		*heading = ay;
		*pitch = ax;
		*bank = az;
	}
	void transform(Vec3<T> &out, const Vec3<T> &in) const
	{
		out.x = in.x * m[0] + in.y * m[4] + in.z * m[8] + m[12];
		out.y = in.x * m[1] + in.y * m[5] + in.z * m[9] + m[13];
		out.z = in.x * m[2] + in.y * m[6] + in.z * m[10] + m[14];
	}

	const bool is_identity() const
	{
		if (!dbl_equals(m[0], 1.f) ||
			!dbl_equals(m[5], 1.f) ||
			!dbl_equals(m[10], 1.f) ||
			!dbl_equals(m[15], 1.f))
			return false;

		for (auto i = 0; i < 4; ++i)
			for (auto j = 0; j < 4; ++j)
				if ((j != i) && (!is_zero(((Matrix) * this)(i, j))))
					return false;

		return true;
	}
};

template <typename T = double_t>
class Quat
{
public:
	union
	{
		T d[4]{0, 0, 0, 0};
		struct
		{
			T w, x, y, z;
		};
	};
	Quat() : w(1), x(0), y(0), z(0) {}
	Quat(T h, T p, T b)
	{
		set(h, p, b);
	}
	Quat(const Matrix<> &m)
	{
		T tr, s, q[4];
		int32_t i, j, k;
		int32_t nxt[3] = {1, 2, 0};
		tr = m.d[0][0] + m.d[1][1] + m.d[2][2];
		// check the diagonal
		if (tr > 0.0)
		{
			s = std::sqrt(tr + 1.0);
			this->w = s / 2.0;
			s = 0.5 / s;
			this->x = (m.d[1][2] - m.d[2][1]) * s;
			this->y = (m.d[2][0] - m.d[0][2]) * s;
			this->z = (m.d[0][1] - m.d[1][0]) * s;
		}
		else
		{
			// diagonal is negative
			i = 0;
			if (m.d[1][1] > m.d[0][0])
				i = 1;
			if (m.d[2][2] > m.d[i][i])
				i = 2;
			j = nxt[i];
			k = nxt[j];
			s = std::sqrt((m.d[i][i] - (m.d[j][j] + m.d[k][k])) + 1.0);
			q[i] = s * 0.5;
			if (s != 0.0)
				s = 0.5 / s;
			q[3] = (m.d[j][k] - m.d[k][j]) * s;
			q[j] = (m.d[i][j] + m.d[j][i]) * s;
			q[k] = (m.d[i][k] + m.d[k][i]) * s;
			this->x = q[0];
			this->y = q[1];
			this->z = q[2];
			this->w = q[3];
		}
	}
	T &operator[](size_t index)
	{
		return d[index];
	}
	void set(T h, T p, T b)
	{
		const T cr = std::cos(p / 2);
		const T cp = std::cos(h / 2);
		const T cy = std::cos(b / 2);
		const T sr = std::sin(p / 2);
		const T sp = std::sin(h / 2);
		const T sy = std::sin(b / 2);
		const T cpcy = cp * cy;
		const T spsy = sp * sy;
		w = cr * cpcy + sr * spsy;
		x = sr * cpcy - cr * spsy;
		y = cr * sp * cy + sr * cp * sy;
		z = cr * cp * sy - sr * sp * cy;
	}
	void set(const Quat &q)
	{
		w = q.w;
		x = q.x;
		y = q.y;
		z = q.z;
	}
	void set(const Polar<> &p)
	{
		const T sina = std::sin(p.rad / 2);
		const T cosa = std::cos(p.rad / 2);
		const T sinlat = std::sin(p.lat);
		const T coslat = std::cos(p.lat);
		const T sinlon = std::sin(p.lon);
		const T coslon = std::cos(p.lon);
		x = sina * coslat * sinlon;
		y = sina * sinlat;
		z = sina * sinlat * coslon;
		w = cosa;
	}
	Quat mult(const Quat<> &q1, const Quat<> &q2) const
	{
		const T A = (q1.w + q1.x) * (q2.w + q2.x);
		const T B = (q1.z - q1.y) * (q2.y - q2.z);
		const T C = (q1.w - q1.x) * (q2.y + q2.z);
		const T D = (q1.y + q1.z) * (q2.w - q2.x);
		const T E = (q1.x + q1.z) * (q2.x + q2.y);
		const T F = (q1.x - q1.z) * (q2.x - q2.y);
		const T G = (q1.w + q1.y) * (q2.w - q2.z);
		const T H = (q1.w - q1.y) * (q2.w + q2.z);
		Quat res;
		res.w = B + (-E - F + G + H) / 2;
		res.x = A - (E + F + G + H) / 2;
		res.y = C + (E - F + G - H) / 2;
		res.z = D + (E - F - G + H) / 2;
		return res;
	}
	Quat &mult(const Quat<> &q2)
	{
		const T A = (w + x) * (q2.w + q2.x);
		const T B = (z - y) * (q2.y - q2.z);
		const T C = (w - x) * (q2.y + q2.z);
		const T D = (y + z) * (q2.w - q2.x);
		const T E = (x + z) * (q2.x + q2.y);
		const T F = (x - z) * (q2.x - q2.y);
		const T G = (w + y) * (q2.w - q2.z);
		const T H = (w - y) * (q2.w + q2.z);
		w = B + (-E - F + G + H) / 2;
		x = A - (E + F + G + H) / 2;
		y = C + (E - F + G - H) / 2;
		z = D + (E - F - G + H) / 2;
		return *this;
	}
	Matrix<> to_matrix()
	{
		Matrix<> m;
		const T x2 = x + x;
		const T y2 = y + y;
		const T z2 = z + z;
		const T xx = x * x2;
		const T xy = x * y2;
		const T xz = x * z2;
		const T yy = y * y2;
		const T yz = y * z2;
		const T zz = z * z2;
		const T wx = w * x2;
		const T wy = w * y2;
		const T wz = w * z2;
		m.d[0][0] = 1.0 - (yy + zz);
		m.d[1][0] = xy - wz;
		m.d[2][0] = xz + wy;
		m.d[3][0] = 0.0;
		m.d[0][1] = xy + wz;
		m.d[1][1] = 1.0 - (xx + zz);
		m.d[2][1] = yz - wx;
		m.d[3][1] = 0.0;
		m.d[0][2] = xz - wy;
		m.d[1][2] = yz + wx;
		m.d[2][2] = 1.0 - (xx + yy);
		m.d[3][2] = 0.0;
		m.d[0][3] = 0;
		m.d[1][3] = 0;
		m.d[2][3] = 0;
		m.d[3][3] = 1;
		return m;
	}
	void apply()
	{
		to_matrix().apply();
	}
	void pitch(T a)
	{
		Quat<> q;
		q.w = std::cos(a / 2);
		q.x = std::sin(a / 2);
		q.y = 0;
		q.z = 0;
		set(mult(*this, q));
	}
	void yaw(T a)
	{
		Quat<> q;
		q.w = std::cos(a / 2);
		q.x = 0;
		q.y = std::sin(a / 2);
		q.z = 0;
		set(mult(*this, q));
	}
	void roll(T a)
	{
		Quat<> q;
		q.w = std::cos(a / 2);
		q.x = 0;
		q.y = 0;
		q.z = std::sin(a / 2);
		set(mult(*this, q));
	}
	void set(Vec3<> v, T a)
	{
		v.normalise();
		v *= sin(a / 2.0);
		x = v.x;
		y = v.y;
		z = v.z;
		w = cos(a / 2.0);
	}
	void get_headings(T *heading, T *attitude, T *bank) const
	{
		*heading = std::atan2(2 * y * w - 2 * x * z, 1 - 2 * y * y - 2 * z * z);
		*attitude = std::asin(2 * x * y + 2 * z * w);
		*bank = std::atan2(2 * x * w - 2 * y * z, 1 - 2 * x * x - 2 * z * z);
	}
	const T length() const
	{
		return std::sqrt(w * w + x * x + y * y + z * z);
	}
	void normalise()
	{
		const T qlength = length();
		if (qlength == 0)
		{
			w = 1.;
			x = 0.;
			y = 0.;
			z = 0.;
		}
		else
		{
			T inv = 1.0f / qlength;
			x = x * inv;
			y = y * inv;
			z = z * inv;
			w = w * inv;
		}
	}
	void conjugate(const Quat<> &q)
	{
		w = q.w;
		x = -q.x;
		y = -q.y;
		z = -q.z;
	}
	const Quat &operator*(const Quat<> &oq)
	{
		Quat<> result = mult(oq);
		x = result.x;
		y = result.y;
		z = result.z;
		w = result.w;
		return *this;
	}
	Quat<> slerp(Quat<> a, Quat<> b, T t)
	{
		T flip = 1.0;
		T cosine = a.w * b.w + a.x * b.x + a.y * b.y + a.z * b.z;
		if (cosine < 0)
		{
			cosine = -cosine;
			flip = -1;
		}
		if ((1 - cosine) < 0.00001)
		{
			a.w *= (1 - t);
			a.x *= (1 - t);
			a.y *= (1 - t);
			a.z *= (1 - t);
			b.w *= (t * flip);
			b.x *= (t * flip);
			b.y *= (t * flip);
			b.z *= (t * flip);
			a.w = a.w + b.w;
			a.x = a.x + b.x;
			a.y = a.y + b.y;
			a.z = a.z + b.z;
			return a;
		}
		const T theta = std::acos(cosine);
		const T sine = std::sin(theta);
		const T beta = std::sin((1 - t) * theta) / sine;
		const T alpha = std::sin(t * theta) / sine * flip;
		a.w *= beta;
		a.x *= beta;
		a.y *= beta;
		a.z *= beta;
		b.w *= alpha;
		b.x *= alpha;
		b.y *= alpha;
		b.z *= alpha;
		a.w = a.w + b.w;
		a.x = a.x + b.x;
		a.y = a.y + b.y;
		a.z = a.z + b.z;
		return Quat(a);
	}
};

typedef Vec3<float_t> Vec3f;
typedef Vec3<double_t> Vec3d;
typedef Vec2<float_t> Vec2f;
typedef Vec2<double_t> Vec2d;
typedef Vec3<int32_t> Vec3i;
typedef Vec2<int32_t> Vec2i;
typedef Matrix<float_t> Matrixf;
typedef Matrix<double_t> Matrixd;

Vec3d vec3i_to_vec3d(const Vec3i &v);
Vec3d vec3f_to_vec3d(const Vec3f &v);

Vec3i vec3d_to_vec3i(const Vec3d &v);
Vec3i vec3f_to_vec3i(const Vec3f &v);

Vec3f vec3i_to_vec3f(const Vec3i &v);
Vec3f vec3d_to_vec3f(const Vec3d &v);

void quadratic_plot(Vec2f startpos, Vec2f controlpos, Vec2f endpos, std::vector<Vec2f> &plotresult, int32_t numsegments = 32);
void cubic_plot(Vec2f startpos, Vec2f controlpos1, Vec2f endpos, Vec2f controlpos2, std::vector<Vec2f> &plotresult, int32_t numsegments = 32);
int32_t get_line_intersection(float_t p0_x, float_t p0_y,
							  float_t p1_x, float_t p1_y,
							  float_t p2_x, float_t p2_y,
							  float_t p3_x, float_t p3_y,
							  float_t &i_x, float_t &i_y);
Vec2f get_line_intersection(Vec2f v0, Vec2f v1, Vec2f v2, Vec2f v3, bool &success);

inline const Vec3f lerp3df(const Vec3f &v1, const Vec3f &v2, const float_t &a)
{
	Vec3f result;
	result.x = lerp(v1.x, v2.x, a);
	result.y = lerp(v1.y, v2.y, a);
	result.z = lerp(v1.z, v2.z, a);
	return result;
}
inline const Vec3d lerp3dd(const Vec3d &v1, const Vec3d &v2, const double_t &a)
{
	Vec3d result;
	result.x = lerp(v1.x, v2.x, a);
	result.y = lerp(v1.y, v2.y, a);
	result.z = lerp(v1.z, v2.z, a);
	return result;
}
inline const Vec2f lerp2df(const Vec2f &v1, const Vec2f &v2, const float_t &a)
{
	Vec2f result;
	result.x = lerp(v1.x, v2.x, a);
	result.y = lerp(v1.y, v2.y, a);
	return result;
}
inline const Vec2d lerp2dd(const Vec2d &v1, const Vec2d &v2, const double_t &a)
{
	Vec2d result;
	result.x = lerp(v1.x, v2.x, a);
	result.y = lerp(v1.y, v2.y, a);
	return result;
}

// AABBox is lifted from Irrlicht (:

// Copyright (C) 2002-2008 Nikolaus Gebhardt
// This file is part of the "Irrlicht Engine".
// For conditions of distribution and use, see copyright notice in irrlicht.h
/*
This file is borrowed from Irrlicht engine.
modified for Simlib use 2008 Dave Rowbotham.
*/

//! Axis aligned bounding box in 3d dimensional space.
/** Has some useful methods used with occlusion culling or clipping.
 */
template <typename T>
class AABBox
{
public:
	//! Default Constructor.
	AABBox() : min_edge(-1, -1, -1), max_edge(1, 1, 1) {}
	//! Constructor with min edge and max edge.
	AABBox(const Vec3<T> &min, const Vec3<T> &max)
		: min_edge(min), max_edge(max) {}
	//! Constructor with only one point.
	AABBox(const Vec3<T> &init) : min_edge(init), max_edge(init) {}
	//! Constructor with min edge and max edge as single values, not vectors.
	AABBox(T minx, T miny, T minz, T maxx, T maxy, T maxz) : min_edge(minx, miny, minz), max_edge(maxx, maxy, maxz) {}
	// operators
	//! Equality operator
	/** \param other box to compare with.
	\return True if both boxes are equal, else false. */
	bool operator==(const AABBox<T> &other) const
	{
		return (min_edge.equals(other.min_edge) && other.max_edge.equals(max_edge));
	}
	//! Inequality operator
	/** \param other box to compare with.
	\return True if both boxes are different, else false. */
	bool operator!=(const AABBox<T> &other) const
	{
		return !(min_edge.equals(other.min_edge) && other.max_edge.equals(max_edge));
	}
	// functions
	//! Adds a point to the bounding box
	/** The box grows bigger, if point was outside of the box.
	\param p: Point to add into the box. */
	void add_internal_point(const Vec3<T> &p)
	{
		add_internal_point(p.x, p.y, p.z);
	}
	//! Adds another bounding box
	/** The box grows bigger, if the new box was outside of the box.
	\param b: Other bounding box to add into this box. */
	void add_internal_box(const AABBox<T> &b)
	{
		add_internal_point(b.max_edge);
		add_internal_point(b.min_edge);
	}
	//! Resets the bounding box to a one-point box.
	/** \param x X coord of the point.
	\param y Y coord of the point.
	\param z Z coord of the point. */
	void reset(T x, T y, T z)
	{
		max_edge.set(x, y, z);
		min_edge = max_edge;
	}
	//! Resets the bounding box.
	/** \param initValue New box to set this one to. */
	void reset(const AABBox<T> &initValue) { *this = initValue; }
	//! Resets the bounding box to a one-point box.
	/** \param initValue New point. */
	void reset(const Vec3<T> &initValue)
	{
		max_edge = initValue;
		min_edge = initValue;
	}
	//! Adds a point to the bounding box
	/** The box grows bigger, if point is outside of the box.
	\param x X coordinate of the point to add to this box.
	\param y Y coordinate of the point to add to this box.
	\param z Z coordinate of the point to add to this box. */
	void add_internal_point(T x, T y, T z)
	{
		if (x > max_edge.x)
			max_edge.x = x;
		if (y > max_edge.y)
			max_edge.y = y;
		if (z > max_edge.z)
			max_edge.z = z;
		if (x < min_edge.x)
			min_edge.x = x;
		if (y < min_edge.y)
			min_edge.y = y;
		if (z < min_edge.z)
			min_edge.z = z;
	}
	//! Determines if a point is within this box.
	/** \param p: Point to check.
	\return True if the point is within the box and false if not */
	bool is_point_inside(const Vec3<T> &p) const
	{
		return (p.x >= min_edge.x && p.x <= max_edge.x &&
				p.y >= min_edge.y && p.y <= max_edge.y &&
				p.z >= min_edge.z && p.z <= max_edge.z);
	};
	//! Determines if a point is within this box and its borders.
	/** \param p: Point to check.
	\return True if the point is within the box and false if not. */
	bool is_point_total_inside(const Vec3<T> &p) const
	{
		return (p.x > min_edge.x && p.x < max_edge.x &&
				p.y > min_edge.y && p.y < max_edge.y &&
				p.z > min_edge.z && p.z < max_edge.z);
	};
	//! Determines if the box intersects with another box.
	/** \param other: Other box to check a intersection with.
	\return True if there is an intersection with the other box,
	otherwise false. */
	bool intersects_with_box(const AABBox<T> &other) const
	{
		return (min_edge <= other.max_edge && max_edge >= other.min_edge);
	}
	//! Determines if the box intersects with another box in parent space.
	/** \param other: Other box to check a intersection with.
	\return True if there is an intersection with the other box,
	otherwise false. */
	bool intersects_with_box_in_parent_space(
		const AABBox<T> &other) const
	{
		bool result = false;
		if (parent != nullptr && other.parent != nullptr)
		{
			Vec3<double> e_min_edge = parent->position + min_edge;
			Vec3<double> e_max_edge = parent->position + max_edge;
			Vec3<double> o_min_edge = other.min_edge + other.parent->position;
			Vec3<double> o_max_edge = other.max_edge + other.parent->position;
			result = e_min_edge <= o_max_edge && e_max_edge >= o_min_edge;
		}
		return result;
	}
	//! Check if this box is completely inside the 'other' box.
	/** \param other: Other box to check against.
	\return True if this box is completly inside the other box,
	otherwise false. */
	bool is_full_inside(const AABBox<T> &other) const
	{
		return min_edge >= other.min_edge && max_edge <= other.max_edge;
	}
	//! Tests if the box intersects with a line
	/** \param line: Line to test intersection with.
	\return True if there is an intersection , else false. */
	bool intersects_with_line(const Line3<T> &line) const
	{
		return intersects_with_line(line.get_middle(), line.get_vector().normalise(),
									(T)(line.length() * 0.5));
	}
	bool intersects_with_line(const Line3<float_t> &line)
	{
		return intersects_with_line(line.get_middle(), line.get_vector().normalise(),
									(float_t)(line.length() * 0.5));
	}
	bool intersects_with_line(const Vec3f &linemiddle,
							  const Vec3f &linevect,
							  float_t halflength) const
	{
		const Vec3f e = get_extent() * 0.5;
		const Vec3f t = get_center() - linemiddle;
		if ((std::fabs(t.x) > e.x + halflength * std::fabs(linevect.x)) ||
			(std::fabs(t.y) > e.y + halflength * std::fabs(linevect.y)) ||
			(std::fabs(t.z) > e.z + halflength * std::fabs(linevect.z)))
			return false;
		float_t r = e.y * std::fabs(linevect.z) + e.z * std::fabs(linevect.y);
		if (std::fabs(t.y * linevect.z - t.z * linevect.y) > r)
			return false;
		r = e.x * std::fabs(linevect.z) + e.z * std::fabs(linevect.x);
		if (std::fabs(t.z * linevect.x - t.x * linevect.z) > r)
			return false;
		r = e.x * std::fabs(linevect.y) + e.y * std::fabs(linevect.x);
		if (std::fabs(t.x * linevect.y - t.y * linevect.x) > r)
			return false;
		return true;
	}
	bool intersects_with_line(const Line3<double_t> &line)
	{
		return intersects_with_line(line.get_middle(), line.get_vector().normalise(),
									(float_t)(line.length() * 0.5));
	}
	bool intersects_with_line(const Vec3d &linemiddle,
							  const Vec3d &linevect,
							  double_t halflength) const
	{
		const Vec3d e = get_extent() * 0.5;
		const Vec3d t = get_center() - linemiddle;
		if ((std::fabs(t.x) > e.x + halflength * std::fabs(linevect.x)) ||
			(std::fabs(t.y) > e.y + halflength * std::fabs(linevect.y)) ||
			(std::fabs(t.z) > e.z + halflength * std::fabs(linevect.z)))
			return false;
		double_t r = e.y * std::fabs(linevect.z) + e.z * std::fabs(linevect.y);
		if (std::fabs(t.y * linevect.z - t.z * linevect.y) > r)
			return false;
		r = e.x * std::fabs(linevect.z) + e.z * std::fabs(linevect.x);
		if (std::fabs(t.z * linevect.x - t.x * linevect.z) > r)
			return false;
		r = e.x * std::fabs(linevect.y) + e.y * std::fabs(linevect.x);
		if (std::fabs(t.x * linevect.y - t.y * linevect.x) > r)
			return false;
		return true;
	}
	//! Tests if the box intersects with a line
	/** \param linemiddle Center of the line.
	\param linevect Vector of the line.
	\param halflength Half length of the line.
	\return True if there is an intersection, else false. */
	bool intersects_with_line(const Vec3<T> &linemiddle,
							  const Vec3<T> &linevect,
							  T halflength)
	{
		const Vec3<T> e = get_extent() * 0.5;
		const Vec3<T> t = get_center() - linemiddle;
		if ((std::fabs(t.x) > e.x + halflength * std::fabs(linevect.x)) ||
			(std::fabs(t.y) > e.y + halflength * std::fabs(linevect.y)) ||
			(std::fabs(t.z) > e.z + halflength * std::fabs(linevect.z)))
			return false;
		T r = e.y * std::fabs(linevect.z) + e.z * std::fabs(linevect.y);
		if (std::fabs(t.y * linevect.z - t.z * linevect.y) > r)
			return false;
		r = e.x * std::fabs(linevect.z) + e.z * std::fabs(linevect.x);
		if (std::fabs(t.z * linevect.x - t.x * linevect.z) > r)
			return false;
		r = e.x * std::fabs(linevect.y) + e.y * std::fabs(linevect.x);
		if (std::fabs(t.x * linevect.y - t.y * linevect.x) > r)
			return false;
		return true;
	}
	//! Classifies a relation with a plane.
	/** \param plane Plane to classify relation to.
	\return Returns ISREL3D_FRONT if the box is in front of the plane,
	ISREL3D_BACK if the box is behind the plane, and
	ISREL3D_CLIPPED if it is on both sides of the plane. */
	EIntersectionRelation3D classify_plane_relation(
		const Plane<T> &plane) const
	{
		Vec3<T> nearPoint(max_edge);
		Vec3<T> farPoint(min_edge);
		if (plane.normal.x > 0.0)
		{
			nearPoint.x = min_edge.x;
			farPoint.x = max_edge.x;
		}
		if (plane.normal.z > 0.0)
		{
			nearPoint.y = min_edge.y;
			farPoint.y = max_edge.y;
		}
		if (plane.normal.z > 0.0)
		{
			nearPoint.z = min_edge.z;
			farPoint.z = max_edge.z;
		}
		if (plane.normal.dot(nearPoint) + plane.D > 0.0)
			return ISREL3D_FRONT;
		if (plane.normal.dot(farPoint) + plane.D > 0.0)
			return ISREL3D_CLIPPED;
		return ISREL3D_BACK;
	}
	//! Get center of the bounding box
	/** \return Center of the bounding box. */
	Vec3<T> get_center() const
	{
		// return (min_edge + max_edge) / 2.0;
		Vec3<T> v;
		v.x = (min_edge.x + max_edge.x) / (T)2.0;
		v.y = (min_edge.y + max_edge.y) / (T)2.0;
		v.z = (min_edge.z + max_edge.z) / (T)2.0;
		return v;
	}
	//! Get extent of the box
	/** \return Extent of the bounding box. */
	Vec3<T> get_extent() const
	{
		// return max_edge - min_edge;
		Vec3<T> v;
		v.x = max_edge.x - min_edge.x;
		v.y = max_edge.y - min_edge.y;
		v.z = max_edge.z - min_edge.z;
		return v;
	}
	//! Stores all 8 edges of the box into an array
	/** \param edges: Pointer to array of 8 edges. */
	void get_edges(Vec3<T> *edges) const
	{
		const Vec3<T> middle = get_center();
		Vec3<T> diag; // = middle - max_edge;
		diag.x = middle.x - max_edge.x;
		diag.y = middle.y - max_edge.y;
		diag.z = middle.z - max_edge.z;
		/*
		Edges are stored in this way:
		Hey, am I an ascii artist, or what? :) niko.
			  /3--------/7
			 /  |      / |
			/   |     /  |
			1---------5  |
			|   2- - -| -6
			|  /      |  /
			|/        | /
			0---------4/
		*/
		edges[0].set(middle.x + diag.x, middle.y + diag.y, middle.z + diag.z);
		edges[1].set(middle.x + diag.x, middle.y - diag.y, middle.z + diag.z);
		edges[2].set(middle.x + diag.x, middle.y + diag.y, middle.z - diag.z);
		edges[3].set(middle.x + diag.x, middle.y - diag.y, middle.z - diag.z);
		edges[4].set(middle.x - diag.x, middle.y + diag.y, middle.z + diag.z);
		edges[5].set(middle.x - diag.x, middle.y - diag.y, middle.z + diag.z);
		edges[6].set(middle.x - diag.x, middle.y + diag.y, middle.z - diag.z);
		edges[7].set(middle.x - diag.x, middle.y - diag.y, middle.z - diag.z);
	}
	//! Check if the box is empty.
	/** This means that there is no space between the min and max
	edge.
	\return True if box is empty, else false. */
	bool is_empty() const { return min_edge.equals(max_edge); }
	//! Repairs the box.
	/** Necessary if for example MinEdge and MaxEdge are swapped. */
	void repair()
	{
		T t;
		if (min_edge.x > max_edge.x)
		{
			t = min_edge.x;
			min_edge.x = max_edge.x;
			max_edge.x = t;
		}
		if (min_edge.y > max_edge.y)
		{
			t = min_edge.y;
			min_edge.y = max_edge.y;
			max_edge.y = t;
		}
		if (min_edge.z > max_edge.z)
		{
			t = min_edge.z;
			min_edge.z = max_edge.z;
			max_edge.z = t;
		}
	}
	//! Calculates a new interpolated bounding box.
	/** \param other: other box to interpolate between
	\param d: value between 0.0f and 1.0f.
	\return Interpolated box. */
	AABBox<T> get_interpolated(const AABBox<T> &other, float_t d) const
	{
		float_t inv = 1.0f - d;
		return AABBox((other.min_edge * inv) + (min_edge * d),
					  (other.max_edge * inv) + (max_edge * d));
	}
	void set_parent(T *p) { parent = p; }
	//! The near edge
	Vec3<T> min_edge;
	//! The far edge
	Vec3<T> max_edge;
	AABBox<T> *parent;
};

typedef AABBox<float_t> AABBoxf;
typedef AABBox<double_t> AABBoxd;
typedef AABBox<int32_t> AABBoxi;

bool box_sphere_intersect(AABBoxd &abox, Vec3d position, double_t radius);
bool box_sphere_intersect(AABBoxf &abox, Vec3f position, float_t radius);
bool box_box_intersect(AABBoxd &abox, AABBoxd &bbox);
bool box_box_intersect(AABBoxf &abox, AABBoxf &bbox);
