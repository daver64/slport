#include "template_3dmath.h"
#include <cstdio>
Vec3f p1{1,4,9};
Vec3f p2{2,8,18};

void print_vector(Vec3f& v)
{
	printf("x=%f,y=%f,z=%f\n",v.x,v.y,v.z);
}
int main()
{
	Vec3f p3=p1+p2;
	print_vector(p1);
	print_vector(p2);
	print_vector(p3);
	return 42;
}
