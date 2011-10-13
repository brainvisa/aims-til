
#include "til/miscTools.h"

#include <cstdarg>		// va_start


// Namespace 

namespace til {

TIL_API ImageAxis operator++(ImageAxis &axis)
{
	int temp = static_cast<int>(axis);
	++temp;
	axis = static_cast<ImageAxis>(temp);
	return axis;
}

TIL_API char * mySprintf(const char * format, ... )
{
	static char *s = new char[1024];
	va_list vl;
	va_start(vl, format);
	vsprintf(s, format, vl);
	return s;
}
} // namespace




