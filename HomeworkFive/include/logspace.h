#ifndef WILL_LOGSPACE_H
#define WILL_LOGSPACE_H

#include <algorithm>
#include <cassert>

namespace will {

template<typename T>
class Logspace
{
    // https://stackoverflow.com/questions/21429294/is-there-something-like-numpy-logspace-in-c
public:
    Logspace(T first, T base) : current(first), base(base)
    {
        assert(base > 0);
        assert(current > 0);
    }

    T operator()()
    {
        T next = current;
        current *= base;
        return next;
    }

private:
    T current, base;
};

} // namespace will

#endif // WILL_LOGSPACE_H
