
#pragma once

#include <new>   // For the new operator.

template <class T>
class NonThrowingVector
{
public:
    NonThrowingVector(int sz) {
        if (sz)
            ptr = new(std::nothrow) T[sz];
        else
            ptr = 0;
    }

    ~NonThrowingVector() {
        dispose();
    }

    operator bool() const {
        return ptr ? true : false;
    }

    // implicit conversion to the underlying pointer
    operator T*() {
        return ptr;
    }

    void dispose() {
        if (ptr)
            delete[] ptr;
        ptr = 0;
    }

    T * data() {
        return ptr;
    }

    T & operator[](int idx) {
        return ptr[idx];
    }

    T operator[](int idx) const {
        return ptr[idx];
    }

private:
    T* ptr;
};

