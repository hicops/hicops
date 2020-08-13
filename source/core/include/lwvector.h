/*
 * Copyright (C) 2020  Muhammad Haseeb, and Fahad Saeed
 * Florida International University, Miami, FL
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 *
 */

#ifndef LWVECTOR_H_
#define LWVECTOR_H_

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <string>
#include <cstring>
#include "config.h"

using namespace std;

using int_t = int;
using VOID = void;
using double_t = double;

template<class T>
class lwvector
{
private:
    T *arr;

    int_t head;
    int_t tail;
    int_t cap;
    int_t sze;

public:
    lwvector()
    {
        const int_t DSIZE = 2 + (100 * 10);
        arr = new T [DSIZE];
        cap = DSIZE;
        head = 0;
        tail = 0;
        sze = 0;
    }

    lwvector(int_t dsiz)
    {
        cap = dsiz;
        arr = new T [cap];
        head = 0;
        tail = 0;
        sze = 0;
        std::memset(arr + head, 0x0, sizeof(T) * cap);
    }

    lwvector(const lwvector &rhs)
    {
        cap = rhs.cap;
        head = rhs.head;
        tail = rhs.tail;
        sze = rhs.sze;

        arr = new T[cap];
        std::memcpy(this->begin(), &rhs.arr[rhs.head], sizeof(T) * sze);
    }

    lwvector(T *p1, T *p2)
    {
        cap = 2 * std::distance(p1, p2);
        arr = new T[cap];
        head = 0;
        tail = 0;
        sze = 0;

        for (auto it = p1; it < p2; it++)
        {
            arr[tail] = *it;
            tail++;
            sze++;
        }
    }

    lwvector(int_t dsiz, T val)
    {
        arr = new T [dsiz];
        cap = dsiz;
        head = 0;
        tail = dsiz;
        sze = dsiz;

        for (int_t l = 0; l < dsiz; l++)
        {
            arr[l] = val;
        }
    }

    VOID setmem(const T &val)
    {
        memset(arr+head, val, sizeof(T) * sze);
    }

    VOID AddRange(int_t stt, int_t enda)
    {
        int_t ss = enda - stt + 1;

        if (sze + ss <= cap)
        {
            for (int_t kk = stt; kk <= enda; kk++)
            {
                arr[tail] = kk;
                tail++;
                sze++;
            }
        }
    }

    VOID MakeRange(int_t stt, int_t enda)
    {
        int_t ss = enda - stt + 1;

        if (ss <= cap)
        {
            /* Make a range */
            head = 0;
            tail = 0;
            sze = ss;

            for (auto kk = stt; kk <= enda; kk++)
            {
                arr[tail] = kk;
                tail++;
            }
        }
    }

    virtual ~lwvector()
    {
        head = tail = 0;
        sze = cap = 0;

        if (arr != NULL)
        {
            delete[] arr;
            arr = NULL;
        }
    }

    T *begin()
    {
        return arr + head;
    }

    T *end()
    {
        return arr + tail;
    }

    VOID Assign(T *p1, T *p2)
    {
        head = 0;
        tail = 0;
        sze = 0;

        for (auto it = p1; it < p2; it++)
        {
            arr[tail] = *it;
            tail++;
            sze++;
        }
    }

    VOID clip(int_t begin, int_t end)
    {
        tail = head + end + 1;
        head += begin;
        sze = end - begin + 1;
    }

    int_t capacity()
    {
        return cap;
    }

    int_t Size()
    {
        return sze;
    }

    VOID clear()
    {
        head = 0;
        tail = 0;
        sze = cap = 0;
        if (arr != nullptr)
        {
            delete[] arr;
            arr = nullptr;
        }
    }

    VOID print()
    {
        cout << "PRINT: " << endl;

        for (auto i = head; i < tail; i++)
        {
            cout << arr[i] << ", ";
        }

        cout << endl;
    }

    VOID print(int_t p1, int_t p2)
    {
        if (p2 > p1)
        {
            cout << "PRINT: head@:" << head << endl;

            for (int_t ikr = head + p1; ikr < head + p2; ikr++)
            {
                cout << ikr << ", ";
            }
            cout << endl;

            for (int_t ikr = head + p1; ikr < head + p2; ikr++)
            {
                cout << arr[head + ikr] << ", ";
            }

            cout << endl;
        }
    }

    VOID Erase()
    {
        head = 0;
        tail = 0;
        memset(arr, 0x0, sizeof(T) * sze);
        sze = 0;
    }

    VOID multiply(const double_t val)
    {
        for (auto ii = head; ii < tail; ii++)
        {
            arr[ii] *= val;
        }
    }

    VOID multiply(const lwvector&rhs)
    {
        if (rhs.sze != this->sze)
        {
            cout << "FATAL: sizes don't match" << endl;
        }
        else
        {
            for (auto ii = 0; ii < sze; ii++)
            {
                arr[head + ii] *= rhs.arr[rhs.head + ii];
            }
        }
    }

    T dot(const lwvector &rhs)
    {
        T vv = 0;

        if (rhs.sze != this->sze)
        {
            cout << "FATAL: sizes don't match" << endl;
        }
        else
        {
            for (auto ii = 0; ii < sze; ii++)
            {
                vv += arr[head + ii] * rhs.arr[rhs.head + ii];
            }
        }

        return vv;
    }

    T sum(T init = 0)
    {
        T vv = init;

        for (auto i = head; i < tail; i++)
        {
            vv += arr[i];
        }

        return vv;
    }

    VOID log10()
    {
        for (auto i = head; i < tail; i++)
        {
            if (arr[i] > 0)
            {
                arr[i] = std::log10(arr[i]);
            }
            else
            {
                cout << "FATAL: log of negative or zero" << endl;
                arr[i] = -INFINITY;
            }
        }

    }

    VOID divide(const double_t ii)
    {
        if (ii == 0)
        {
            cout << "FATAL: Divide by zero not allowed" << endl;
            return;
        }

        for (auto i = head; i < tail; i++)
        {
            arr[i] /= ii;
        }
    }

    VOID add(const lwvector&rhs)
    {
        if (sze >= rhs.sze)
        {
            for (auto ii = 0; ii < rhs.sze; ii++)
            {
                arr[head + ii] += rhs.arr[rhs.head + ii];
            }
        }
        else
        {
            if (head + rhs.sze <= cap)
            {
                auto ii = 0;

                for (auto ii = 0; ii < sze; ii++)
                {
                    arr[head + ii] += rhs.arr[rhs.head + ii];
                }
                for (; ii < rhs.sze; ii++)
                {
                    arr[head + ii] = rhs.arr[rhs.head + ii];
                }

                sze = rhs.sze;
                tail = sze;
            }
            else
            {

            }
        }
    }

    T *data()
    {
        return arr + head;
    }

    VOID AddData(T *p1, T *p2)
    {
        for (auto it = p1; it < p2; it++)
        {
            arr[tail] = (T)*it;
            tail++;
            sze++;
        }
    }

    VOID add(const double_t val)
    {
        for (auto ii = head; ii < tail; ii++)
        {
            arr[ii] += val;
        }
    }

    VOID negative()
    {
        for (auto ii = head; ii < tail; ii++)
        {
            arr[ii] = -arr[ii];
        }
    }

    VOID OneMinusX()
    {
        for (auto ii = head; ii < tail; ii++)
        {
            arr[ii] = 1 - arr[ii];
        }
    }

    T& operator[](int_t index)
    {
        auto actual = head;

        /*if (index > 0)
        {*/
            actual += index;

            if (actual >= tail)
            {
                actual -= tail;
            }
        /*}

        if (index < 0)
        {
            actual = tail - index;

            if (actual < head)
            {
                actual += head;
            }
        }
*/
        return arr[actual];
    }

};

#endif /* LWVECTOR_H_ */
